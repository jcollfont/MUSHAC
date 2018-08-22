### IMPORTS
# standard python
import numpy as np 
import nrrd
import os
import shutil
from subprocess import call
import sys
import argparse

# DIPY
from dipy.core.gradients import gradient_table

# local
from MUSHACreconstruction import MUSHACreconstruction
sys.path.insert(0, '/home/ch199899/Documents/Research/DWI/MFMestimation/python/')
from loadDWIdata import saveNRRDwithHeader, getListofDWIFiles, loadDWIdata

# GLOBAL
refHeader = '/fileserver/projects6/jcollfont/refHeader.nhdr'


#%% run Benoit's denoising algo
def runDenoising( inputData, outputName ):

    if not os.path.exists( outputName ): 
        functionDenoising = '/home/ch137122/Software/mrtrix/mrtrix3/release/bin/dwidenoise'
        functionMRconvert = '/home/ch137122/Software/mrtrix/mrtrix3/release/bin/mrconvert'

        tmpdir = '/tmp/tmp%d'  %(np.random.randint(1e6))
        os.makedirs(tmpdir)

        print 'Convert to FSL'
        out = call([ 'crlDWIConvertNHDRForFSL', '-i', inputData, '--data', tmpdir + '/data.nii' ,\
                            '--bvals',  tmpdir + '/bvals', '--bvecs', tmpdir + '/bvecs'])
        print 'Convert to MIF'
        out = call([ functionMRconvert, '-fslgrad', tmpdir + '/bvecs', tmpdir + '/bvals', tmpdir + '/data.nii',  tmpdir + '/data.mif' ])

        print 'run denoise'
        out = call([ functionDenoising, tmpdir + '/data.mif',  tmpdir + '/data-denoised.nii' ])

        print 'convert back to NHDR'
        out = call([ 'crlDWIConvertFSLToNHDR', '-i', tmpdir + '/data-denoised.nii', '-o', outputName ])

        shutil.rmtree(tmpdir)
    else:
        print '%s already exists. NOt computing denoising algorithm.' %(outputName)

#%% run DIAMOND on benchmark data
def runDIAMOND( dataNHDR, inputFolder, mask, outputName, targetDWI=None, numThreads=50):


    # define which DIAMOND code to call
    functionName = '/home/ch137122/bin/crlDCIEstimateJaume'
    
    # run DIAMOND
    if not os.path.exists( outputName[:-5] + '_predicted.nhdr'):
        print 'Computing ' + outputName[:-5] + '_t0.nrrd'
        try:
            out = call([functionName, '-i',  inputFolder + dataNHDR, '-m', mask,'-o', outputName,  \
                        '--residuals' , '-n 3', '-p ' + str(numThreads) , '--automose aicu',\
                        '--fascicle diamondNCcyl', '--fractions_sumto1 1', '--estimateDisoIfNoFascicle 1',\
                        '--predictedsignalscheme',targetNHDR,'--predictedsignal',outputName[:-5]+'_predicted.nhdr'] )
        except :
            print 'Could not run DIAMOND on ' + dataNHDR
    else:
        print outputName[:-5] + '_t0.nrrd'  + ' already computed'

    return outputName[:-5]+'_predicted.nhdr'


def readInBvecsAndBvals( inputFolder, bvecsFile='dwi.bvec', bvalsFile='dwi.bval' ):

    # b values
    fbvals = open( inputFolder + bvalsFile )
    lines = fbvals.readlines()
    fbvals.close()
    splitLine = lines[0].split(' ')
    val = [ float(vv) for vv in splitLine if not vv=='' ]
    bvals = np.array( val )

    # b vectors
    fbvals = open( inputFolder + bvecsFile )
    lines = fbvals.readlines()
    fbvals.close()
    bvecs = []
    for ll in lines:
        splitLine = ll.split(' ')
        vecs = [ float(vv) for vv in splitLine if (not vv=='')&(not vv=='\n') ]

        bvecs.append( np.array( vecs ) )

    # return bvecs, bvals
    return bvals, np.array(bvecs).T


#%%
#
#
#
#
def correctDWIbetweenOrigins( origNHDR, targetNHDR, newNHDR):
    
    # read original NHDR
    coordsOrig, gradOrig = extractOriginMatrix(origNHDR)
    u,s,v = np.linalg.svd(coordsOrig)
    UVOrig = u.dot(v)

    # read target NHDR
    coordsTarget, gradTarget = extractOriginMatrix(targetNHDR)
    u,s,v = np.linalg.svd(coordsTarget)
    UVTarget = u.dot(v)

    # move gradients from origNHDR to coordinates in targetNHDR
    newGrads = UVTarget.dot( np.linalg.inv(UVOrig) ).dot( gradOrig.T ).T
    newGradsList = newGrads.tolist()

    # write gradients in new NHDR file
    fo = open(origNHDR)
    lines = fo.readlines()
    fo.close()

    fo = open(newNHDR, 'w+')
    for ll in lines:
        if ll.find('DWMRI_gradient')>-1:
            grad = newGradsList.pop(0)
            fo.writelines(ll.split(':=')[0] + ':= %0.6f %0.6f %0.6f\n' %(grad[0], grad[1], grad[2]) )
        elif ll.find('space directions:')>-1:
            fo.writelines(ll.split(':')[0] + ': (%0.6f,%0.6f,%0.6f) (%0.6f,%0.6f,%0.6f) (%0.6f,%0.6f,%0.6f) none\n' \
                                    %(coordsTarget[0,0],coordsTarget[1,0],coordsTarget[2,0],\
                                    coordsTarget[0,1],coordsTarget[1,1],coordsTarget[2,1],\
                                    coordsTarget[0,2],coordsTarget[1,2],coordsTarget[2,2]) )
        else:
            fo.writelines(ll)
    fo.close()

    nhdrFiles = getListofDWIFiles(origNHDR)
    for ff in nhdrFiles:
        shutil.copy( os.path.dirname(origNHDR) + '/' + ff, os.path.dirname(newNHDR) + '/' + ff )

    return newGrads

def extractOriginMatrix( NHDR ):

    fo = open(NHDR)
    lines = fo.readlines()
    fo.close()

    coords = range(3)
    grads = []
    for ll in lines:
        if ll.find('space directions:') > -1:
            vec = ll.split(':')[-1].split('(')
            for d1 in range(3):
                vstr = vec[d1+1].split(')')[0].split(',')
                coords[d1] = np.array([ float(vstr[0]),float(vstr[1]),float(vstr[2]) ] )

        if ll.find('DWMRI_gradient')>-1:
            vstr = ll.split(':=')[-1].split(' ')
            grads.append( np.array( [ float(vstr[1]),float(vstr[2]),float(vstr[3]) ] ) )

    return np.array(coords).T, np.array(grads)

#%%
#
#
#
#
def regenerateDWI( baseName, targetName, outputName):

    call(['/home/ch137122/bin/crlDCIRegenerateDWI', '-i', baseName + '_t0.nrrd', '-i', baseName + '_t1.nrrd', '-i', baseName + '_t2.nrrd' \
                        , '--fractions', baseName + '_fractions.nrrd'\
                        ,  '--kappa', baseName + '_kappa.nrrd'\
                        , '--refdwi', targetName \
                        , '-o', outputName ])
                        

# ------------------ MAIN ----------------------- #

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Benchmark code for DWI data')
    parser.add_argument('-d', '--dir', required=True,
                        help='path to subject folders')
    parser.add_argument('-s', '--subj', default='all',
                        help='subject ID')
    parser.add_argument('-seq', '--seq', default='all',
                        help='sequence ID')
    parser.add_argument('-seq_ref', '--seq_ref', default='prisma',
                        help='sequence for refrence ID')
    parser.add_argument('-res_ref', '--res_ref', default='st',
                        help='resolution for refrence ID')
    parser.add_argument('-p', '--threads', default='50',
                        help='number of threads')
    parser.add_argument('--report', '--report', default=None,
                        help='File where to write the report')
    args = parser.parse_args()

    if args.subj == 'all': 
        subjects = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O']
    else:
        subjects = [args.subj]

    if args.seq == 'all':
        sequences = ['prisma','connectom']
    else:
        sequences = [args.seq]


    # for all subjects
    for subj in subjects:

        print 'Subject ' + subj 

        # input data names and paths
        refInputFolder = args.dir + subj + '/' + args.seq_ref + '/' + args.res_ref + '/' 
        dataNHDR = subj + '_' + args.seq_ref + '_' + args.res_ref + '_dwi.nhdr'
        dataDenoisedNHDR = subj + '_' + args.seq_ref + '_' + args.res_ref + '_denoised_dwi.nhdr'
        mask = refInputFolder + 'mask.nrrd'


        # ------------------ DENOISE DATA ----------------------- #
        print '------------------------  RUNNING DENOISING  ------------------------'
        runDenoising( refInputFolder + 'dwi/' + dataNHDR,  refInputFolder + 'dwi/' + dataDenoisedNHDR )

        # ------------------ UPSAMPLE  ----------------------- # 
        print '------------------------  RUNNING UPSAMPLE  ------------------------'
        upsampledNHDR = refInputFolder + 'dwi_iso1mm/' +  dataDenoisedNHDR[:-5] + '_iso1mm.nhdr'
        if not os.path.exists(os.path.dirname(upsampledNHDR)):
                os.mkdir(os.path.dirname(upsampledNHDR))
        if not os.path.exists(upsampledNHDR):
            call(['crlResampler2', '--voxelsize', '1.0,1.0,1.0', '-i', refInputFolder + 'dwi/' + dataDenoisedNHDR, \
                                    '-o' , upsampledNHDR, '--interp sinc', '-p', args.threads])

        # ------------------ UPSAMPLE MASK  ----------------------- # 
        upsampledMask = refInputFolder + 'mask_iso1mm.nrrd'
        if not os.path.exists(upsampledMask):
            call(['crlResampler2', '-g', upsampledNHDR, '-i', mask, \
                                '-o' , upsampledMask, '--interp nearest', '-p', args.threads])

        # ------------------ DIAMOND model  ----------------------- # 
        # print 'Preparing DIAMOND model for extrapolation'
        # diamondModel = MUSHACreconstruction( refInputFolder,'DIAMOND/' ,\
        #                                 'dwi/' + dataNHDR, maskPath='mask.nrrd')
        

        # ------------------ COMPARISON WITH OTHER MODELS ----------------------- #
        for seq in sequences:
            
            for res in ['sa','st']:

                # # set paths
                inputFolder = args.dir + subj + '/' + seq + '/' + res + '/' 
                outputFolder = inputFolder +  'fullDIAMOND_denoised_iso1mm/'
                targetNHDR = inputFolder + 'dwi/' + subj + '_' + seq + '_' + res + '_dwi.nhdr'

                recNHDR = 'diamondREC_' + subj + '_' + args.seq_ref + '_' +  args.res_ref +  '_denoised_iso1mm_2_' + seq + '_' + res + '.nrrd'

                print 'Mapping: ' + subj + '_' + args.seq_ref + '_' +  args.res_ref +  '_denoised_iso1mm'  + \
                        ' to: ' + subj + '_' + seq + '_' + res

                try:
                    os.makedirs( outputFolder )
                except:
                    print '\t-Folder ' + outputFolder


                # ------------------ regenerate DWI AND DIAMOND model ----------------------- # 
                print '------------------------  RUNNING DIAMOND  ------------------------'
                diamondOutput = outputFolder + recNHDR

                # register to target file 
                # call(['crlResampler2', '-g', upsampledNHDR, '-i', targetNHDR , \
                #                         '-o', outputFolder + 'tmp/' + os.path.basename(targetNHDR)[:-5] + '_registered.nhdr',\
                #                         '--tinv','-t',outputFolder + 'tmp/transform' ])
                newTarget = outputFolder + 'tmp/nweTarget.nhdr'
                if not os.path.exists(outputFolder + 'tmp/'):
                    os.makedirs(outputFolder + 'tmp/')
                correctDWIbetweenOrigins( targetNHDR, upsampledNHDR, newTarget)


                # run diamond and apply reconstruction
                predictedNHDR = runDIAMOND( os.path.basename(upsampledNHDR) , \
                                        os.path.dirname(upsampledNHDR)+ '/', \
                                        upsampledMask, \
                                        diamondOutput, \
                                        newTarget, \
                                        numThreads=args.threads )


                print '------------------------  RUNNING RESAMPLER  ------------------------'
                # register to target file 
                # call(['crlBlockMatchingRegistration', '-r', targetNHDR, '-f',  predictedNHDR, \
                #                                     '-o',  predictedNHDR[:-5] + '_registered.nhdr' , \
                #                                     '-i', outputFolder + 'tmp/transform.tfm', \
                #                                     '-N', '-n 0', '-e 0', '-s 0' ,'-p 1', '-t affine'])
                call(['crlResampler2', '-g', targetNHDR, '-i', predictedNHDR , '--interp sinc',\
                                        '-o', predictedNHDR[:-5] + '_registered.nhdr'])

                ##### TODO SANITY CHECK, the registere NHDR should have the same gradients as the target NHDR





                # # evaluate results
                # bvals, bvecs = readInBvecsAndBvals( inputFolder )                                           # define gradients from target file

                # # # ------------------ get FA and MD from DTI model  ----------------------- # 
                # print '\t-computing MD and FA from DTI model'
                # # if not os.path.exists( outputFolder + recNHDR + 'tensStick0.nrrd' ):
                # diamondModel.computeParamsFromSingleTensorFromDWI( \
                #                         bvecs=bvecs, bvals=bvals, recSignal=None, \
                #                         recNHDR= outputFolder + recNHDR + '.nhdr', \
                #                         outputName= outputFolder + recNHDR + '_DTI', \
                #                         anatMask= refInputFolder + 'mask'+nhdrResol+'.nrrd')
