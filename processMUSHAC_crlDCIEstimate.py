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
        functionGibbs = '/home/ch199899/links/DWI/external/mrtrix3/bin/mrdegibbs'
        functionMRconvert = '/home/ch137122/Software/mrtrix/mrtrix3/release/bin/mrconvert'

        tmpdir = '/tmp/tmp%d'  %(np.random.randint(1e6))
        os.makedirs(tmpdir)

        print 'Convert to FSL'
        out = call([ 'crlDWIConvertNHDRForFSL', '-i', inputData, '--data', tmpdir + '/data.nii' ,\
                            '--bvals',  tmpdir + '/bvals', '--bvecs', tmpdir + '/bvecs'])
        print 'Convert to MIF'
        out = call([ functionMRconvert, '-fslgrad', tmpdir + '/bvecs', tmpdir + '/bvals', tmpdir + '/data.nii',  tmpdir + '/data.mif' ])

        print 'run denoise'
        out = call([ functionDenoising, tmpdir + '/data.mif',  tmpdir + '/data-gibbs.mif'])

        print 'run mrdegibbs'
        out = call([ functionGibbs,  tmpdir + '/data-gibbs.mif', tmpdir + '/data-denoised.nii'  ])

        print 'convert back to NHDR'
        out = call([ 'crlDWIConvertFSLToNHDR', '-i', tmpdir + '/data-denoised.nii', '-o', outputName ])

        shutil.rmtree(tmpdir)
    else:
        print '%s already exists. NOt computing denoising algorithm.' %(outputName)

#%% run DIAMOND on benchmark data
def runDIAMOND( dataNHDR, inputFolder, mask, outputName, targetNHDR=None, numThreads=50):


    # define which DIAMOND code to call
    functionName = '/home/ch137122/bin/crlDCIEstimateJaume'
    
    # run DIAMOND
    # if not os.path.exists( outputName[:-5] + '_predicted.nhdr'):
    print 'Computing ' + outputName[:-5] + '_t0.nrrd'
    try:
        out = call([functionName, '-i',  inputFolder + dataNHDR, '-m', mask,'-o', outputName,  \
                    '--residuals' , '-n 3', '-p ' + str(numThreads) , '--automose aicu',\
                    '--fascicle diamondcyl', '--fractions_sumto1 1', '--estimateDisoIfNoFascicle 1',\
                    '--predictedsignalscheme',targetNHDR,'--predictedsignal',outputName[:-5]+'_predicted.nhdr' \
                    ])#,'--bbox 0,0,80,229,229,3'] )
    except :
        print 'Could not run DIAMOND on ' + dataNHDR
            # else:
        #     print outputName[:-5] + '_t0.nrrd'  + ' already computed'
    # else:
    #     print 'DIAMOND already computed. Skipping'

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
        subjects = args.subj.split(',')

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
        if not os.path.exists(refInputFolder + 'dwi_GIBBS/'):
                os.mkdir(refInputFolder + 'dwi_GIBBS/')
        print '------------------------  RUNNING DENOISING  ------------------------'
        runDenoising( refInputFolder + 'dwi/' + dataNHDR,  refInputFolder + 'dwi_GIBBS/' + dataDenoisedNHDR )

        # ------------------ UPSAMPLE  ----------------------- # 
        print '------------------------  RUNNING UPSAMPLE  ------------------------'
        upsampledNHDR = refInputFolder + 'dwi_GIBBS_iso1mm/' +  dataDenoisedNHDR[:-5] + '_iso1mm.nhdr'
        if not os.path.exists(os.path.dirname(upsampledNHDR)):
                os.mkdir(os.path.dirname(upsampledNHDR))
        if not os.path.exists(upsampledNHDR):
            call(['crlResampler2', '--voxelsize', '1.0,1.0,1.0', '-i', refInputFolder + 'dwi_GIBBS/' + dataDenoisedNHDR, \
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
        
        print '------------------ JOIN ALL TARGETS -----------------------'
        fullPrediction =  refInputFolder + 'fullPrediction_GIBBS/fullPredictionDWI.nhdr'
        fullPrediction_flipped =  refInputFolder + 'fullPrediction/fullPredictionDWI_flipped.nhdr'
        if not os.path.exists(os.path.dirname(fullPrediction)):
                os.mkdir(os.path.dirname(fullPrediction))
        if not os.path.exists(os.path.dirname(fullPrediction)+ '/tmp/'):
                os.mkdir(os.path.dirname(fullPrediction)+ '/tmp/')

        print '------------------ RESAMPLE ALL TARGETS -----------------------'
        fullCombineInput = []
        scanNumbers = {'prisma_st':[],'prisma_sa':[],'connectom_st':[],'connectom_sa':[]}
        prevPTR = 0
        for seq in sequences:
            for res in ['sa','st']:

                # # set paths
                inputFolder = args.dir + subj + '/' + seq + '/' + res + '/' 
                targetNHDR = inputFolder + 'dwi/' + subj + '_' + seq + '_' + res + '_dwi.nhdr'
                
                # read in number of DWI images
                fo = open(targetNHDR)
                lines = fo.readlines()
                fo.close()
                for ll in lines:
                    if ll.find('sizes:') > -1:
                        numScans = int(ll.split(' ')[-1][:-1])
                        print 'Num scans:' + str(numScans)
                        scanNumbers[ seq + '_' + res ] = range(prevPTR,(prevPTR+numScans))
                        prevPTR += numScans

                # resample
                if not os.path.exists(os.path.dirname(fullPrediction) + '/tmp/' + os.path.basename(targetNHDR)):
                    call(['crlResampler2', '-g', refInputFolder + 'dwi/' + dataNHDR, \
                                            '-i', targetNHDR, \
                                            '-o' , os.path.dirname(fullPrediction) + '/tmp/' + os.path.basename(targetNHDR), \
                                            '--interp nearest', '-p', args.threads])

                fullCombineInput.append( os.path.dirname(fullPrediction) + '/tmp/' + os.path.basename(targetNHDR) )

        print '------------------ COMBINE ALL TARGETS -----------------------'
        if not os.path.exists(fullPrediction):
            call(['crlDWICombineAcquisitions', '-i', fullCombineInput[0],\
                                        '-i', fullCombineInput[1],\
                                        '-i', fullCombineInput[2],\
                                        '-i', fullCombineInput[3],\
                                        '-o',  fullPrediction,'--nonormalize'])

        if not os.path.exists(fullPrediction_flipped):
            call(['/home/ch137122/bin/crlDWIPrepare', '--open', fullPrediction, '--write',fullPrediction_flipped, '--flipAxis y -applytogradients 1 -applytoimage 0' ])


        print '------------------------  RUNNING DIAMOND  ------------------------'
        diamondOutput = os.path.dirname(fullPrediction) + '/' + dataDenoisedNHDR[:-5] + '_diamondcylDIAMOND3T.nrrd'


        # run diamond and apply reconstruction
        predictedNHDR = runDIAMOND( os.path.basename(upsampledNHDR) , \
                                os.path.dirname(upsampledNHDR)+ '/', \
                                upsampledMask, \
                                diamondOutput, \
                                fullPrediction_flipped, \
                                numThreads=args.threads )

        predictedDWIFiles = getListofDWIFiles(fullPrediction_flipped)







        # ------------------ regenerate DWI for each case ----------------------- # 
        print '------------------ regenerate DWI for each case -----------------------'
        for seq in sequences:
            
            for res in ['sa','st']:

                ## TODO separate files
                # relevant file paths
                inputFolder = args.dir + subj + '/' + seq + '/' + res + '/' 
                outputFolder = inputFolder +  'fulldiamondcylDIAMOND_allpred_denoised_GIBBS_iso1mm/'
                targetNHDR = inputFolder + 'dwi/' + subj + '_' + seq + '_' + res + '_dwi.nhdr'
                recNHDR = 'diamondREC_' + subj + '_' + args.seq_ref + '_' +  args.res_ref +  'diamondcyl_GIBBS_denoised_iso1mm_2_' + seq + '_' + res + '.nhdr'

                try:
                    os.makedirs( outputFolder )
                except:
                    print '\t-Folder ' + outputFolder
                try:
                    os.makedirs( outputFolder + 'tmp/' )
                except:
                    print '\t-Folder ' + outputFolder + 'tmp/ already exists'

                print 'Mapping: ' + subj + '_' + args.seq_ref + '_' +  args.res_ref +  '_GIBBS_denoised_iso1mm'  + \
                        ' to: ' + subj + '_' + seq + '_' + res

                print '------------------------  REMOVE GRADIENTS  ------------------------'
                print 'Total num scans:' + str(prevPTR)
                listofRMGradients = ''
                boolIx = np.ones([prevPTR])
                boolIx[scanNumbers[ seq + '_' + res ]] = 0
                for ii in range(prevPTR):
                    if boolIx[ii] == 1:
                        listofRMGradients += str(ii)+', '
                listofRMGradients = listofRMGradients[:-2]
                # print 'Remove gradients: ' + listofRMGradients
                # if not os.path.exists( outputFolder + recNHDR[:-5] + '_registered.nhdr' ):
                call(['crlDWIRemoveGradientImage','-i',diamondOutput[:-5] + '_predicted.nhdr' ,\
                                                '-o',outputFolder + 'tmp/' + recNHDR,\
                                                '-r', listofRMGradients ])





                print '------------------------  RUNNING RESAMPLER  ------------------------'
                call(['crlResampler2', '-g', targetNHDR, '-i', outputFolder + 'tmp/' +recNHDR , '--interp sinc',\
                                        '-o', outputFolder + recNHDR[:-5] + '_registered.nhdr'])


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
