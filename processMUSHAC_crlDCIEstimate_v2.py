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



def convertFSLToNHDR( fslData, bvecs, bvals, outputNHDR ):
    ## This back and forth conversion from FSL to NHDR makes sure we are not screwing up the gradient diretctions...


    # create temporary dir
    tmpdir = '/tmp/tmp%d'  %(np.random.randint(1e6))
    os.makedirs(tmpdir)

    # convert once 
    out = call([ 'crlDWIConvertFSLToNHDR', '-i', fslData, '-o', outputNHDR ,'--bvals', bvals, '--bvecs', bvecs])
    
    # convert back to FSL
    out = call([ 'crlDWIConvertNHDRForFSL', '-i', outputNHDR, '--data', tmpdir + '/data.nii' ,\
                        '--bvals',  tmpdir + '/bvals', '--bvecs', tmpdir + '/bvecs'])

    # convert twice
    out = call([ 'crlDWIConvertFSLToNHDR', '-i', tmpdir + '/data.nii', '-o', outputNHDR ,\
                        '--bvals', tmpdir + '/bvals', '--bvecs', tmpdir + '/bvecs'])

    # remove temporary folder
    shutil.rmtree(tmpdir)

    return outputNHDR



#%% run Benoit's denoising algo
def runDenoising( inputData, bvals, bvecs, outputName ):

    if not os.path.exists( outputName ): 
        functionDenoising = '/home/ch137122/Software/mrtrix/mrtrix3/release/bin/dwidenoise'
        functionGibbs = '/home/ch199899/links/DWI/external/mrtrix3/bin/mrdegibbs'
        functionMRconvert = '/home/ch137122/Software/mrtrix/mrtrix3/release/bin/mrconvert'

        tmpdir = '/tmp/tmp%d'  %(np.random.randint(1e6))
        os.makedirs(tmpdir)

        if inputData.endswith('.nrrd'):
            print 'Convert to FSL'
            out = call([ 'crlDWIConvertNHDRForFSL', '-i', inputData, '--data', tmpdir + '/data.nii' ,\
                            '--bvals',  tmpdir + '/bvals', '--bvecs', tmpdir + '/bvecs'])
            inputData =  tmpdir + '/data.nii' 
            bvecs =  tmpdir + '/bvecs'
            bvals =  tmpdir + '/bvals'

        print 'Convert to MIF'
        out = call([ functionMRconvert, '-fslgrad', bvecs, bvals, inputData,  tmpdir + '/data.mif' ])

        print 'run denoise'
        out = call([ functionDenoising, tmpdir + '/data.mif',  tmpdir + '/data-gibbs.mif'])

        print 'run mrdegibbs'
        out = call([ functionGibbs,  tmpdir + '/data-gibbs.mif', tmpdir + '/data-denoised.nii'  ])

        print 'convert back to NHDR'
        convertFSLToNHDR( tmpdir + '/data-denoised.nii', bvecs, bvals, outputName )
        

        shutil.rmtree(tmpdir)
    else:
        print '%s already exists. Not computing denoising algorithm.' %(outputName)

#%% run DIAMOND on benchmark data
def runDIAMOND( dataNHDR,  mask, outputName, targetNHDR=None, numThreads=50, diamondType='diamondNCcyl', sumToOneFlag='1'):


    # define which DIAMOND code to call
    functionName = '/home/ch137122/bin/crlDCIEstimateJaume'
    
    # run DIAMOND
    if not os.path.exists( outputName[:-5] + '_predicted.nhdr'):
        print 'Computing ' + outputName[:-5] + '_t0.nrrd'
        try:
            out = call([functionName, '-i',  dataNHDR, '-m', mask,'-o', outputName,  \
                        '--residuals' , '-n 3', '-p ' + str(numThreads) , '--automose aicu',\
                        '--fascicle', diamondType, '--fractions_sumto1', sumToOneFlag, '--estimateDisoIfNoFascicle 1',\
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
    parser.add_argument('-denoise', '--denoise', default='gibbs',
                        help='Binary option. Should we use gibs sampler (default is true)')
    parser.add_argument('-diamond', '--diamond', default='diamondNCcyl',
                        help='Which type of diamond is to be used?')
    parser.add_argument('-sumToOne', '--sumToOneFlag', default='1',
                        help='Use the sum-to-one option (default=1)?')
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

    if args.denoise == 'gibbs':
        denoised_Tag = 'denoised_GIBBS'
    else:
        denoised_Tag = 'denoised'

    
    
    diamond_Tag = args.diamond 
    sumToOneFlag = args.sumToOneFlag


    # for all subjects
    for subj in subjects:

        print 'Subject ' + subj 

        # input data names and paths
        trainFolder = args.dir + subj + '/' + args.seq_ref + '/' + args.res_ref + '/' 
        trainDataFSL = trainFolder + 'dwi.nii.gz'
        bvalsFSL = trainFolder + 'dwi.bval'
        bvecsFSL = trainFolder + 'dwi.bvec'

        maskFSL = trainFolder +'mask.nii.gz'
        maskNRRD = trainFolder + 'mask.nrrd'

        # ------------------ DENOISE DATA ----------------------- #
        folder_denoised_data = trainFolder + 'dwi_' + denoised_Tag + '/'
        dataDenoisedNHDR = folder_denoised_data + subj + '_' + args.seq_ref + '_' + args.res_ref + '_' + denoised_Tag + '_dwi.nhdr'
    

        if not os.path.exists(folder_denoised_data):
                os.mkdir(folder_denoised_data)
        print '------------------------  RUNNING DENOISING  ------------------------'
        runDenoising( trainDataFSL, bvalsFSL, bvecsFSL, dataDenoisedNHDR )
        # input is in FSL format, output is in NRRD

        

        # ------------------ UPSAMPLE to 1mm**3 ----------------------- # 
        print '------------------------  RUNNING UPSAMPLE  ------------------------'
        folder_upsampled = trainFolder + 'dwi_' + denoised_Tag + '_iso1mm/'
        upsampledNHDR = folder_upsampled  +  subj + '_' + args.seq_ref + '_' + args.res_ref + '_' + denoised_Tag + '_iso1mm_dwi.nhdr'
        if not os.path.exists(folder_upsampled):
                os.mkdir(folder_upsampled)

        # ------------------ UPSAMPLE IMAGE  ----------------------- #
        if not os.path.exists(upsampledNHDR):
            call(['crlResampler2', '--voxelsize', '1.0,1.0,1.0', '-i', dataDenoisedNHDR, \
                                    '-o' , upsampledNHDR, '--interp sinc', '-p', args.threads])

        # ------------------ UPSAMPLE MASK  ----------------------- # 
        upsampledMask = folder_upsampled + 'mask_iso1mm.nrrd'
        if not os.path.exists(upsampledMask):
            call(['crlResampler2', '-g', upsampledNHDR, '-i', maskNRRD, \
                            '-o' , upsampledMask, '--interp nearest', '-p', args.threads])





        # ------------------ prepare prediction targets  ----------------------- # 
        # after this section, each scan type of a subject shoould have an upsampled version (1mm**3) of the original data in nrrd
        print '------------------ JOIN ALL TARGETS -----------------------'
        folder_predictionTarget =  trainFolder + 'predictionTarget/'
        predictionTargetNHDR = folder_predictionTarget + 'predictionTarget.nhdr'

        if not os.path.exists(folder_predictionTarget):
                os.mkdir(folder_predictionTarget)
        if not os.path.exists(folder_predictionTarget + '/tmp/'):
                os.mkdir(folder_predictionTarget + '/tmp/')



        print '------------------ CONVERT TO NRRD AND RESAMPLE ALL TARGETS -----------------------'
        fullCombineInput = []
        scanNumbers = {'prisma_st':[],'prisma_sa':[],'connectom_st':[],'connectom_sa':[]}
        prevPTR = 0
        for seq in sequences:
            for res in ['sa','st']:

                # # set paths
                targetFSL = args.dir + subj + '/' + seq + '/' + res + '/' \
                                 + 'dwi.nii.gz'
                targetBvals = os.path.dirname(targetFSL) + '/dwi.bval'
                targetBvecs = os.path.dirname(targetFSL) + '/dwi.bvec'
                
                # convert to NRRD
                folder_targetNRRD = args.dir + subj + '/' + seq + '/' + res + '/dwi/'
                targetNRRD = folder_targetNRRD + subj + '_' + seq + '_' + res + '_dwi.nhdr'
                if not os.path.exists(folder_targetNRRD):
                    os.mkdir(folder_targetNRRD)
                
                if not os.path.exists(targetNRRD):
                    convertFSLToNHDR( targetFSL, targetBvecs, targetBvals, targetNRRD )


                # read in number of DWI images (to be used later)
                fo = open(targetNRRD)
                lines = fo.readlines()
                fo.close()
                for ll in lines:
                    if ll.find('sizes:') > -1:
                        numScans = int(ll.split(' ')[-1][:-1])
                        print 'Num scans:' + str(numScans)
                        scanNumbers[ seq + '_' + res ] = range(prevPTR,(prevPTR+numScans))
                        prevPTR += numScans



                # resample
                folder_target_upsampledNRRD = args.dir + subj + '/' + seq + '/' + res + '/dwi_iso1mm/'
                target_upsampledNRRD = folder_targetNRRD + subj + '_' + seq + '_' + res + '_iso1mm_dwi.nhdr'
                if not os.path.exists(folder_target_upsampledNRRD):
                    os.mkdir(folder_target_upsampledNRRD)

                if not os.path.exists(target_upsampledNRRD):
                    call(['crlResampler2', '-g', dataDenoisedNHDR, \
                                            '-i', targetNRRD, \
                                            '-o' , target_upsampledNRRD, \
                                            '--interp nearest', '-p', args.threads])

                fullCombineInput.append( target_upsampledNRRD )



        # ------------------ COMBINE TARGET AQUISITIONS  ----------------------- # 
        print '------------------ COMBINE ALL TARGETS -----------------------'
        if not os.path.exists(predictionTargetNHDR):
            call(['/opt/el7/pkgs/crkit/release-current/bin/crlDWICombineAcquisitions', '-i', fullCombineInput[0],\
                                        '-i', fullCombineInput[1],\
                                        '-i', fullCombineInput[2],\
                                        '-i', fullCombineInput[3],\
                                        '-o',  predictionTargetNHDR,'--nonormalize'])





        # ------------------ RUN DIAMOND  ----------------------- #
        print '------------------------  RUNNING DIAMOND  ------------------------'
        diamondOutputFolder = trainFolder  + '/diamond_' + diamond_Tag + '_sumToOne' + sumToOneFlag +  '/'
        diamondOutputName = diamondOutputFolder + subj + '_' + args.seq_ref + '_' + args.res_ref + '_' + \
                                    denoised_Tag +  '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_iso1mm_DIAMOND3T.nrrd'

        if not os.path.exists(diamondOutputFolder):
                    os.mkdir(diamondOutputFolder)

        # run diamond and apply reconstruction
        predictedNHDR = runDIAMOND( upsampledNHDR , \
                                upsampledMask, \
                                diamondOutputName, \
                                predictionTargetNHDR, \
                                args.threads,\
                                diamondType=diamond_Tag,
                                sumToOneFlag=sumToOneFlag)

        predictedDWIFiles = getListofDWIFiles(predictionTargetNHDR)





        # ------------------ regenerate DWI for each case (separate the predicted results ----------------------- # 
        print '------------------ regenerate DWI for each case -----------------------'
        for seq in sequences:
            
            for res in ['sa','st']:

                print 'Mapping: ' + subj + '_' + args.seq_ref + '_' +  args.res_ref +  '_' + denoised_Tag + '_' + diamond_Tag  + '_sumToOne' + sumToOneFlag + \
                        ' to: ' + subj + '_' + seq + '_' + res

                diamondOutput_gradient_removed = folder_target_upsampledNRRD + '/tmp/' + \
                                                        subj + '_' + seq + '_' + res + '_' + \
                                                        denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_iso1mm_dwi.nhdr'
                try:
                    os.makedirs( os.path.dirname(diamondOutput_gradient_removed) )
                except:
                    print '\t-Folder ' + os.path.dirname(diamondOutput_gradient_removed)

                
                print '------------------------  REMOVE GRADIENTS  ------------------------'
                
                # get list of gradients to remove from main file
                print 'Total num scans:' + str(prevPTR)
                listofRMGradients = ''
                boolIx = np.ones([prevPTR])
                boolIx[scanNumbers[ seq + '_' + res ]] = 0
                for ii in range(prevPTR):
                    if boolIx[ii] == 1:
                        listofRMGradients += str(ii)+','
                listofRMGradients = listofRMGradients[:-2]

                # remove gradients from main file
                call(['crlDWIRemoveGradientImage','-i',diamondOutputName[:-5] + '_predicted.nhdr' ,\
                                                '-o', diamondOutput_gradient_removed,\
                                                '-r', listofRMGradients ])




                print '------------------------  RUNNING RESAMPLER  ------------------------'

                # Resample target folder + file
                folder_targetNRRD = args.dir + subj + '/' + seq + '/' + res + '/dwi/'
                targetNRRD = folder_targetNRRD + subj + '_' + seq + '_' + res + '_dwi.nhdr'

                # output ressampled folder + file
                outputFolder = args.dir + subj + '/' + seq + '/' + res + '/predicted_' + denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_dwi/'
                outputPredictionNRRD = outputFolder + subj + '_' + seq + '_' + res + 'predicted_' + denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_dwi.nhdr'

                try:
                    os.makedirs( outputFolder )
                except:
                    print '\t-Folder ' + outputFolder

                call(['crlResampler2',  '-g', targetNRRD, \
                                        '-i', diamondOutput_gradient_removed, '--interp sinc',\
                                        '-o', outputPredictionNRRD ])


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
