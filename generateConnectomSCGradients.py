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

#%%
#
#
#
#
def regenerateDWI( xmlFile, svGradients, targetFile, outputName):
    call(['/home/ch137122/bin/crlDCIGenerateSignal','-i', xmFlile, \
                '-o', outputName, \
                '--scheme',targetFile,\
                '--svgradients', svGradients,\
                '-p 100'])
                        

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
        sequences = ['connectom']
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



        # ------------------ RUN DIAMOND  ----------------------- #
        print '------------------------  RUNNING DIAMOND  ------------------------'
        diamondOutputFolder = trainFolder  + '/diamond_' + diamond_Tag + '_sumToOne' + sumToOneFlag +  '/'
        diamondOutputName = diamondOutputFolder + subj + '_' + args.seq_ref + '_' + args.res_ref + '_' + \
                                    denoised_Tag +  '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_iso1mm_DIAMOND3T.nrrd'

        diamondXML = diamondOutputName[:-5] + '.xml'

        # ------------------ regenerate DWI for each case (separate the predicted results ----------------------- # 
        print '------------------ regenerate DWI for each case -----------------------'
        for seq in sequences:
            
            for res in ['sa','st']:

                # Resample target folder + file
                targetFolder = args.dir + subj + '/' + seq + '/' + res + '/'
                folder_targetNRRD =  targetFolder + '/dwi/'
                targetNRRD = folder_targetNRRD + subj + '_' + seq + '_' + res + '_dwi.nhdr'
                folder_target_upsampledNRRD = args.dir + subj + '/' + seq + '/' + res + '/dwi_iso1mm/'
                target_upsampledNRRD = folder_targetNRRD + subj + '_' + seq + '_' + res + '_iso1mm_dwi.nhdr'

                if not os.path.exists(folder_target_upsampledNRRD):
                    os.makedirs(folder_target_upsampledNRRD)

                print 'Mapping: ' + subj + '_' + args.seq_ref + '_' +  args.res_ref +  '_' + denoised_Tag + '_' + diamond_Tag  + '_sumToOne' + sumToOneFlag + \
                        ' to: ' + subj + '_' + seq + '_' + res

                print '------------------------  RUNNING predictor  ------------------------'
                svGradients = targetFolder + 'dwi_modb_res.nrrd'
                outputName_SVGPrediction = folder_target_upsampledNRRD + '/tmp/' + \
                                        subj + '_' + seq + '_' + res + '_' + \
                                        denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_iso1mm_dwi.nhdr'
                regenerateDWI( diamondXML, svGradients, target_upsampledNRRD, outputName_SVGPrediction)

                print '------------------------  RUNNING RESAMPLER  ------------------------'

                

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
