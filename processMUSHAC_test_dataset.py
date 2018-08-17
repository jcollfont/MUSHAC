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
def runDenoising( inputDataList, outputNameList ):

    for inputData in inputDataList:

        outputName = outputNameList.pop(0)

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
def runDIAMOND( dataNHDRList, inputFolderList, mask, outputNameList, numThreads=50):

    for dataNHDR in dataNHDRList:

        inputFolder = inputFolderList.pop(0)
        outputName = outputNameList.pop(0)

        # define which DIAMOND code to call
        functionName = 'crlDCIEstimate'
        if np.any( os.uname()[1] == [ 'galahad','lancelot','percival' ] ):
            functionName = '/home/ch199899/bins/crlDCIEstimate_March2018_XeonPhi'

        # run DIAMOND
        if not os.path.exists( outputName[:-5] + '_t0.nrrd' ):
            print 'Computing ' + outputName[:-5] + '_t0.nrrd'
            try:
                out = call([functionName, '-i',  inputFolder + dataNHDR, '-m', mask,'-o', outputName,  \
                            '--residuals' , '-n 3', '-p ' + str(numThreads) , '--automose aicu', '--fascicle diamondcyl'] )
            except :
                print 'Could not run DIAMOND on ' + dataNHDR
        else:
            print outputName[:-5] + '_t0.nrrd'  + ' already computed'


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

def resample2HighRes( headerFileList, inputPathList, outputPathList, maskPath, resolution ,nThreads=50):

    outName = []
    for headerFile in headerFileList:

        inputPath = inputPathList.pop(0)
        outputPath = outputPathList.pop(0)

        # update header file
        fo = open( inputPath + headerFile )
        lines = fo.readlines()
        fo.close()

        # retrieve files reated to the given header
        dwiFiles = []
        readFiles = False
        for ll in lines:
            if  readFiles:
                dwiFiles.append( ll[:-1] )
            if 'data file: LIST' in ll:
                readFiles = True

        # resample files
        newDWIfiles = []
        for ff in dwiFiles:

            # update resolution of individual nrrds
            resolutionSTR = '%0.2f,%0.2f,%0.2f' %(resolution, resolution, resolution)
            newDWIfiles.append( ff[:-9] + 'isores%dmm_' %(resolution*10) + ff[-9:] )

            if not os.path.exists( outputPath +  newDWIfiles[-1] ):
                out = call(['crlResampler2', '--voxelsize', resolutionSTR, '-i', inputPath + ff, \
                                '-o' , outputPath +  newDWIfiles[-1], '--interp sinc', '-p', str(nThreads)])

         
        # extract new size from an example file
        fr = open( outputPath + newDWIfiles[0])
        rlines = fr.readlines()
        fr.close()
        refSize = [ ll for ll in rlines if 'sizes:' in ll ][0][:-1]
        refSpD = [ ll for ll in rlines if 'space directions' in ll ][0][:-1]
        refOrg = [ ll for ll in rlines if 'space origin:' in ll ][0][:-1]

        # write 
        for ll in range(len(lines)):
            if 'sizes' in lines[ll]:
                numBval = lines[ll].split(' ')[-1]
                lines[ll] = refSize + ' ' + numBval 
            if 'space directions' in lines[ll]:
                lines[ll] = refSpD + ' none\n'
            if 'space origin:' in  lines[ll]:
                lines[ll] = refOrg + '\n'

            if 'data file: LIST' in lines[ll]:
                for l2 in range(ll+1, len(lines)):
                    lines[l2] = newDWIfiles.pop(0) + '\n'
                break

        outName.append( headerFile[:-5] + '_isores%dmm.nhdr' %(resolution*10) )
        fo = open( outputPath + outName[-1] , 'w+')
        fo.writelines( lines )
        fo.close()
        

    # resample mask 
    outMask = maskPath[:-5] + '_isores%dmm.nrrd' %(resolution*10)
    if not os.path.exists( outMask ):
        out = call(['crlResampler2', '--voxelsize', resolutionSTR, '-i', maskPath, '-o' , outMask, '--interp nearest', '-p', str(nThreads)])

    return outName, outMask


def registerDWI( newDWINHDR, refDWINHDR, outputTransform ):

    newDWI = getListofDWIFiles( newDWINHDR )
    regTarget = os.path.dirname(newDWINHDR) + '/' + newDWI[0]
    regRef = os.path.dirname(refDWINHDR) + '/' + getListofDWIFiles( refDWINHDR )[0]

    regSTR = os.path.basename(outputTransform)
    
    # compute registration transform from first dwi images (hopefully B0)
    if not os.path.exists( outputTransform + '_FinalS.tfm' ):
        call(['crlBlockMatchingRegistration', '-r', regRef, '-f',  regTarget, '-o', outputTransform, '-t affine'])#\
                                        # '-s 4', '-e 0', '-k 0.8', '--sig 2.5', '--mv 0.0', '-n 10', \
                                        # '--bh 20', '--blv 2', '--rs 1', '--ssi cc', '-I linear', '-p 200', '-t affine'])

    # register all DWI images
    # call(['crlResampler2', '-i', refDWINHDR, '-o', refDWINHDR[:-5] + '.nhdr', \
    #                     '--interp', 'sinc' , '-g', regTarget, '-t', outputTransform + '_FinalS.tfm' , '-p', '50'])
    regFiles = []
    for fi in newDWI:
        regTarget = os.path.dirname(newDWINHDR) + '/' + fi
        if not os.path.exists( regTarget[:-5] + '_registered_'+ regSTR + '_FinalS.nrrd' ):
            call(['crlBlockMatchingRegistration', '-r', regRef, '-f',  regTarget, '-o',  regTarget[:-5] + '_registered_'+ regSTR , \
                                            '-i', outputTransform + '_FinalS.tfm', '-N', '-n 0', '-e 0', '-s 0' ,'-p 1', '-t affine'])
        regFiles.append( os.path.basename( regTarget[:-5] ) + '_registered_'+ regSTR + '_FinalS.nrrd\n' )

    fi = open(refDWINHDR)
    lines = fi.readlines()
    fi.close()

    newLines = lines[:-len(regFiles)]
    newLines += regFiles

    fi = open(newDWINHDR[:-5] + '_registered_'+ regSTR +'.nhdr', 'w+')
    fi.writelines(newLines)
    fi.close()

    return newDWINHDR[:-5] + '_registered_'+ regSTR +'.nhdr'



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

    # set temp folder
    # tmpdir = '/tmp/tmp%d'  %(np.random.randint(1e6))
    # try:
    #     os.makedirs(tmpdir)
    # except:
    #     print tmpdir + ' already exists'
    # print 'temp save at: ' + tmpdir

    # for all subjects
    for subj in subjects:

        print 'Subject ' + subj 

        # input data names and paths
        refInputFolder = args.dir + subj + '/' + args.seq_ref + '/' + args.res_ref + '/' 
        dataNHDR = subj + '_' + args.seq_ref + '_' + args.res_ref + '_dwi.nhdr'
        dataDenoisedNHDR = subj + '_' + args.seq_ref + '_' + args.res_ref + '_denoised_dwi.nhdr'
        mask = refInputFolder + 'mask.nrrd'


        # ------------------ DENOISE DATA ----------------------- #
        runDenoising( [refInputFolder + 'dwi/' + dataNHDR],  [refInputFolder + 'dwi/' + dataDenoisedNHDR] )

        
        # ------------------ RUNN DIAMOND DATA ----------------------- #
        if not os.path.exists(refInputFolder + 'DIAMOND/' ):
            os.makedirs(refInputFolder + 'DIAMOND/' )
        # set Paths for reference file
        dataFilesList =     [ dataNHDR, \
                            dataDenoisedNHDR ]
        refInputFolderList = [ refInputFolder + 'dwi/',
                            refInputFolder + 'dwi/' ]
        diamondFilesList =  [refInputFolder + 'DIAMOND/' + dataNHDR[:-9] + '_DIAMOND3T.nrrd', \
                            refInputFolder + 'DIAMOND/' +  dataDenoisedNHDR[:-9]  + '_DIAMOND3T.nrrd']

        # run DIAMOND on the reference sequence
        print 'Running DIAMOND on reference data'
        runDIAMOND( dataFilesList[:], refInputFolderList[:], mask, diamondFilesList[:],  numThreads=args.threads )


        # ------------------ UPSAMPLE  ----------------------- # 
        for resol in [1.5, 1.2]:

            print 'Resampling data to %0.1fmm' %(resol)
            if not os.path.exists(refInputFolder + 'dwi_iso%dmm/' %(resol*10)):
                os.mkdir(refInputFolder + 'dwi_iso%dmm/' %(resol*10))
            ref15ResolutionFolder = [ refInputFolder + 'dwi_iso%dmm/' %(resol*10) for ii in range(len(dataFilesList)) ]
            dataNHDRupsample, maskupsample = resample2HighRes( dataFilesList[:], refInputFolderList[:], ref15ResolutionFolder[:],  \
                                                                    mask, resol ,nThreads=int(args.threads))
        
            print 'Running DIAMOND on data at resolution %0.1fmm' %(resol)
            if not os.path.exists(refInputFolder + 'DIAMOND_iso%dmm/' %(resol*10)):
                os.mkdir(refInputFolder + 'DIAMOND_iso%dmm/' %(resol*10))
            diamondFilesList = [ refInputFolder + 'DIAMOND_iso%dmm/' %(resol*10) + dd[:-9] + '_iso%dmm_DIAMOND3T.nrrd' %(resol*10) \
                                    for dd in dataNHDRupsample ]
            runDIAMOND( dataNHDRupsample[:], ref15ResolutionFolder[:], maskupsample, diamondFilesList[:],  args.threads )

        for resol in [1, 2]: #[0,1,2]:

            if resol == 0:
                resolSTR = '_origRes'
                folderResol = ''
                nhdrResol = ''
            elif resol == 1:
                resolSTR = '_iso12mmRes'
                folderResol = '_iso12mm'
                nhdrResol = '_isores12mm'
            elif resol ==2:
                resolSTR = '_iso15mmRes'
                folderResol = '_iso15mm'
                nhdrResol = '_isores15mm'


            for denoise in [0,1]:

                if denoise == 0:
                    # ------------------ DIAMOND model  ----------------------- # 
                    print 'Preparing DIAMOND model for extrapolation'
                    NHDR = dataNHDR[:-5] + nhdrResol + '.nhdr'
                    diamondModel = MUSHACreconstruction( refInputFolder,'DIAMOND'+ folderResol + '/' ,\
                                                        'dwi'+folderResol + '/' + NHDR, maskPath='mask'+nhdrResol+'.nrrd')
                    denoiseSTR = '_dwi'
                    refHeader =  refInputFolder + 'dwi'+folderResol + '/' + NHDR
                elif denoise == 1:
                    # ------------------ DIAMOND model  ----------------------- # 
                    print 'Preparing DIAMOND model for extrapolation'
                    NHDR = dataNHDR[:-5] + nhdrResol + '.nhdr'
                    diamondModel = MUSHACreconstruction( refInputFolder,'DIAMOND'+ folderResol + '/',\
                                                        'dwi'+folderResol + '/' + NHDR, maskPath='mask'+nhdrResol+'.nrrd')
                    denoiseSTR = '_dwi_denoised'
                    refHeader =  refInputFolder + 'dwi' +folderResol + '/' + NHDR


                # ------------------ COMPARISON WITH OTHER MODELS ----------------------- #
                for seq in sequences:
                    
                    for res in ['sa','st']:
                            
                        # # set paths
                        inputFolder = args.dir + subj + '/' + seq + '/' + res + '/' 
                        outputFolder = inputFolder + 'reconstructed' +  denoiseSTR + resolSTR + '/'

                        recNHDR = 'recDWI_' + subj + '_' + args.seq_ref + '_' +  args.res_ref +  denoiseSTR + resolSTR + '_2_' + seq + '_' + res

                        print 'Mapping: ' + subj + '_' + args.seq_ref + '_' +  args.res_ref +  denoiseSTR + resolSTR  + \
                                ' to: ' + subj + '_' + seq + '_' + res

                        try:
                            os.makedirs( outputFolder )
                        except:
                            print '\t-Folder ' + outputFolder
                        
                        bvals, bvecs = readInBvecsAndBvals( inputFolder )                                           # define gradients from target file
                        
                        # ------------------ register DWI from DIAMOND model to new data ----------------------- # 
                        print '\t-generating reconstructed signal...'
                        # if not os.path.exists( outputFolder + recNHDR + '.nhdr' ):
                        diamondReconstruction = diamondModel.generateDWIdata( bvecs, bvals )                        # generate new data
                        saveNRRDwithHeader( diamondReconstruction, refHeader, outputFolder, recNHDR , bvals, bvecs )    # save
                        
                    
                        # ------------------ get FA and MD from DIAMOND model  ----------------------- # 
                        print '\t-computing MD and FA on from DIAMOND model'
                        # if not os.path.exists( outputFolder + recNHDR + 'tensStick0.nrrd' ):
                        sticks, diffusionVec, meanDiffusivity, fractionalAnisotropy, radialDiffusivity, axialDiffusivity = \
                                    diamondModel.computeDIAMONDTensorParams( outputDir=outputFolder, outputName=recNHDR )

                        # ------------------ get FA and MD from DTI model  ----------------------- # 
                        print '\t-computing MD and FA from DTI model'
                        # if not os.path.exists( outputFolder + recNHDR + 'tensStick0.nrrd' ):
                        diamondModel.computeParamsFromSingleTensorFromDWI( \
                                                bvecs=bvecs, bvals=bvals, recSignal=None, \
                                                recNHDR= outputFolder + recNHDR + '.nhdr', \
                                                outputName= outputFolder + recNHDR + '_DTI', \
                                                anatMask= refInputFolder + 'mask'+nhdrResol+'.nrrd')
