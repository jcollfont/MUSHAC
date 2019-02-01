import os
import shutil
from subprocess import call
import sys
import argparse
import numpy as np

def getListofDWIFiles( nhdrFile ):

    dwiFiles = []

    # read in all data file
    dataFi = open(nhdrFile,'r')
    dataLines = dataFi.readlines()
    dataFi.close()

    # reopen to write final maks file
    readFiles = False
    for ll in dataLines:
        if readFiles:
            dwiFiles.append( ll[:-1] )
        if ll.find('data file: LIST') > -1:
           readFiles = True


    return dwiFiles

def computeParamsFromSingleTensorFromDWI(recNHDR, bvecs=[], bvals=[], outputName=None, anatMask=None):

        # tempSave new signal
        if outputName == None:
            tmpdir = '/tmp/tmp%d/'  %(np.random.randint(1e6))
            baseName = ''
        else:
            tmpdir = os.path.dirname(outputName)
            baseName = os.path.basename(outputName)

        try:
            os.makedirs(tmpdir)
        except:
            print tmpdir + ' already exists'

        # compute single tensor on data
        # if not os.path.exists(tmpdir + '/' + baseName + '_recDTI.nrrd'):
        call(['tend', 'estim', '-i', recNHDR, '-o', tmpdir + '/' + baseName + '_recDTI.nrrd' ,  '-B', 'kvp', '-knownB0', 'false'  ])

        # compute FA, MD
        # if not os.path.exists(tmpdir +'/'+ baseName + '_fracAnis.nrrd'):
        call(['crlTensorScalarParameter', '-m',  tmpdir +'/'+ baseName + '_meanDiff.nrrd', '-f',  tmpdir +'/'+ baseName + '_fracAnis.nrrd', tmpdir + '/' + baseName + '_recDTI.nrrd'])

        # load results
        meanDiff =  tmpdir +'/'+  baseName + '_meanDiff.nrrd'
        fracAnisotropy = tmpdir +'/'+ baseName + '_fracAnis.nrrd'

        # # load DWI data
        # dwiData = loadDWIdata( os.path.dirname(recNHDR) + '/', os.path.basename(recNHDR) )[0]

        # if not anatMask == None:
        #     mask = nrrd.read(anatMask)[0]
        #     dwiData = dwiData* np.tile( mask , [dwiData.shape[-1],1,1,1]).transpose(1,2,3,0)
        #
        # compute Mean Kurtosis (MK)
        # if not os.path.exists(tmpdir +'/'+ baseName + '_meanKurtosis.nrrd'):
        #     dkiModel = dki.DiffusionKurtosisModel(self.gtab)
        #     dkifit = dkiModel.fit(dwiData)
        #     meanKurtosis = dkifit.mk()
        #     nrrd.write( tmpdir +'/'+ baseName + '_meanKurtosis.nrrd', meanKurtosis )
        # else:
        #     meanKurtosis = nrrd.read( tmpdir +'/'+ baseName + '_meanKurtosis.nrrd')[0]
        #
        # Compute Return to Origin Probability (RTOP)
        # if not os.path.exists(tmpdir +'/'+ baseName + '_rtop.nrrd'):
        #     dsmodel = DiffusionSpectrumModel( self.gtab )
        #     rtop = dsmodel.fit(  dwiData / np.mean( dwiData[:,:,:, self.gtab.b0s_mask ] ,axis=3, keepdims=True)  ).rtop_pdf()
        #     nrrd.write( tmpdir +'/'+ baseName + '_rtop.nrrd', rtop )
        # else:
        #     rtop = nrrd.read( tmpdir +'/'+ baseName + '_rtop.nrrd')[0]

        if outputName == None:
            shutil.rmtree(tmpdir)


        return meanDiff, fracAnisotropy#, meanKurtosis, rtop

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Benchmark code for DWI data')
    parser.add_argument('-d', '--dir', required=True,
                        help='path to subject folders')
    parser.add_argument('-o', '--outdir', required=True,
                        help='path to save folders')
    parser.add_argument('-s', '--subj', required=True,
                        help='subject ID')
    parser.add_argument('-denoise', '--denoise', default='gibbs',
                        help='Binary option. Should we use gibs sampler (default is true)')
    parser.add_argument('-diamond', '--diamond', default='diamondNCcyl',
                        help='Which type of diamond is to be used?')
    parser.add_argument('-sumToOne', '--sumToOneFlag', default='1',
                        help='Use the sum-to-one option (default=1)?')
    args = parser.parse_args()

    if not os.path.exists(args.dir +'/oriented/'):
        os.mkdir( args.dir +'/oriented/')

    for scanModel in ['prisma', 'connectom']:
        for scanType in ['st', 'sa']:

            # get reconstructed NHDR files
            if args.denoise == 'gibbs':
                denoised_Tag = 'denoised_GIBBS'
            else:
                denoised_Tag = 'denoised'

            baseDir = args.dir + '/'
            subj = args.subj
            diamond_Tag = args.diamond 
            sumToOneFlag = args.sumToOneFlag


            diamondFolder = baseDir + subj + '/' + scanModel + '/' + scanType + '/predicted_' + denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_dwi/'
            diamondName = diamondFolder + subj + '_' + scanModel + '_' + scanType + 'predicted_' + denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_dwi.nhdr'


            # # evaluate FA and MD on image reconstruction
            # outputDTIFolder = baseDir + '/' + subj + '/' + scanModel + '/' + scanType + '/dtiResults_' + denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_dwi/'
            # outputDTIName = outputDTIFolder + subj + '_' + scanModel + '_' + scanType + 'dtiResults_' + denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_dwi.nhdr'
            
            # if not os.path.exists(outputDTIFolder):
            #     os.makedirs(outputDTIFolder)
            
            # meanDiff, fracAnisotropy = computeParamsFromSingleTensorFromDWI(diamondName, outputName=outputDTIName)
            
            # save to FSL format
            outputFSLFolder = args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'
            outputFSLName = args.subj+'_'+ scanModel +'_'+ scanType +'_dwi.nii'
            
            if not os.path.exists(outputFSLFolder):
                os.makedirs(outputFSLFolder)

            call(['crlDWIConvertNHDRForFSL', '-i', diamondName,\
                                        '--data', outputFSLFolder + outputFSLName,\
                                        '--bvecs', outputFSLFolder + 'dwi.bvecs',\
                                        '--bvals', outputFSLFolder + 'dwi.bvals'])


            if (scanModel == 'connectom'):

                diamondFolder = baseDir + subj + '/' + scanModel + '/' + scanType + '/predicted_' + denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_svGrad_dwi/'
                diamondName = diamondFolder + subj + '_' + scanModel + '_' + scanType + 'predicted_' + denoised_Tag + '_' + diamond_Tag + '_sumToOne' + sumToOneFlag + '_svGrad_dwi.nhdr'

                # save to FSL format
                outputFSLFolder = args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'
                outputFSLName = args.subj+'_'+ scanModel +'_'+ scanType +'_svGrad_dwi.nii'
                
                if not os.path.exists(outputFSLFolder):
                    os.makedirs(outputFSLFolder)

                call(['crlDWIConvertNHDRForFSL', '-i', diamondName,\
                                            '--data', outputFSLFolder + outputFSLName,\
                                            '--bvecs', outputFSLFolder + 'dwi_svGrad.bvecs',\
                                            '--bvals', outputFSLFolder + 'dwi_svGrad.bvals'])



















    ## OLD VERSION
    #         baseFolder = args.dir + '/' + args.subj + '/'+ scanModel +'/'+ scanType +'/'
    #         loadFolder = baseFolder +'fullDIAMOND_allpred_denoised_iso1mm/'
    #         dwiFiles = getListofDWIFiles( loadFolder + 'diamondREC_'+args.subj+'_prisma_st_denoised_iso1mm_2_'+ scanModel +'_'+ scanType +'_registered.nhdr')
    #         for fi in dwiFiles:
    #             call([ 'crlOrientImage', '-r', args.dir + args.subj + '/prisma/st/dwi/' + args.subj + '_prisma_st_denoised_dwi_0000.nrrd', \
    #                                 loadFolder + fi, args.dir +'/oriented/' +fi ])


    #         shutil.copy( loadFolder + 'diamondREC_'+args.subj+'_prisma_st_denoised_iso1mm_2_'+ scanModel +'_'+ scanType +'_registered.nhdr', \
    #                        args.dir +'/oriented/' + 'oriented.nhdr' )


    #         # modify NHDR
    #         fo = open(args.dir +'/oriented/targetInfo.txt','w+')
    #         call(['crlImageInfo', baseFolder +'dwi.nii.gz'], stdout=fo)
    #         fo.close()
    #         fo = open(args.dir +'/oriented/targetInfo.txt','rb')
    #         infolines = fo.readlines()
    #         fo.close()

    #         spacing = np.diag(np.array(np.matrix( infolines[4].split('[')[-1].split(']')[0] )).ravel()[:-1])
    #         M = spacing.dot( np.matrix( infolines[7] + ';' + infolines[8] + ';' +infolines[9]  )[:,:-1] )

    #         fo = open(loadFolder + 'diamondREC_'+args.subj+'_prisma_st_denoised_iso1mm_2_'+ scanModel +'_'+ scanType +'_registered.nhdr', 'rb')
    #         lines = fo.readlines()
    #         fo.close()
    #         print M

    #         fo = open(args.dir +'/oriented/' + 'oriented.nhdr', 'w+')
    #         for ll in lines:
    #             if ll.find('space directions:') > -1:
    #                 fo.writelines( 'space directions: (%0.6f,%0.6f,%0.6f) (%0.6f,%0.6f,%0.6f) (%0.6f,%0.6f,%0.6f) none\n' %( M[0,0],M[0,1],M[0,2],\
    #                                                                                                                     M[1,0],M[1,1],M[1,2],\
    #                                                                                                                     M[2,0],M[2,1],M[2,2] ))
    #             elif ll.find('space origin:') > -1:
    #                 fo.writelines( 'space origin: (' + infolines[5].split('[')[-1].split(']')[0][:-3] + ')\n')
    #             elif ll.find('line skip:') > -1:
    #                 fo.writelines('line skip: 14\n')
    #             else:
    #                 fo.writelines(ll)
    #         fo.close()


    #         if not os.path.exists(args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'):
    #             os.makedirs(args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/')

    #         call(['crlDWIConvertNHDRForFSL', '-i',  args.dir +'/oriented/' + 'oriented.nhdr',\
    #                                     '--data', args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'+args.subj+'_'+ scanModel +'_'+ scanType +'_dwi.nii',\
    #                                     '--bvecs', args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'+'dwi.bvecs',\
    #                                     '--bvals', args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'+'dwi.bvals'])


    # # shutil.rmtree(args.dir +'/oriented/')


    # # call(['crlDWIConvertNHDRForFSL', '-i', args.dir + '/' + args.subj + '/prisma/sa/fullDIAMOND_allpred_denoised_GIBBS_iso1mm/diamondREC_'+args.subj+'_prisma_st_GIBBS_denoised_iso1mm_2_prisma_sa_registered.nhdr',\
    # #                             '--data', args.dir + args.outdir + args.subj+'/prisma/sa/'+args.subj+'_prisma_sa_dwi.nii'])

    # # call(['crlDWIConvertNHDRForFSL', '-i', args.dir + '/' + args.subj + '/connectom/st/fullDIAMOND_allpred_denoised_GIBBS_iso1mm/diamondREC_'+args.subj+'_prisma_st_GIBBS_denoised_iso1mm_2_connectom_st_registered.nhdr',\
    # #                             '--data', args.dir + args.outdir + args.subj+'/connectom/st/'+args.subj+'_connectom_st_dwi.nii'])

    # # call(['crlDWIConvertNHDRForFSL', '-i', args.dir + '/' + args.subj + '/connectom/sa/fullDIAMOND_allpred_denoised_GIBBS_iso1mm/diamondREC_'+args.subj+'_prisma_st_GIBBS_denoised_iso1mm_2_connectom_sa_registered.nhdr',\
    # #                             '--data', args.dir +  args.outdir + args.subj+'/connectom/sa/'+args.subj+'_connectom_sa_dwi.nii'])

