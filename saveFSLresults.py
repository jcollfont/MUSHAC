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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Benchmark code for DWI data')
    parser.add_argument('-d', '--dir', required=True,
                        help='path to subject folders')
    parser.add_argument('-o', '--outdir', required=True,
                        help='path to save folders')
    parser.add_argument('-s', '--subj', required=True,
                        help='subject ID')
    args = parser.parse_args()

    if not os.path.exists(args.dir +'/oriented/'):
        os.mkdir( args.dir +'/oriented/')

    for scanModel in ['prisma', 'connectom']:
        for scanType in ['st', 'sa']:

            baseFolder = args.dir + '/' + args.subj + '/'+ scanModel +'/'+ scanType +'/'
            loadFolder = baseFolder +'fullDIAMOND_allpred_denoised_iso1mm/'
            dwiFiles = getListofDWIFiles( loadFolder + 'diamondREC_'+args.subj+'_prisma_st_denoised_iso1mm_2_'+ scanModel +'_'+ scanType +'_registered.nhdr')
            for fi in dwiFiles:
                call([ 'crlOrientImage', '-r', args.dir + args.subj + '/prisma/st/dwi/' + args.subj + '_prisma_st_denoised_dwi_0000.nrrd', \
                                    loadFolder + fi, args.dir +'/oriented/' +fi ])


            shutil.copy( loadFolder + 'diamondREC_'+args.subj+'_prisma_st_denoised_iso1mm_2_'+ scanModel +'_'+ scanType +'_registered.nhdr', \
                           args.dir +'/oriented/' + 'oriented.nhdr' )


            # modify NHDR
            fo = open(args.dir +'/oriented/targetInfo.txt','w+')
            call(['crlImageInfo', baseFolder +'dwi.nii.gz'], stdout=fo)
            fo.close()
            fo = open(args.dir +'/oriented/targetInfo.txt','rb')
            infolines = fo.readlines()
            fo.close()

            spacing = np.diag(np.array(np.matrix( infolines[4].split('[')[-1].split(']')[0] )).ravel()[:-1])
            M = spacing.dot( np.matrix( infolines[7] + ';' + infolines[8] + ';' +infolines[9]  )[:,:-1] )

            fo = open(loadFolder + 'diamondREC_'+args.subj+'_prisma_st_denoised_iso1mm_2_'+ scanModel +'_'+ scanType +'_registered.nhdr', 'rb')
            lines = fo.readlines()
            fo.close()
            print M

            fo = open(args.dir +'/oriented/' + 'oriented.nhdr', 'w+')
            for ll in lines:
                if ll.find('space directions:') > -1:
                    fo.writelines( 'space directions: (%0.6f,%0.6f,%0.6f) (%0.6f,%0.6f,%0.6f) (%0.6f,%0.6f,%0.6f) none\n' %( M[0,0],M[0,1],M[0,2],\
                                                                                                                        M[1,0],M[1,1],M[1,2],\
                                                                                                                        M[2,0],M[2,1],M[2,2] ))
                elif ll.find('space origin:') > -1:
                    fo.writelines( 'space origin: (' + infolines[5].split('[')[-1].split(']')[0][:-3] + ')\n')
                elif ll.find('line skip:') > -1:
                    fo.writelines('line skip: 14\n')
                else:
                    fo.writelines(ll)
            fo.close()


            if not os.path.exists(args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'):
                os.makedirs(args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/')

            call(['crlDWIConvertNHDRForFSL', '-i',  args.dir +'/oriented/' + 'oriented.nhdr',\
                                        '--data', args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'+args.subj+'_'+ scanModel +'_'+ scanType +'_dwi.nii',\
                                        '--bvecs', args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'+'dwi.bvecs',\
                                        '--bvals', args.outdir + args.subj+'/'+ scanModel +'/'+ scanType +'/'+'dwi.bvals'])


    # shutil.rmtree(args.dir +'/oriented/')


    # call(['crlDWIConvertNHDRForFSL', '-i', args.dir + '/' + args.subj + '/prisma/sa/fullDIAMOND_allpred_denoised_GIBBS_iso1mm/diamondREC_'+args.subj+'_prisma_st_GIBBS_denoised_iso1mm_2_prisma_sa_registered.nhdr',\
    #                             '--data', args.dir + args.outdir + args.subj+'/prisma/sa/'+args.subj+'_prisma_sa_dwi.nii'])

    # call(['crlDWIConvertNHDRForFSL', '-i', args.dir + '/' + args.subj + '/connectom/st/fullDIAMOND_allpred_denoised_GIBBS_iso1mm/diamondREC_'+args.subj+'_prisma_st_GIBBS_denoised_iso1mm_2_connectom_st_registered.nhdr',\
    #                             '--data', args.dir + args.outdir + args.subj+'/connectom/st/'+args.subj+'_connectom_st_dwi.nii'])

    # call(['crlDWIConvertNHDRForFSL', '-i', args.dir + '/' + args.subj + '/connectom/sa/fullDIAMOND_allpred_denoised_GIBBS_iso1mm/diamondREC_'+args.subj+'_prisma_st_GIBBS_denoised_iso1mm_2_connectom_sa_registered.nhdr',\
    #                             '--data', args.dir +  args.outdir + args.subj+'/connectom/sa/'+args.subj+'_connectom_sa_dwi.nii'])

