
import argparse
import os
from subprocess import call


###  MAIN CODE
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Benchmark code for DWI data')
    parser.add_argument('-d', '--dir', required=True,
                        help='path to the files')
    parser.add_argument('-s', '--subj', required=True,
                        help='path to the dwi file')
    parser.add_argument('--scan', '--scan', default='',
                        help='path to the dwi file')
    parser.add_argument('-i', '--dwi', default='dwi.nii.gz',
                        help='path to the dwi file')
    parser.add_argument('-m', '--mask', default='mask.nii.gz',
                        help='path to the dwi file')
    parser.add_argument('-bvec', '--bvec', default='dwi.bvec',
                        help='path to the bvecs file')
    parser.add_argument('-bval', '--bval', default='dwi.bval',
                        help='path to the bvecs file')
    args = parser.parse_args()

    if args.scan == '':
        scan = ['connectom', 'prisma']
    else:
        scan = [args.scan]

    for sc in scan:

        for res in ['sa', 'st']:

            inputFolder = args.dir + '/' + args.subj + '/' + sc + '/' + res + '/'
            
            # read in bvals file
            fbvals = open( inputFolder + args.bval )
            lines = fbvals.readlines()
            fbvals.close()

            os.rename(inputFolder + args.bval, inputFolder + args.bval +'.bkup' )

            singleLine = ''
            for ll in lines:
                singleLine += ll.replace( '\t',' ' )

            fbvals = open( inputFolder + args.bval ,'w')
            fbvals.write(singleLine)
            fbvals.close()

            numBval = len(singleLine.split(' '))
        

            # reformat bvecs
            fbvec = open( inputFolder + args.bvec )
            lines = fbvec.readlines()
            fbvec.close()

            # read in B-vectors (formated in colums)
            bvecs = []
            for ll in lines:
                splitLine =  ll.split('\t')[:-1]
                if len(splitLine) == 0:
                    splitLine =  ll.split(' ')[:-1]
                bvecs.append( splitLine )
            
            print 'Num Bvals: ' + str(numBval)
            print 'Num rows: ' + str(len(bvecs))
            print 'Num cols: ' + str(len(bvecs[0]))

            if len(bvecs) > 3:     # if there are as many rows as B-values (vectors are stored in rows)
                # prepare new lines for bvecs (now in rowss)
                newLines = range(3)
                for dd in range(3):
                    newLines[dd] = ''
                    for bb in bvecs:
                        newLines[dd] += bb[dd]
                    newLines[dd] +='\n'
            else:               # the vectors are stored in columns (only clean tabs)
                newLines = []
                for ll in lines:
                    newLines.append( ll.replace( '\t',' ' ) )

            os.rename(inputFolder + args.bvec, inputFolder + args.bvec +'.bkup' )

            # write
            fo = open(inputFolder + args.bvec,'w+')
            fo.writelines(newLines)
            fo.close()
        

            # create NHDR files
            # if not os.path.exists(inputFolder + 'dwi_new/' ):
            #     os.makedirs(inputFolder + 'dwi_new/' )
            # # if not os.path.exists( inputFolder + 'dwi_new/' + args.subj + '_' + sc + '_' + res + '_dwi.nhdr' ):
            # call([ 'crlDWIConvertFSLToNHDR', '-i', inputFolder + args.dwi \
            #                                 ,'-o', inputFolder  + 'dwi_new/' + args.subj + '_' + sc + '_' + res + '_dwi.nhdr' \
            #                                 ,'--bvals', inputFolder + args.bval \
            #                                 ,'--bvecs', inputFolder + args.bvec\
            #                                 ,'--mirrorGY'])

            if (sc == 'prisma')&(res=='st'):
                # reformat masks
                functionCRLConvertBetweenFileFormats = '/opt/el7/pkgs/crkit/nightly/20160503/bin/crlConvertBetweenFileFormats'
                if not os.path.exists( inputFolder + 'mask.nrrd' ):
                    call([ functionCRLConvertBetweenFileFormats, '-in', inputFolder + args.mask ,'-out', inputFolder + 'mask.nrrd'  ])

            # # reformat MPRAGE
            # if not os.path.exists( args.dir + '/'  + args.subj + '/' + 'mprage.nrrd' ):
            #     call([ functionCRLConvertBetweenFileFormats, '-in', args.dir + '/'  + args.subj + '/' + 'mprage.nii.gz' ,'-out', args.dir + 'mprage.nrrd'  ])
    