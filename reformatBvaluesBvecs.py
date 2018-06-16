
import argparse
import os
from subprocess import call

###  MAIN CODE
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Benchmark code for DWI data')
    parser.add_argument('-i', '--dwi', required=True,
                        help='path to the dwi file')
    parser.add_argument('-m', '--mask', required=True,
                        help='path to the dwi file')
    parser.add_argument('-bvec', '--bvec', required=True,
                        help='path to the bvecs file')
    parser.add_argument('-bval', '--bval', required=True,
                        help='path to the bvecs file')
    parser.add_argument('-d', '--dir', required=True,
                        help='path to the files')
    args = parser.parse_args()


    for res in ['sa', 'st']:

        inputFolder = args.dir + res + '/'
            
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

        # reformat bvecs
        fbvec = open( inputFolder + args.bvec )
        lines = fbvec.readlines()
        fbvec.close()

        os.rename(inputFolder + args.bvec, inputFolder + args.bvec +'.bkup' )

        singleLine = []
        for ll in lines:
            singleLine.append( ll.replace( '\t',' ' ) )

        fbvec = open( inputFolder + args.bvec ,'w')
        fbvec.writelines(singleLine)
        fbvec.close()

        # create NHDR files
        if not os.path.exists(inputFolder + 'dwi/' ):
            os.makedirs(inputFolder + 'dwi/' )
        if not os.path.exists( inputFolder + 'dwi/dwi.nhdr' ):
            call([ 'crlDWIConvertFSLToNHDR', '-i', inputFolder + args.dwi, '-o', inputFolder + 'dwi/dwi.nhdr' , '--bvals', inputFolder + args.bval, '--bvecs', inputFolder + args.bvec])

        # reformat masks
        if not os.path.exists( inputFolder + 'mask.nrrd' ):
            call([ '/opt/el7/pkgs/crkit/nightly/20170601/crkit/bin/crlConvertBetweenFileFormats', '-in', inputFolder + args.mask ,'-out', inputFolder + 'mask.nrrd'  ])

        # reformat MPRAGE
        if not os.path.exists( args.dir + 'mprage.nrrd' ):
            call([ '/opt/el7/pkgs/crkit/nightly/20170601/crkit/bin/crlConvertBetweenFileFormats', '-in', args.dir + 'mprage.nii.gz' ,'-out', args.dir + 'mprage.nrrd'  ])