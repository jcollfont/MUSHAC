#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:09:42 2019

@author: ch199899
"""




### IMPORTS
# standard python
import numpy as np 
import nrrd
import os
import shutil
from subprocess import call
import sys
import argparse
from joblib import Parallel, delayed

def applyResampler( referenceImage, inputFile, outputFile ):
    call(['crlResampler2','-g',  referenceImage,\
                '-i',  inputFile, '-o', outputFile ,'-p 40'])
    return  outputFile


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Benchmark code for DWI data')
    parser.add_argument('-d', '--dir', required=True,
                        help='path to subject folders')
    args = parser.parse_args()

    inputDir = args.dir

    
    for axis in ['x','y','z']:
            
        # retrieve all files
        svGradFile = 'dwi_modb_'+axis+'.nii.gz'

        # create output folder
        tempDir = inputDir + 'tempSCG_'+axis+'/'
        if not os.path.exists(tempDir):
            os.makedirs(tempDir)

        # convert to N3 series
        call(['/home/ch137122/bin/crlConvert4DToN3D','-i', inputDir + svGradFile, '-o', tempDir + svGradFile])
        
        # resample
        files = os.listdir(tempDir)
        # joinFiles = []
        # for ff in files:
        #     call(['crlResampler2','-g',  inputDir + 'N_prisma_st_denoised_GIBBS_diamondNCcyl_sumToOne1_iso1mm_DIAMOND3T_b0.nrrd',\
        #             '-i',  tempDir + ff, '-o', tempDir + ff[:-7] + '_1mm.nrrd' ,'-p 40'])
        #     joinFiles += ['-i', tempDir + ff[:-7] + '_1mm.nrrd'  ]
        b0File = inputDir + 'N_prisma_st_denoised_GIBBS_diamondNCcyl_sumToOne1_iso1mm_DIAMOND3T_b0.nrrd'
        outFiles = Parallel(n_jobs=40)( delayed(applyResampler) \
                (  b0File,\
                    tempDir + ff, \
                    tempDir + ff[:-7] + '_1mm.nrrd' ) for ff in files)

        joinFiles = []
        for ff in outFiles:
            joinFiles += ['-i', ff]

        # move back to N4 dim
        # call(['/home/ch137122/bin/crlConvertN3DTo4D'] + joinFiles + ['-o', inputDir + 'dwi_modb_res_x.nrrd'])
        call(['crlConstructVectorImage'] + outFiles + [inputDir + 'dwi_modb_res_'+axis+'.nrrd'])