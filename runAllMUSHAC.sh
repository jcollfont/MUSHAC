#! /bin/bash

# source python
source ~/anaconda2_startup.sh

# reformat B-values
python reformatBvaluesBvecs.py -d $2/MUSHAC/ -s $1

# diamondNCcyl sumToOne=0
# python processMUSHAC_crlDCIEstimate_v2.py -d $2/MUSHAC/ -s $1 -denoise gibbs -diamond diamondNCcyl -sumToOne 0 -p 200
# python saveFSLresults.py -d $2/MUSHAC/ -s $1 -denoise gibbs -diamond diamondNCcyl -sumToOne 0 -o $2/MUSHAC_results/MUSHAC_results_denoised_GIBBS_diamondNCcyl_sumToOne0/

# diamondNCcyl sumToOne=1
python processMUSHAC_crlDCIEstimate_v2.py -d $2/MUSHAC/ -s $1 -denoise gibbs -diamond diamondNCcyl -sumToOne 1 -p 200
python saveFSLresults.py -d $2/MUSHAC/ -s $1 -denoise gibbs -diamond diamondNCcyl -sumToOne 1 -o $2/MUSHAC_results/MUSHAC_results_denoised_GIBBS_diamondNCcyl_sumToOne1/

# diamondcyl sumToOne=0
# python processMUSHAC_crlDCIEstimate_v2.py -d $2/MUSHAC/ -s $1 -denoise gibbs -diamond diamondcyl -sumToOne 0 -p 200
# python saveFSLresults.py -d $2/MUSHAC/ -s $1 -denoise gibbs -diamond diamondcyl -sumToOne 0 -o $2/MUSHAC_results/MUSHAC_results_denoised_GIBBS_diamondcyl_sumToOne0/

# # diamondcyl sumToOne=1
# python processMUSHAC_crlDCIEstimate_v2.py -d $2/MUSHAC/ -s $1 -denoise gibbs -diamond diamondcyl -sumToOne 1 -p 200
# python saveFSLresults.py -d $2/MUSHAC/ -s $1 -denoise gibbs -diamond diamondcyl -sumToOne 1 -o $2/MUSHAC_results/MUSHAC_results_denoised_GIBBS_diamondcyl_sumToOne1/
