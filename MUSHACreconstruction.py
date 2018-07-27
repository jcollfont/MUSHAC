### IMPORTS
# standard python
import numpy as np 
import nrrd
import os
import shutil
from subprocess import call
import sys

# DIPY
from dipy.core.gradients import gradient_table
import dipy.reconst.dki as dki
from dipy.reconst.dsi import DiffusionSpectrumModel

# local
sys.path.insert(0, os.path.expanduser("~") + '/Documents/Research/DWI/MFMestimation/python/')
from loadDWIdata import saveNRRDwithHeader, loadDWIdata, dwiClass


# GLOBAL
refHeader = '/fileserver/projects6/jcollfont/refHeader.nhdr'

class MUSHACreconstruction():

    # ------------------ MEMBERS ----------------------- #
    # file paths
    paths = None
    loadedFiles = None

    # inputs
    imgSize = (0,0,0)
    anatMask = np.arange(0)
    numTensors = 0
    diamondTensor = None
    diamondFractions = np.zeros([0])
    diamondMose = np.zeros([0])
    diamondB0 = np.zeros([0])
    diamondFreeWaterDiff = 0

    # simulation
    dwiData = np.zeros([0])
    gtab = []

    # ------------------ METHODS ----------------------- #
    #
    #       - __init__                  -> constructor
    #       - __setTensorMatrixfrom6D   -> create full tensor matrices from 6D vectors (upper triang)
    #       - generateDWIdata           -> generate DWI data from DIAMOND model
    #       - computeDIAMONDTensorParams -> compute scalar params from tensor model
    #       - computeParamsFromSingleTensorFromDWI -> compute scalar params from tensor model (MUSHAC)
    #
    #

    #%% constructor
    def __init__(self, basePath='', path2Diamond='', path2DWI='', refName='', maskPath='', diamondFreeWaterDiff=3e-3):

        # set path
        self.paths = {'base':'', 'diamond':'','dwi':''}
        self.paths['base'] = basePath
        self.paths['dwi'] = self.paths['base'] + path2DWI
        self.paths['diamond'] = self.paths['base'] + path2Diamond
        files = os.listdir( self.paths['diamond'] )

        # get MOSE mask
        self.loadedFiles = {'mask':'', 'b0':'', 'mose':'', 'fractions':'', 'tensors':[]}
        self.loadedFiles['mose'] = self.paths['diamond'] + [ff for ff in files \
                                                                if ('_mosemap.nrrd' in ff)&( refName in ff )&( 'mtm' not in ff )][0]
        self.diamondMose = nrrd.read( self.loadedFiles['mose'] )[0]
        self.numTensors = np.max(self.diamondMose.ravel())
        self.imgSize = self.diamondMose.shape
        print '\t-MOSE: %s' %(self.loadedFiles['mose'])

        # get anatomical mask
        self.loadedFiles['mask'] = self.paths['base'] + maskPath
        try:
            mask = nrrd.read( self.loadedFiles['mask'] )[0]
            self.anatMask = np.where(mask.ravel())[0]
            print '\t-Mask: %s' %(self.loadedFiles['mask'])
        except:
            self.anatMask = np.arange(np.prod(self.imgSize))
            print '\t-Mask: No mask found! Using all voxels'

        # reshape mose to mask format
        self.diamondMose = self.diamondMose.ravel()[self.anatMask]

        # get fractions
        self.loadedFiles['fractions'] = self.paths['diamond'] + [ff for ff in files \
                                                                    if ('_fractions.nrrd' in ff)&( refName in ff )&( 'mtm' not in ff )][0]
        self.diamondFractions = nrrd.read( self.loadedFiles['fractions'] )[0]
        self.diamondFractions = self.diamondFractions[:self.numTensors+1,:,:,:]\
                                                .reshape(self.numTensors+1, np.prod(self.imgSize)).T[self.anatMask,:]
        print '\t-Fractions: %s' %(self.loadedFiles['fractions'])

        # get B0
        self.loadedFiles['b0'] = self.paths['diamond'] + [ff for ff in files \
                                                            if ('_b0.nrrd' in ff)&( refName in ff )&( 'mtm' not in ff )][0]
        self.diamondB0 = nrrd.read( self.loadedFiles['b0'] )[0].ravel()[self.anatMask]
        print '\t-Bo: %s' %(self.loadedFiles['b0'])

        # get tensors
        self.diamondTensor = range(self.numTensors)
        tensorList = list()
        for tt in range(self.numTensors):
            
            tensorList +=  [self.paths['diamond'] + [ff for ff in files \
                                                                        if ('_t%d.nrrd' %(tt) in ff)&( refName in ff )&( 'mtm' not in ff )][0]]

            tensor6D = nrrd.read( tensorList[-1] )[0].reshape( (6, np.prod(self.imgSize)) ).T[self.anatMask,:]
            self.__setTensorMatrixfrom6D( tensor6D, tt )

            print '\t-Tensor %d: %s' %(tt, tensorList[-1])

        self.loadedFiles['tensors'] = tensorList
        print self.loadedFiles['tensors'] 

        # set water fraction
        self.diamondFreeWaterDiff = diamondFreeWaterDiff

    #%% redo tensor shape
    def __setTensorMatrixfrom6D(self, tens6D, tt):

        self.diamondTensor[tt] = np.zeros( ( tens6D.shape[0], 3,3 ))
        self.diamondTensor[tt][:,0,:] = tens6D[:,:3]
        self.diamondTensor[tt][:,1,1:] = tens6D[:,3:5]
        self.diamondTensor[tt][:,2,2] = tens6D[:,5]
        self.diamondTensor[tt][:,1,0] = tens6D[:,1]
        self.diamondTensor[tt][:,2,0] = tens6D[:,2]
        self.diamondTensor[tt][:,2,1] = tens6D[:,4]

    #%% 
    #
    #   This function generates new DWI data from the DIAMOND model
    #
    #   INPUTS:
    #       - bvecs - <G,3>double - normalized gradient vectors
    #       - bvals - <G,1>double - B-values that match gradient vetors
    #       - outputPath - str - path where the new DWI data will be saved.
    #                              (if '' is provided, the algorithm only returns the raw data)
    #       
    #
    #   OUTPUT:
    #       - reconstructedDWI - <imgSize,G> double - signal generated from the diamond model.
    #                               (if outputPath != '')
    #                           - dwiClass - dwi class object containing the signal generated 
    #                                       from the diamond model. (if outputPath == '') 
    #
    def generateDWIdata(self, bvecs, bvals, outputPath=''):

        # set gradient table
        self.gtab = gradient_table(bvals, bvecs)

        # adjust S0 according to TE and TR
        S0 = self.diamondB0

        # prealocate
        numGrad = self.gtab.bvals.size
        signal = np.zeros([np.prod(self.imgSize), numGrad])
        # for every bvec and bval
        for bb in range( numGrad ):

            # for every tensor generate decaying exp
            expData = np.zeros([ self.anatMask.size, self.numTensors +1 ])
            for tt in range(self.numTensors):
                expData[:,tt] = np.exp( -1.0*self.gtab.bvals[bb] \
                                * self.gtab.bvecs[bb,:].dot( self.diamondTensor[tt] ).dot( self.gtab.bvecs[bb,:].T ) )

            # generate decaying exp for water fraction
            expData[:,-1] = np.tile( np.exp( -1.0*self.gtab.bvals[bb] * self.diamondFreeWaterDiff \
                                                * self.gtab.bvecs[bb,:].dot( self.gtab.bvecs[bb,:].T ) ), \
                                                [ self.anatMask.size ])

            # generate data
            signal[self.anatMask,bb] = np.sum( self.diamondFractions * expData ,axis=1) * S0
        

        if outputPath == '':
            reconstructedDWI = signal.reshape( self.imgSize + (numGrad,) )
        else:
            saveNRRDwithHeader( signal.reshape( self.imgSize + (numGrad,) ), refHeader, \
                            outputPath, os.path.basename(self.paths['dwi'])[:-5], \
                            bvals, bvecs )
                            
            reconstructedDWI = dwiClass( outputPath + os.path.basename(self.paths['dwi']), \
                                    maskPath=self.loadedFiles['mask'], outputPath=outputPath)
            
        
        return reconstructedDWI


    #%% 
    #
    #   This function computes the basic parameter estimates from the tensors in DIAMOND
    #
    def computeDIAMONDTensorParams(self, outputDir='', outputName=''):

        if outputDir == '':
            outputDir = '/tmp/tmp%d/'  %(np.random.randint(1e6))
            try:
                os.makedirs(outputDir)
            except:
                print outputDir + ' already exists'

        # for all tensors
        diffusionVec = []
        sticks = []
        meanDiffusivity = []
        fractionalAnisotropy = []
        radialDiffusivity = []
        axialDiffusivity = []
        for tt in range(self.numTensors):
            
            # compue sticks
            call([ 'crlDCIToPeaks', '-i', self.loadedFiles['tensors'][tt], '-n 1', \
                                        '-o',  outputDir + outputName + 'tensStick%d.nrrd'%(tt) ])

            # compute params
            call(['crlTensorScalarParameter', '-m', outputDir + outputName + 'tensMD%d.nrrd'%(tt), \
                                                '-f', outputDir + outputName + 'tensFA%d.nrrd'%(tt), \
                                                '-1', outputDir + outputName + 'tensEIG0%d.nrrd'%(tt), \
                                                '-2', outputDir + outputName + 'tensEIG1%d.nrrd'%(tt), \
                                                '-3', outputDir + outputName + 'tensEIG2%d.nrrd'%(tt), \
                                                '-r', outputDir + outputName + 'tensRD%d.nrrd'%(tt), \
                                                '-a', outputDir + outputName + 'tensAD%d.nrrd'%(tt), \
                                            self.loadedFiles['tensors'][tt] ])

            # load diffusion
            eig0 = nrrd.read( outputDir + outputName + 'tensEIG0%d.nrrd'%(tt) )[0]
            eig1 = nrrd.read( outputDir + outputName + 'tensEIG1%d.nrrd'%(tt) )[0]
            eig2 = nrrd.read( outputDir + outputName + 'tensEIG2%d.nrrd'%(tt) )[0]
            diffusionVec.append( np.concatenate( [eig0, eig1, eig2], axis=-1).reshape( np.prod(self.imgSize), 3 )  )
            # load stick
            sticks.append( nrrd.read(outputDir + outputName + 'tensStick%d.nrrd'%(tt))[0].reshape( np.prod(self.imgSize), 3 ) )
            # load MD
            meanDiffusivity.append(  nrrd.read(outputDir + outputName + 'tensMD%d.nrrd'%(tt))[0].ravel()  ) 
            #load FA
            fractionalAnisotropy.append(  nrrd.read(outputDir + outputName + 'tensFA%d.nrrd'%(tt))[0].ravel() ) 
            # load RD
            radialDiffusivity.append(  nrrd.read(outputDir + outputName + 'tensRD%d.nrrd'%(tt))[0].ravel()  ) 
            # load AD
            axialDiffusivity.append(  nrrd.read(outputDir + outputName + 'tensAD%d.nrrd'%(tt))[0].ravel()  ) 

        return sticks, diffusionVec, meanDiffusivity, fractionalAnisotropy, radialDiffusivity, axialDiffusivity

        






    #%% 
    def computeParamsFromSingleTensorFromDWI(self, bvecs=[], bvals=[], recSignal=np.zeros([0]), recNHDR='', outputName=None, anatMask=None):

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
        if recNHDR== '':
            saveNRRDwithHeader( recSignal, refHeader, tmpdir, '/recDWI' , bvals, bvecs )
            recNHDR = tmpdir+'/recDWI.nhdr'

        # compute single tensor on data
        if not os.path.exists(tmpdir + '/recDTI.nrrd'):
            call(['tend', 'estim', '-i', recNHDR, '-o', tmpdir + '/' + baseName + '_recDTI.nrrd' ,  '-B', 'kvp', '-knownB0', 'false'  ])

        # compute FA, MD
        if not os.path.exists(tmpdir +'/'+ baseName + '_fracAnis.nrrd'):
            call(['crlTensorScalarParameter', '-m',  tmpdir +'/'+ baseName + '_meanDiff.nrrd', '-f',  tmpdir +'/'+ baseName + '_fracAnis.nrrd', tmpdir + '/recDTI.nrrd'])

        # load results
        meanDiff = nrrd.read( tmpdir +'/'+  baseName + '_meanDiff.nrrd' )[0]
        fracAnisotropy = nrrd.read( tmpdir +'/'+ baseName + '_fracAnis.nrrd' )[0]


        # load DWI data
        dwiData = loadDWIdata( os.path.dirname(recNHDR) + '/', os.path.basename(recNHDR) )[0]

        if not anatMask == None:
            mask = nrrd.read(anatMask)[0]
            dwiData = dwiData* np.tile( mask , [dwiData.shape[-1],1,1,1]).transpose(1,2,3,0)

        # compute Mean Kurtosis (MK)
        if not os.path.exists(tmpdir +'/'+ baseName + '_meanKurtosis.nrrd'):
            dkiModel = dki.DiffusionKurtosisModel(self.gtab)
            dkifit = dkiModel.fit(dwiData)
            meanKurtosis = dkifit.mk()
            nrrd.write( tmpdir +'/'+ baseName + '_meanKurtosis.nrrd', meanKurtosis )
        else:
            meanKurtosis = nrrd.read( tmpdir +'/'+ baseName + '_meanKurtosis.nrrd')[0]

        # Compute Return to Origin Probability (RTOP)
        # if not os.path.exists(tmpdir +'/'+ baseName + '_rtop.nrrd'):
        #     dsmodel = DiffusionSpectrumModel( self.gtab )
        #     rtop = dsmodel.fit(  dwiData / np.mean( dwiData[:,:,:, self.gtab.b0s_mask ] ,axis=3, keepdims=True)  ).rtop_pdf()
        #     nrrd.write( tmpdir +'/'+ baseName + '_rtop.nrrd', rtop )
        # else:
        #     rtop = nrrd.read( tmpdir +'/'+ baseName + '_rtop.nrrd')[0]


        if outputName == None:
            shutil.rmtree(tmpdir)


        return meanDiff, fracAnisotropy, meanKurtosis#, rtop

#%% check whether to compute file
def __fileComputeCheck(self, newFile, oldFile='./'):
    computeBool = False
    if not os.path.exists( newFile ):
            computeBool = True
    elif os.path.getmtime(oldFile) > os.path.getmtime(newFile):
            computeBool = True     
    return computeBool