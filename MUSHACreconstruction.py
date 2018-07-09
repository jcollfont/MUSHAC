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

# local
sys.path.insert(0, '/home/ch199899/Documents/Research/DWI/MFMestimation/python/')
from loadDWIdata import saveNRRDwithHeader


# GLOBAL
refHeader = '/fileserver/projects6/jcollfont/refHeader.nhdr'

class MUSHACreconstruction():

    # ------------------ INSTANCES ----------------------- #
    # file paths
    paths = {'base':'', 'diamond':'','dwi':''}
    loadedFiles = {'mask':'', 'b0':'', 'mose':'', 'fractions':'', 'tensors':[]}

    # inputs
    imgSize = [0,0,0]
    anatMask = np.arange(0)
    numTensors = 0
    diamondTensor = []
    diamondFractions = np.zeros([0])
    diamondMose = np.zeros([0])
    diamondB0 = np.zeros([0])
    diamondFreeWaterDiff = 0

    # simulation
    dwiData = np.zeros([0])
    gtab = []

    # ------------------ METHODS ----------------------- #


    #%% constructor
    def __init__(self, basePath='', path2Diamond='', path2DWI='', refName='', maskPath='', diamondFreeWaterDiff=3e-3):

        # set path
        self.paths['base'] = basePath
        self.paths['dwi'] = self.paths['base'] + path2DWI
        self.paths['diamond'] = self.paths['base'] + path2Diamond
        files = os.listdir( self.paths['diamond'] )
        print 'From path %s load:' %(self.paths['diamond'])

        # get MOSE mask
        self.loadedFiles['mose'] = self.paths['diamond'] + [ff for ff in files if ('_mosemap.nrrd' in ff)&( refName in ff )&( 'mtm' not in ff )][0]
        self.diamondMose = nrrd.read( self.loadedFiles['mose'] )[0]
        self.numTensors = np.max(self.diamondMose.ravel())
        self.imgSize = self.diamondMose.shape
        print '\t-%s' %(self.loadedFiles['mose'])

        # get anatomical mask
        self.loadedFiles['mask'] = self.paths['base'] + maskPath
        try:
            mask = nrrd.read( self.loadedFiles['mask'] )[0]
            self.anatMask = np.where(mask.ravel())[0]
            print '\t-%s' %(self.loadedFiles['mask'])
        except:
            self.anatMask = np.arange(np.prod(self.imgSize))
            print 'No mask found! Using all voxels'

        # reshape mose to mask format
        self.diamondMose = self.diamondMose.ravel()[self.anatMask]

        # get fractions
        self.loadedFiles['fractions'] = self.paths['diamond'] + [ff for ff in files if ('_fractions.nrrd' in ff)&( refName in ff )&( 'mtm' not in ff )][0]
        self.diamondFractions = nrrd.read( self.loadedFiles['fractions'] )[0]
        self.diamondFractions = self.diamondFractions.reshape(self.numTensors+1, np.prod(self.imgSize)).T[self.anatMask,:]
        print '\t-%s' %(self.loadedFiles['fractions'])

        # get B0
        self.loadedFiles['b0'] = self.paths['diamond'] + [ff for ff in files if ('_b0.nrrd' in ff)&( refName in ff )&( 'mtm' not in ff )][0]
        self.diamondB0 = nrrd.read( self.loadedFiles['b0'] )[0].ravel()[self.anatMask]
        print '\t-%s' %(self.loadedFiles['b0'])

        # get tensors
        for tt in range(self.numTensors):
            self.diamondTensor.append( -1 )
            
            self.loadedFiles['tensors'] +=  [self.paths['diamond'] + [ff for ff in files if ('_t%d.nrrd' %(tt) in ff)&( refName in ff )&( 'mtm' not in ff )][0]]
            tensor6D = nrrd.read( self.loadedFiles['tensors'][-1] )[0].reshape( (6, np.prod(self.imgSize)) ).T[self.anatMask,:]
            self.__setTensorMatrixfrom6D( tensor6D, tt )

            print '\t-%s' %(self.loadedFiles['tensors'][-1])

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
    def generateDWIdata(self, bvecs, bvals):

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
                expData[:,tt] = np.exp( -1.0*self.gtab.bvals[bb] * self.gtab.bvecs[bb,:].dot( self.diamondTensor[tt] ).dot( self.gtab.bvecs[bb,:].T ) )

            # generate decaying exp for water fraction
            expData[:,-1] = np.tile( np.exp( -1.0*self.gtab.bvals[bb] * self.diamondFreeWaterDiff * self.gtab.bvecs[bb,:].dot( self.gtab.bvecs[bb,:].T ) ), [ self.anatMask.size ])

            # generate data
            signal[self.anatMask,bb] = np.sum( self.diamondFractions * expData ,axis=1) * S0
        

        return signal.reshape( self.imgSize + (numGrad,) )


    #%% 
    def computeParamsFromSingleTensorFromDWI(self, bvecs=[], bvals=[], recSignal=np.zeros([0]), recNHDR=''):

        # tempSave new signal
        tmpdir = '/tmp/tmp%d'  %(np.random.randint(1e6))
        try:
            os.makedirs(tmpdir)
        except:
            print tmpdir + ' already exists'
        if recNHDR== '':
            saveNRRDwithHeader( recSignal, refHeader, tmpdir, '/recDWI' , bvals, bvecs )
            recNHDR = tmpdir+'/recDWI.nhdr'

        # compute single tensor on data
        call(['tend', 'estim', '-i', recNHDR, '-o', tmpdir + '/recDTI.nrrd' ,  '-B', 'kvp', '-knownB0', 'false'  ])

        # compute FA, MD
        call(['crlTensorScalarParameter', '-m',  tmpdir + '/meanDiff.nrrd', '-f',  tmpdir + '/fracAnis.nrrd', tmpdir + '/recDTI.nrrd'])

        # load results
        meanDiff = nrrd.read( tmpdir + '/meanDiff.nrrd' )[0]
        fracAnisotropy = nrrd.read( tmpdir + '/fracAnis.nrrd' )[0]

        shutil.rmtree(tmpdir)


        return meanDiff, fracAnisotropy