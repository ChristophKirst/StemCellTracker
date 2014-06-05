'''Illastic Classification Interface'''

import numpy as np
import sys, os

from IlastikConfig import ilastikpath

print "ilastik path: " + ilastikpath

sys.path.insert(1, ilastikpath)
#sys.path.insert(1, ilastikpath + '/lib')
#sys.path.insert(1, ilastikpath + '/lib/python2.7/site-packages')


# Import vigra
try:
   import vigra
except ImportError, vigraImport:
   print """vigra import: failed to import the vigra library. Please follow the instructions on "http://hci.iwr.uni-heidelberg.de/vigra/" to install vigra"""
   raise vigraImport

print "viagra version: " + vigra.version


# Import h5py
try:
   import h5py
except ImportError, h5pyImport:
   print """h5py import: failed to import the h5py library."""
   raise h5pyImport

print "h5py version: " + h5py.version.version

   
# Import ilastik 

#old_stdout = sys.stdout
try:
   #sys.stdout = sys.stderr = open(os.devnull, "w")
   from ilastik.core.dataMgr import DataMgr, DataItemImage
   from ilastik.modules.classification.core.featureMgr import FeatureMgr
   from ilastik.modules.classification.core.classificationMgr import ClassificationMgr
   from ilastik.modules.classification.core.features.featureBase import FeatureBase
   from ilastik.modules.classification.core.classifiers.classifierRandomForest import ClassifierRandomForest
   from ilastik.modules.classification.core.classificationMgr import ClassifierPredictThread
   from ilastik.core.volume import DataAccessor
   #sys.stdout = old_stdout
except ImportError, ilastikImport:
   #sys.stdout = old_stdout
   print """ilastik import: failed to import the ilastik library. Please follow the instructions on "http://www.ilastik.org" to install ilastik"""
   raise ilastikImport


class IlastikClassifier():
   
   def __init__(self):
      
      #self.image_name = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/hESCells_DAPI.tif'
      #self.classifier_name = '/home/ckirst/Desktop/' + 'classifier.h5' 
      
      self.classifiers = [];
      self.features = [];
      self.image = None;


   def load_image(self, filename = None):
      #if filename is None:
      #    filename = self.image_name
      
      # get input image
      self.image_name = filename
      self.set_image(vigra.readImage(filename))
      #if self.image is None:
      #   self.image_name = None
    
   def set_image(self, image):
      # Transform input image to ilastik conventions
      # 3d = (time,x,y,z,channel) 
      # 2d = (time,1,x,y,channel)
      self.image = image
      
      if len(self.image.shape) == 3:
         self.image.shape = (1,1) + self.image.shape
      elif len(self.image_.shape) == 4:
         self.image.shape = (1,) + self.image.shape
      else:
        print """Error: set image: failed: image must be 2 or 3d."""
        self.image = None
 
      print "image info: shape: " + str(self.image.shape) + " file: " + self.image_name

   def load_randomforest(self, filename = None):
      #if filename is None:
      #   filename = self.classifier_name
         
      prefix = 'classifiers'
      hf = h5py.File(filename,'r')
      if prefix in hf:
         cids = hf[prefix].keys()
         hf.close()
         del hf
         
         classifiers = []
         for cid in cids:
               classifiers.append(ClassifierRandomForest.loadRFfromFile(filename, str(prefix + '/' + cid)))   
         self.classifiers = classifiers
      else:
         raise ImportError('No Classifiers in prefix')

   
   def load_features(self, filename = None):
      #if filename is None:
      #   filename = self.classifier_name
      
      featureItems = []
      hf = h5py.File(filename,'r')
      for fgrp in hf['features'].values():
            featureItems.append(FeatureBase.deserialize(fgrp))
      hf.close()
      del hf
      self.features = featureItems
             
             
   def load_classifier(self, filename = None):
      #if filename is None:
      #   filename = self.classifier_name
      
      self.load_randomforest(filename)
      self.load_features(filename)
  
        
   def run(self, imagename = None):
      if imagename is not None:
         self.load_image(imagename)
         
      if self.image is None:
         raise ImportError('Image cannot be loaded')
      
      # Create ilastik dataMgr
      dataMgr = DataMgr()

      di = DataItemImage('')
      di.setDataVol(DataAccessor(self.image))
      dataMgr.append(di, alreadyLoaded=True)
      dataMgr.module["Classification"]["classificationMgr"].classifiers = self.classifiers

      # Create FeatureMgr
      fm = FeatureMgr(dataMgr, self.features)
      fm.prepareCompute(dataMgr)
      fm.triggerCompute()
      fm.joinCompute(dataMgr)
      
      # Predict with loaded classifier
      classificationPredict = ClassifierPredictThread(dataMgr)
      classificationPredict.start()
      classificationPredict.wait()

      #del dataMgr

      return classificationPredict._prediction[0]






