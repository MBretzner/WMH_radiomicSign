
# coding: utf-8

# # Radiomics tests

# In[1]:


from __future__ import print_function
import sys
import collections
import csv
import os
import logging
import six
import pywt
import SimpleITK as sitk
from radiomics import featureextractor, getFeatureClasses
import radiomics


# In[2]:


# Regulate verbosity with radiomics.setVerbosity
logger = radiomics.logger
logger.setLevel(logging.INFO)

# Write out all log entries to a file
handler = logging.FileHandler(filename='bug_log.txt', mode='w')
formatter = logging.Formatter('%(levelname)s:%(name)s: %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


# Input : CSV batch called "batch.csv". First column named "Image", second column "Mask", rows = paths
# Output : CSV "Radiomics_features.csv"

# In[3]:


outPath = r''

inputCSV = os.path.join(outPath, 'paths.csv')
outputFilepath = os.path.join(outPath, 'MRIgenie_no_coordinates_Resamp_2D_116_FlairIntres_binW5_NAWM.csv')
progress_filename = os.path.join(outPath, 'progress_log.txt')
params = os.path.join(outPath, 'extraction_params.yaml')

# Configure logging
rLogger = logging.getLogger('radiomics')

# Set logging level
# rLogger.setLevel(logging.INFO)  # Not needed, default log level of logger is INFO

# Create handler for writing to log file
handler = logging.FileHandler(filename=progress_filename, mode='w')
handler.setFormatter(logging.Formatter('%(levelname)s:%(name)s: %(message)s'))
rLogger.addHandler(handler)

# Initialize logging for batch log messages
logger = rLogger.getChild('batch')

# Set verbosity level for output to stderr (default level = WARNING)
radiomics.setVerbosity(logging.INFO)

logger.info('pyradiomics version: %s', radiomics.__version__)
logger.info('Loading CSV')

flists = []
try:
  with open(inputCSV, 'r') as inFile:
    cr = csv.DictReader(inFile, lineterminator='\n')
    flists = [row for row in cr]
except Exception:
  logger.error('CSV READ FAILED', exc_info=True)

logger.info('Loading Done')
logger.info('Patients: %d', len(flists))

if os.path.isfile(params):
  extractor = featureextractor.RadiomicsFeatureExtractor(params)

logger.info('Enabled input images types: %s', extractor.enabledImagetypes)
logger.info('Enabled features: %s', extractor.enabledFeatures)
logger.info('Current settings: %s', extractor.settings)

headers = None

for idx, entry in enumerate(flists, start=1):

  logger.info("(%d/%d) Processing Patient (Image: %s, Mask: %s)", idx, len(flists), entry['Image'], entry['Mask'])

  imageFilepath = entry['Image']
  maskFilepath = entry['Mask']
  label = entry.get('Label', None)

  if str(label).isdigit():
    label = int(label)
  else:
    label = None

  if (imageFilepath is not None) and (maskFilepath is not None):
    featureVector = collections.OrderedDict(entry)
    featureVector['Image'] = os.path.basename(imageFilepath)
    featureVector['Mask'] = os.path.basename(maskFilepath)

    try:
      featureVector.update(extractor.execute(imageFilepath, maskFilepath, label))

      with open(outputFilepath, 'a') as outputFile:
        writer = csv.writer(outputFile, lineterminator='\n')
        if headers is None:
          headers = list(featureVector.keys())
          writer.writerow(headers)

        row = []
        for h in headers:
          row.append(featureVector.get(h, "N/A"))
        writer.writerow(row)

    except Exception:
      logger.error('FEATURE EXTRACTION FAILED', exc_info=True)
