#!/usr/bin/env python3
'''Example script for a dataset provider to submit metadata and data location records to the SOLARNET Virtual Observatory (SVO) RESTful API'''
#import pdb
import os
import sys
import argparse
import logging
from datetime import datetime
from urllib.parse import urljoin
from dateutil.parser import parse, ParserError
from slumber import API
from astropy.io import fits

# Set with proper base URL of data, or set file URL explicitly through script argument
BASE_FILE_URL = None

# The default hdu name or index to use for extracting the metadata from the FITS file (can be specified here to avoid passing it by parameter to the script)
DEFAULT_FITS_HDU = 0

# If data should be offline by default (can be specified here to avoid passing it by parameter to the script)
DEFAULT_OFFLINE = False

# Default dataset name to use (can be specified here to avoid passing it by parameter to the script)
DEFAULT_DATASET = None

# Default username and API key of the user in the SVO owning the data (can be specified here to avoid passing it by parameter to the script)
DEFAULT_USERNAME = None
DEFAULT_API_KEY = None

# Keyword to use to generate a default oid
DATE_KEYWORD = 'date_obs'


# URL of the SVO RESTful API
# Don't change this
API_URL = 'https://solarnet.oma.be/service/api/svo'

class SvoApi(API):
  '''RESTful API interface for the SVO'''
  def __init__(self, api_url, username, api_key):
    self.username = username
    self.api_key = api_key
    super().__init__(api_url, auth=self.api_key_auth)
  
  def api_key_auth(self, request):
    '''Sets the API key authentication in the request header'''
    request.headers['Authorization'] = 'ApiKey %s:%s' % (self.username, self.api_key)
    return request
  
  def __call__(self, resource_uri):
    '''Return a resource from a resource URI'''
    return getattr(self, resource_uri)

class Record:
  def __init__(self, fits_file, fits_hdu = DEFAULT_FITS_HDU, file_url = None, file_path = None, thumbnail_url = None, offline = DEFAULT_OFFLINE, oid = None, dataset = DEFAULT_DATASET, username = DEFAULT_USERNAME, api_key = DEFAULT_API_KEY, **kwargs):
    self.fits_file = fits_file
    self.fits_hdu = fits_hdu
    self.file_url = file_url
    self.file_path = file_path
    self.thumbnail_url = thumbnail_url
    self.offline = offline
    self.dataset = dataset
    self.oid = oid
    self.api = SvoApi(API_URL, username, api_key)
  
  def get_file_url(self):
    '''Override to return the proper URL for the file'''
    if self.file_url:
      return self.file_url
    elif BASE_FILE_URL:
      return urljoin(BASE_FILE_URL, self.get_file_path())
    else:
      raise ValueError('file_url must be provided or BASE_FILE_URL must be set')
  
  def get_file_path(self):
    '''Override to return the proper relative file path for the file'''
    if self.file_path:
      return self.file_path
    else:
      return self.fits_file.lstrip('./')
  
  def get_thumbnail_url(self):
    '''Override to return the proper URL for the thumbnail'''
    return self.thumbnail_url
  
  def get_oid(self, metadata = None):
    '''Override to return the proper OID for the metadata'''
    if self.oid:
      return self.oid
    elif not metadata:
      raise ValueError('oid or metadata must be provided')
    else:
      try:
        return parse(metadata[DATE_KEYWORD]).strftime('%Y%m%d%H%M%S')
      except ParserError as why:
        raise ValueError('Cannot parse keyword "%s" into an oid: %s' % (DATE_KEYWORD, why)) from why
      except KeyError as why:
        raise ValueError('Keyword "%s" missing in FITS header, cannot generate default oid' % DATE_KEYWORD) from why
  
  def get_data_location(self):
    '''Return the data_location URI or info for creating a new metadata record'''
    # If a data_location record for this file already exist, we reuse it
    # Else we need to create a new one, with all the necessary info
    
    # The pair dataset/file_url must be unique in the database so we can search by it
    data_locations = self.api.data_location.get(dataset__name = self.dataset, file_url = self.get_file_url())
    
    if data_locations['objects']:
      data_location = data_locations['objects'][0]['resource_uri']
    else:
      dataset = self.api.dataset(self.dataset).get()
      
      data_location = {
        'dataset': dataset['resource_uri'],
        'file_url': self.get_file_url(),
        'file_size': os.path.getsize(self.fits_file),
        'file_path': self.get_file_path(),
        'thumbnail_url': self.get_thumbnail_url(),
        'offline': self.offline,
      }
    
    return data_location
  
  def get_metadata(self):
    '''Return the metadata info for creating a new metadata record'''
    
    with fits.open(self.fits_file) as hdus:
      fits_header = hdus[self.fits_hdu].header
    
    # Create a new metadata record from the FITS file header
    metadata = {
      'fits_header': fits_header.tostring().strip()
    }
    
    for keyword in self.api.keyword.get(limit=0, dataset__name=self.dataset)['objects']:
      try:
        metadata[keyword['name']] = fits_header[keyword['verbose_name']]
      except KeyError:
        logging.warning('Missing keyword "%s" in header, skipping!', keyword['verbose_name'])
      else:
        logging.debug('Header keyword "%s" = "%s"', keyword['verbose_name'], metadata[keyword['name']])
    
    metadata['oid'] = self.get_oid(metadata)
    
    return metadata
  
  def create(self):
    '''Create the metadata and data_location records in the API'''
    #pdb.set_trace()
    # Retrieve the metadata resource URI from the dataset
    dataset = self.api.dataset(self.dataset).get()
    resource_uri = dataset['metadata']['resource_uri']
    
    # Get the data to send to the API
    metadata = self.get_metadata()
    metadata['data_location'] = self.get_data_location()
    
    # To create a new record, POST to the resource URI of the metadata record
    result = self.api(resource_uri).post(metadata)
    return result



if __name__ == "__main__":

  # Get the arguments
  parser = argparse.ArgumentParser(description='Submit metadata from a FITS file to the SVO')
  parser.add_argument('--debug', '-d', action='store_true', help='Set the logging level to debug')
  parser.add_argument('fits_file', metavar = 'FITS FILE', help='The FITS file to submit to the SVO')
  parser.add_argument('--fits-hdu', metavar = 'FITS HDU', default=argparse.SUPPRESS, help='The FITS HDU index or name from which to extract the metadata to submit to the SVO')
  parser.add_argument('--file-url', metavar = 'FILE URL', default=argparse.SUPPRESS, help='The URL of the file')
  parser.add_argument('--file-path', metavar = 'FILE PATH', default=argparse.SUPPRESS, help='The relative path of the file')
  parser.add_argument('--thumbnail-url', metavar = 'THUMBNAIL URL', default=argparse.SUPPRESS, help='The URL of the thumbnail')
  offline_parser = parser.add_mutually_exclusive_group(required=False)
  offline_parser.add_argument('--offline', dest='offline', default=argparse.SUPPRESS, action='store_true', help='Set the record as offline')
  offline_parser.add_argument('--online', dest='offline', default=argparse.SUPPRESS, action='store_false', help='Set the record as online (or not offline)')
  parser.add_argument('--dataset', default=argparse.SUPPRESS, help='The name of the dataset in the SVO')
  parser.add_argument('--oid', default=argparse.SUPPRESS, help='The unique observation ID of the metadata')
  parser.add_argument('--username', '-u', default=argparse.SUPPRESS, help='The username (email) of the user owning the data')
  parser.add_argument('--api-key', '-k', default=argparse.SUPPRESS, help='The API key of the user owning the data')
  
  args = parser.parse_args()
  
  # Setup the logging
  if args.debug:
    logging.basicConfig(level = logging.DEBUG, format = '%(levelname)-8s: %(funcName)s %(message)s')
  else:
    logging.basicConfig(level = logging.INFO, format = '%(levelname)-8s: %(message)s')
  
  try:
    record = Record(**vars(args))
    result = record.create()
  except Exception as why:
    logging.critical('Could not submit data: %s', why)
    if getattr(why, 'content', False):
      logging.error(why.content)
    sys.exit(1)
  else:
    logging.info('Succesfully submited data')
    logging.debug(result)

