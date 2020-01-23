#!/usr/bin/env python3
import os
import sys
import argparse
import logging
from datetime import datetime
from urllib.parse import urljoin
from slumber import API
from astropy.io import fits
from requests.auth import AuthBase

# Set with proper base directory of data
BASE_FILE_DIRECTORY = '/ABSOLUTE/PATH/TO/DATA/DIRECTORY'

# Set with proper base URL of data
BASE_FILE_URL = 'https://dubshen.astro.su.se/data/'

# Default username and API key of the user in the SOLARNET Data Archive owning the data
DEFAULT_USERNAME = ''
DEFAULT_API_KEY = ''

# Crisp and Chromis data must always be offline to disable download
DEFAULT_OFFLINE = True

# URL of the RESTFull API of the SOLARNET Data Archive
# Don't change this
SOLARNET_API = 'http://solarnet.oma.be/SDA/api/v1'


class TastypieApiKeyAuth(AuthBase):
	'''Sets the appropriate authentication headers for the RESTFull API'''
	def __init__(self, username, api_key):
		self.username = username
		self.api_key = api_key
	
	def __call__(self, request):
		request.headers['Authorization'] = 'ApiKey %s:%s' % (self.username, self.api_key)
		return request

class Record:
	def __init__(self, dataset, fits_file, file_url = None, file_path = None, thumbnail_url = None, oid = None, offline = DEFAULT_OFFLINE, username = DEFAULT_USERNAME, api_key = DEFAULT_API_KEY, **kwargs):
		self.dataset = dataset
		self.fits_file = fits_file
		self.file_url = file_url
		self.file_path = file_path
		self.thumbnail_url = thumbnail_url
		self.oid = oid
		self.offline = offline
		self.api = API(SOLARNET_API, auth = TastypieApiKeyAuth(username, api_key))
		self.HDUs = fits.open(self.fits_file)
		self.fits_header = self.HDUs[0].header
	
	def get_file_url(self):
		'''Override to return the proper URL for the file'''
		if self.file_url:
			return self.file_url
		else:
			return urljoin(BASE_FILE_URL, self.get_file_path())
	
	def get_file_path(self):
		'''Override to return the proper relative file path for the file'''
		if self.file_path:
			return self.file_path
		else:
			return os.path.relpath(self.fits_file, BASE_FILE_DIRECTORY)
	
	def get_thumbnail_url(self):
		'''Override to return the proper URL for the thumbnail'''
		return self.thumbnail_url
	
	def get_oid(self):
		'''Override to return the proper OID for the metadata'''
		if self.oid:
			return self.oid
		else:
			return datetime.strptime(self.fits_header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S').strftime('%Y%m%d%H%M%S')
	
	def get_data_location(self):
		# If a data_location record for this file already exist, we reuse it
		# Else we need to create a new one, with all the necessary info
		
		# The pair dataset/file_url must be unique in the database so we can search by it
		data_locations = self.api.data_location.get(dataset = self.dataset, file_url = self.get_file_url())
		
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
			
			if self.file_path:
				data_location['file_path'] = self.file_path
			else:
				data_location['file_path'] = os.path.basename(self.fits_file)
			
			if self.thumbnail_url:
				data_location['thumbnail_url'] = self.thumbnail_url
		
		return data_location
	
	def get_metadata(self):
		# Create a new metadata record from the FITS file header
		metadata = {
			'oid' : self.get_oid(),
			'fits_header': self.fits_header.tostring()
		}
		
		for keyword in self.api.keyword.get(limit=0, dataset=self.dataset)['objects']:
			try:
				metadata[keyword['db_column']] = self.fits_header[keyword['name']]
			except KeyError:
				logging.warning('Missing keyword "%s" in header, skipping!', keyword['name'])
			else:
				logging.debug('Header keyword "%s" = "%s"', keyword['name'], metadata[keyword['db_column']])
		
		return metadata
	
	def submit(self):
		# Send the metadata and data_location records to the server
		metadata = self.get_metadata()
		metadata['data_location'] = self.get_data_location()
		# Always PUT to the detail URI of the metadata record
		# to allow creating new record or updating existing one
		# The detail URI for the metadata record is metadata/{dataset}/{oid}
		result = self.api.metadata(self.dataset)(metadata['oid']).put(metadata)
		return result.get('objects')



if __name__ == "__main__":

	# Get the arguments
	parser = argparse.ArgumentParser(description='Submit meta-data from a FITS file to the SOLARNET Data Archive')
	parser.add_argument('--debug', '-d', action='store_true', help='Set the logging level to debug')
	parser.add_argument('dataset', help='The dataset ID in the SOLARNET Data Archive')
	parser.add_argument('fits_file', metavar = 'FITS FILE', help='The FITS file to submit to the SOLARNET Data Archive')
	parser.add_argument('--file-url', default=argparse.SUPPRESS, help='The URL of the file')
	parser.add_argument('--file-path', default=argparse.SUPPRESS, help='The relative path of the file')
	parser.add_argument('--thumbnail-url', default=argparse.SUPPRESS, help='The URL of the thumbnail')
	parser.add_argument('--oid', default=argparse.SUPPRESS, help='The unique observation ID of the metadata')
	parser.add_argument('--username', '-u', default=argparse.SUPPRESS, help='The username of the user owning the data')
	parser.add_argument('--api-key', '-k', default=argparse.SUPPRESS, help='The API key of the user owning the data')
	
	args = parser.parse_args()
	
	# Setup the logging
	if args.debug:
		logging.basicConfig(level = logging.DEBUG, format='%(levelname)-8s: %(message)s')
	else:
		logging.basicConfig(level = logging.INFO, format='%(levelname)-8s: %(message)s')
	
	try:
		record = Record(**vars(args))
		result = record.submit()
	except Exception as why:
		logging.critical('Could not submit data: %s', why)
		if getattr(why, 'content', False):
			logging.error(why.content)
		sys.exit(1)
	else:
		logging.info('Succesfully submited data')
		logging.debug(result)
