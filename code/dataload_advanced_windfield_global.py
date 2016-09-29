#### function to load hurricane data
#### USE with programme advanced_windfield_global.pyplot

from __future__ import division
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import csv
import os, glob
from datetime import datetime
from time import time
from netCDF4 import Dataset
import dimarray as da

#### load data on latitude, longitude, pressure, windspeed, radius of max-winds, poci (pressure of outer-closed isobar)
#### load IBTrACS data
def TC_load_IBTrACS(infile,lonmin,lonmax,hemisphere,basin,provider, provider_alt,provider_alt2):
	# yea = str(year)
	
	stop=False
	
	ncfile=Dataset(infile)
	stname=infile[-31:-18]
	if stname[7]=='N':
		placeholder='0'
	else:
		placeholder='1'
	stmid=stname
	stmid=int(stname[:7]+placeholder+stname[8:])
	# preselect all tracks by basin       
	if (stname[7]==hemisphere):
		if basin != 'SP':
			lonorig=int(stname[-3:])
			if lonorig > 180:
				lonorig=lonorig -360
			if (lonorig >= lonmin/10. ) & (lonorig <= lonmax/10. ):
				print stname, 'passed'
				pass
			else:
				stop=True
		else:
			lonorig=int(stname[-3:])
			if (lonorig >= 135 ) & (lonorig <= 290 ):
				print stname, 'passed'
				pass
			else:
				stop=True
	else:
		stop=True
	
	# load variables

	try:
		latsall=ncfile.variables[provider+'_lat'][:]
		lonsall=ncfile.variables[provider+'_lon'][:]
		# 1-min sustained wind
		windall=ncfile.variables[provider+'_wind'][:].astype('int')
		presall=ncfile.variables[provider+'_pres'][:].astype('int')
	except(KeyError):
		try:
			print 'provider {0} used for basin {1}'.format(provider_alt, basin)
			latsall=ncfile.variables[provider_alt+'_lat'][:]
			lonsall=ncfile.variables[provider_alt+'_lon'][:]
			# 1-min sustained wind
			windall=ncfile.variables[provider_alt+'_wind'][:].astype('int')
			presall=ncfile.variables[provider_alt+'_pres'][:].astype('int')
		except(KeyError):
			if provider_alt2 != '':
				try:
					print 'provider {0} used for basin {1}'.format(provider_alt2, basin)
					latsall=ncfile.variables[provider_alt2+'_lat'][:]
					lonsall=ncfile.variables[provider_alt2+'_lon'][:]
					# 1-min sustained wind
					windall=ncfile.variables[provider_alt2+'_wind'][:].astype('int')
					presall=ncfile.variables[provider_alt2+'_pres'][:].astype('int')
				except(KeyError):
					print 'data ..._for_mapping used for basin {1}'.format(provider_alt2, basin)
					latsall=ncfile.variables['lat_for_mapping'][:]
					lonsall=ncfile.variables['lon_for_mapping'][:]
					# 10 min-sustained winds (average over many providers) # not very crucial as few storms and mostly one provider
					windall=ncfile.variables['wind_for_mapping'][:].astype('int')
					presall=ncfile.variables['pres_for_mapping'][:].astype('int')
			else:
				print 'data ..._for_mapping used for basin {1}'.format(provider_alt2, basin)
				latsall=ncfile.variables['lat_for_mapping'][:]
				lonsall=ncfile.variables['lon_for_mapping'][:]
				# 10 min-sustained winds (average over many providers) # not very crucial as few storms and mostly one provider
				windall=ncfile.variables['wind_for_mapping'][:].astype('int')
				presall=ncfile.variables['pres_for_mapping'][:].astype('int')
	lonsall=(lonsall*10).astype('int')
	latsall=(latsall*10).astype('int')
	if basin !='EP':
		try:
			try:
				rmaxall=ncfile.variables[provider+'_rmw'][:].astype('int')
				poci=ncfile.variables[provider+'_poci'][:].astype('int')
			except (KeyError):
				print 'alternative provider needed'
				rmaxall=ncfile.variables[provider_alt+'_rmw'][:].astype('int')
				poci=ncfile.variables[provider_alt+'_poci'][:].astype('int')
		except (KeyError):        
			rmaxall=np.zeros(len(latsall))
			poci=np.zeros(len(latsall))
	else:
		try:
			try:
				rmaxall=ncfile.variables[provider_alt+'_rmw'][:].astype('int')
				poci=ncfile.variables[provider_alt+'_poci'][:].astype('int')
			except (KeyError):
				print 'alternative provider needed'
				rmaxall=ncfile.variables[provider_alt2+'_rmw'][:].astype('int')
				poci=ncfile.variables[provider_alt2+'_poci'][:].astype('int')
		except (KeyError):        
			rmaxall=np.zeros(len(latsall))
			poci=np.zeros(len(latsall))
			
	return stop, stname, stmid, lonsall, latsall, windall, presall, rmaxall, poci
	
	
#### load IBTrACS data
def TC_load_emanuel(infile,lonmin,lonmax,hemisphere,basin,ibeg,stm, gcmno):
	
	stop=False; here=False
	
	latsall=infile['latstore'][ibeg]
	latsall=(latsall*10).astype('int')
	if (hemisphere=='N'):
		latsall=latsall[(np.where(latsall > 0))]
	else:
		latsall=latsall[(np.where(latsall < 0))]
	
	if (hemisphere=='N') & ((latsall > 0).any()):
		pass
	elif (hemisphere=='S') & ((latsall < 0).any()):
		#here=True
		pass
	else:
		stop=True
	
	lonsall=infile['longstore'][ibeg]
	lonsall=(lonsall*10).astype('int')
	if (hemisphere=='N'):
		lonsall=lonsall[np.where(latsall>0)]
	else:
		lonsall=lonsall[np.where(latsall<0)]
	lonsall[lonsall > 1800]=lonsall[lonsall > 1800]-3600
	
	if ((lonsall > 0).any()) & ((lonsall < 0).any()):
		lonspos=lonsall[lonsall >= 0]
		lonsneg=lonsall[lonsall < 0]
		#here=True
		if (((lonspos <= lonmax).any()) & ((lonspos >= lonmin).any())) & (((lonsneg <= lonmax).any()) & ((lonsneg >= lonmin).any())):
			pass
		else:
			stop=True
	else: 
		#print ((lonsall >= lonmin).any())
		#print np.min(lonsall) , np.max(lonsall)
		if ((lonsall >= lonmin).any()) & ((lonsall <= lonmax).any()):
			here=True
			pass
		else:
			stop=True
	
	windall=infile['vstore'][ibeg].astype('int')
	if (hemisphere=='N'):
		windall=windall[np.where(latsall>0)]
	else:
		windall=windall[np.where(latsall<0)]
	
	rmaxall=infile['rmstore'][ibeg].astype('int')
	if (hemisphere=='N'):
		rmaxall=rmaxall[np.where(latsall>0)]
	else:
		rmaxall=rmaxall[np.where(latsall<0)]
	# transform to nautical miles
	rmaxall=(rmaxall / 1.852 ).astype('int')
	
	presall=infile['pstore'][ibeg].astype('int')
	if (hemisphere=='N'):
		presall=presall[np.where(latsall>0)]
	else:
		presall=presall[np.where(latsall<0)]
	
	poci=np.zeros(len(latsall))
	
	#tempf=mat['freqyear'].flatten()
	#tempf=tempf.astype('int')
	#tempf=tempf[no]
	tempd=infile['daystore'][ibeg][0]
	tempm=infile['monthstore'][ibeg][0]
	tempy=infile['yearstore'].flatten()[ibeg]
	
	# define stmid==stname
	if tempd < 10:
		tempd='0'+str(tempd)
	if tempm < 10:
		tempm='0'+str(tempm)
	if stm < 100:
		if stm < 10:
			stm='00' + str(stm)
		else:
			stm='0' + str(stm)

	stmid=str(tempy) + str(tempm) + str(tempd) +str(gcmno) +str(stm)
	stname=stmid
	
	return here,stop, stname, stmid, lonsall, latsall, windall, presall, rmaxall, poci
