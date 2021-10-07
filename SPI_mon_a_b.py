import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime
from numpy import arange, dtype 
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.stats as stats 
from datetime import date, timedelta
from scipy.stats import gamma
import warnings
warnings.filterwarnings("ignore")
#define the corners of the target rectangular
llclat=24.  
llclon=-124. +360.
urclat=50. 
urclon=-67. +360.

yst=1993 #start of SPI calibration (year)
yfin=2010 # end of SPI calibration (year)
#Parameters to go from gamma to standard distribution function for precip (Magyari-Saska, 2009) 

monnum=4 # month (here for testing)



fnamebase='G:/steppe_climate_data/seasonal_can/' # path to historical hindcasts
foutbase='G:/steppe_climate_data/SPI_for_web/out/' 

		
#cut gridded data by rectangular and from start to end year
for im in range(0,1):	
	imm=im+monnum
	ny=0
	precip_s_all=[]
#Read hindcasts year by year and cut by rectangular
	for iy in range(yst, yfin+1):
		ny=ny+1
		filename = (fnamebase+'cansips_hindcast_raw_latlon1.0x1.0_PRATE_SFC_0_%s-0%s_allmembers.nc' %(iy,imm))
		ncfile = Dataset(filename, 'r') 
		latitude_s = ncfile.variables['lat'][:]
		longitude_s = ncfile.variables['lon'][:]
		time = ncfile.variables['time'][:] 
		precip_s0= ncfile.variables['prate'][:]	
		ncfile.close()
		a1=precip_s0.shape[1]
		a2=precip_s0.shape[2]
		long=a1*a2
		precip_s=np.reshape(precip_s0, (1, precip_s0.shape[0], long))
		precip_s_all=np.reshape(np.append(precip_s_all, precip_s), (ny*1, precip_s0.shape[0], long))
	filenameo= (foutbase+'glob_a_b_%s.nc' %imm)
	ncfileo = Dataset(filenameo, 'w') 
	lat_dim = ncfileo.createDimension('lat', len(latitude_s)) # latitude axis
	lon_dim = ncfileo.createDimension('lon', len(longitude_s)) # longitude axis
	time_dim = ncfileo.createDimension('time', None) 
	lat = ncfileo.createVariable('lat', np.float32, ('lat',))
	lat.units = 'degrees_north'
	lat.long_name = 'latitude'
	lon = ncfileo.createVariable('lon', np.float32, ('lon',))
	lon.units = 'degrees_east'
	lon.long_name = 'longitude'
	timeo = ncfileo.createVariable('time', np.float64, ('time',))
	timeo.long_name = 'time'
	outA= ncfileo.createVariable('A_coeff',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
	outA.standard_name = 'A_coeff' 	
	outB= ncfileo.createVariable('B_coeff',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
	outB.standard_name = 'B_coeff' 	
	outQ= ncfileo.createVariable('q',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
	outQ.standard_name = 'q' 		
	lat[:]=latitude_s[:]
	lon[:]=longitude_s[:]
	timeo=time[0]
# SPI simulation for historical forecasts (hindcasts) from Canadian Prediction Center
	fit_alpha_s=np.zeros(long)
	fit_loc_s=np.zeros(long)
	fit_beta_s=np.zeros(long)
	qs=np.zeros(long)	
	for ig in range(0, long):
		precip_s_amean=np.mean(precip_s_all[:,:,ig], axis=1)
		prs=precip_s_amean[:]
		qs[ig]=len(prs[prs<=0.])/len(prs)
		precip_snonz=prs[prs>0.]
		fit_alpha_s[ig], fit_loc_s[ig], fit_beta_s[ig] = stats.gamma.fit(precip_snonz, floc=0 ) 
		
	A=np.reshape(fit_alpha_s, (a1, a2))	
	B=np.reshape(fit_beta_s, (a1, a2))
	qq=np.reshape(qs, (a1, a2))
	
	outA[0,:,:] = A[:, :]
	outB[0,:,:] = B[:, :]
	outQ[0,:,:] = qq[:, :]
	ncfileo.close()		