import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime
from numpy import arange, dtype 
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.stats as stats 
from datetime import date, timedelta,  datetime
from scipy.stats import gamma
import warnings
import xarray as xr 
import rioxarray as rio 
import cfgrib
import urllib.request

warnings.filterwarnings("ignore")
#define the corners of the target rectangular
llclat=24.  
llclon=-124. +360.
urclat=50. 
urclon=-67. +360.

#Parameters to go from gamma to standard distribution function for precip (Magyari-Saska, 2009) 
c0=2.515517
c1=0.802853
c2=0.010328
d1=1.432788 
d2=0.189269 
d3=0.001308

thresh=0.9
noval=-9999

falfabeta='G:/steppe_climate_data/SPI_for_web/out/' 
fforename='G:/steppe_climate_data/SPI_for_web/forecast/' 

#read system time

date=datetime.today()
mon=date.month
yr=date.year

#Connect to url and load the latest forecast
if (mon>=10):
	urllib.request.urlretrieve("https://dd.weather.gc.ca/ensemble/cansips/grib2/forecast/raw/%s/%s/cansips_forecast_raw_latlon1.0x1.0_PRATE_SFC_0_%s-%s_allmembers.grib2" %(yr, mon, yr, mon), "cansips_forecast_raw_latlon1.0x1.0_PRATE_SFC_0_%s-%s_allmembers.grib2" %(yr, mon))
else:	
	urllib.request.urlretrieve("https://dd.weather.gc.ca/ensemble/cansips/grib2/forecast/raw/%s/0%s/cansips_forecast_raw_latlon1.0x1.0_PRATE_SFC_0_%s-0%s_allmembers.grib2" %(yr, mon, yr, mon), "cansips_forecast_raw_latlon1.0x1.0_PRATE_SFC_0_%s-0%s_allmembers.grib2" %(yr, mon))

filename2= ("cansips_forecast_raw_latlon1.0x1.0_PRATE_SFC_0_%s-0%s_allmembers.grib2" %(yr, mon)) 
ds = xr.load_dataset(filename2, engine='cfgrib')
prec = ds['prate']

#Load global predifined gamme fit

filename= (falfabeta+'glob_a_b_%s.nc' %mon)
ncfile = Dataset(filename, 'r') 
latitudes = ncfile.variables['lat'][:]
longitudes = ncfile.variables['lon'][:]
alfa= ncfile.variables['A_coeff'][:]
beta= ncfile.variables['B_coeff'][:]
qq= ncfile.variables['q'][:]

ncfile.close()
a1=prec.shape[1]
a2=prec.shape[2]

long=a1*a2
fit_alpha=np.reshape(alfa, long)
fit_beta=np.reshape(beta, long)
q=np.reshape(qq, long)
precip_s=np.reshape(prec, (prec.shape[0], long))

#Calculate SPI globally  using predefined gamma fit
SPI_s=np.zeros((long))
SPI_s[:]=noval
for ig in range(0, long):	
	precip_mean=np.mean(precip_s[:,ig], axis=0)	
	qs=q[ig]
	G=gamma.cdf(precip_mean[ig], fit_alpha[ig], scale=fit_beta[ig])
	H=qs+(1-qs)*G
	if (qs<thresh):
		if (H<0.5):
			t=np.sqrt(np.log(1./(H*H)))
			SPI_s[ig]=-(t-((c0+c1*t+c2*t*t)/(1+d1*t+d2*t*t+d3*t*t*t)))
		else:						
			t=np.sqrt(np.log(1./((1.-H)*(1-H))))
			SPI_s[ig]=t-((c0+c1*t+c2*t*t)/(1+d1*t+d2*t*t+d3*t*t*t))
	else:
		SPI_s[ig]=noval
SPI_out_fore_glob=np.reshape(SPI_s, (a1, a2))			

# Cut by rectangular
dlat=np.abs(latitudes[:]-llclat)
dlon=np.abs(longitudes[:]-llclon)
all=np.where(dlat==np.min(dlat))
bll=np.where(dlon==np.min(dlon))
dlat2=np.abs(latitudes[:]-urclat)
dlon2=np.abs(longitudes[:]-urclon)
aur=np.where(dlat2==np.min(dlat2))
bur=np.where(dlon2==np.min(dlon2))
SPI_out_fore=SPI_out_fore_glob[all[0][0]:aur[0][0]+1, bll[0][0]:bur[0][0]+1]
lat1=latitudes[all[0][0]:aur[0][0]+1]
lon1=longitudes[bll[0][0]:bur[0][0]+1]

#Write to netcdf
if (im<10):
	filenameo=(fforename+"SPI_forecast_%s_0%s.nc" %(yr,mon))
else:
	filenameo=(fforename+"SPI_forecast_%s_%s.nc" %(yr,mon))
ncfileo = Dataset(filenameo, 'w') 
lat_dim = ncfileo.createDimension('lat', len(lat1)) # latitude axis
lon_dim = ncfileo.createDimension('lon', len(lon1)) # longitude axis
time_dim = ncfileo.createDimension('time', None) 
lat = ncfileo.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfileo.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
timeo = ncfileo.createVariable('time', np.float64, ('time',))
timeo.units = 'hours since 1800-01-01'
timeo.long_name = 'time'
outSPI= ncfileo.createVariable('SPI',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
outSPI.standard_name = 'SPI' 
outSPI.valid_range = (-4.0, 4.0)
lat[:]=lat1[:]
lon[:]=lon1[:]
timeo=date
outSPI[0,:,:] = SPI_out_fore[:, :]

#Write to geotiff
SPI_out_fore = SPI_out_fore.rio.set_spatial_dims('lon', 'lat')
SPI_out_fore.rio.crs
SPI_out_fore.rio.set_crs("epsg:4326")
if (mon<10):
	fname=(fforename+"SPI_%s_0%s.tiff" %(yr,mon))
else:
	fname=(fforename+"SPI_%s_%s.tiff" %(yr,mon))
SPI_out_fore.rio.to_raster(fname)
