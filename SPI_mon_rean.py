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
start = date(1800,1,1)
#Parameters to go from gamma to standard distribution function for precip (Magyari-Saska, 2009) 
c0=2.515517
c1=0.802853
c2=0.010328
d1=1.432788 
d2=0.189269 
d3=0.001308
monnum=1 # month (here for testing)
thresh=0.9
noval=-9999

filename2 = 'G:/steppe_climate_data/precip/precip.mon.mean.1x1.nc' # path to historical gridded data
fnamebase='G:/steppe_climate_data/SPI_for_web/out/' # path to output files

# Read gridded observation data
ncfile2 = Dataset(filename2, 'r') 
latitudes = ncfile2.variables['lat'][:]
longitudes = ncfile2.variables['lon'][:]
time = ncfile2.variables['time'][:]   
precip0= ncfile2.variables['precip'][:]
ncfile2.close()

#define target corners
dlat=np.abs(latitudes[:]-llclat)
dlon=np.abs(longitudes[:]-llclon)
all=np.where(dlat==np.min(dlat))
bll=np.where(dlon==np.min(dlon))
dlat2=np.abs(latitudes[:]-urclat)
dlon2=np.abs(longitudes[:]-urclon)
aur=np.where(dlat2==np.min(dlat2))
bur=np.where(dlon2==np.min(dlon2))

#Define the time slice from datetime
mm = np.ones(len(time))
yyyy = np.ones(len(time))
dd = np.ones(len(time))

for j in range(0,len(time)):
		delta = timedelta(days=np.int(time[j]/24))
		o = start + delta
		mm[j] = o.month
		yyyy[j] = o.year
		dd[j] = o.day
		
#cut gridded data by rectangular and from start to end year
for im in range(0,1):	
	imm=im+monnum
	ta=np.where((yyyy>=yst) )
	tb=np.where(yyyy==yfin)	
	mms=mm[(yyyy>=yst) & (yyyy<yfin+1)]
	precip=precip0[ta[0][0]:tb[0][0]+12, aur[0][0]:all[0][0]+1, bll[0][0]:bur[0][0]+1]
	time1=time[ta[0][0]:tb[0][0]+12]
	ytest=yyyy[ta[0][0]:tb[0][0]+12]
	a1=precip.shape[1]
	a2=precip.shape[2]
	lat1=latitudes[aur[0][0]:all[0][0]+1]
	lon1=longitudes[bll[0][0]:bur[0][0]+1]
	long=a1*a2
	precip=np.reshape(precip, (precip.shape[0], long))
	precip=precip[(mms==im+monnum)]
	time_out=time1[(mms==im+monnum)]
	
	precip_s_all=[]
	ny=0
	#Prepare output netcdf files
	filenameo= (fnamebase+'SPI_rean_%s.nc' %imm)
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
	timeo[:]=time_out[:]	

# SPI calculation for historical gridded data		
	fit_alpha=np.zeros(long)
	fit_loc=np.zeros(long)
	fit_beta=np.zeros(long)
	SPI=np.zeros((precip.shape[0], long))
	SPI[:,:]=noval
	for ig in range(0, long):		
		pr=precip[:,ig]
		q=len(pr[pr<=0.])/len(pr)
		prnonz=pr[pr>0.]
		fit_alpha[ig], fit_loc[ig], fit_beta[ig] = stats.gamma.fit(prnonz	, floc=0 ) 
		if (q<thresh):
			for it in range(0,precip.shape[0]):
				G=gamma.cdf(precip[it,ig], fit_alpha[ig], scale=fit_beta[ig])
				H=q+(1-q)*G			
				if (H<0.5):
					t=np.sqrt(np.log(1./(H*H)))
					SPI[it,ig]=-(t-((c0+c1*t+c2*t*t)/(1+d1*t+d2*t*t+d3*t*t*t)))	
				else:	
					t=np.sqrt(np.log(1./((1.-H)*(1-H))))
					SPI[it,ig]=t-((c0+c1*t+c2*t*t)/(1+d1*t+d2*t*t+d3*t*t*t))
		else:
			SPI[it,ig]=noval
	SPI_out_rean=np.reshape(SPI, (SPI.shape[0], a1, a2))		
	outSPI[:,:,:] = SPI_out_rean[:, :, :]
	ncfileo.close()
	