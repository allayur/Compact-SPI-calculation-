import xarray as xr 
import rioxarray as rio 
fnamebase='G:/steppe_climate_data/SPI_for_web/out/' # path to output files
yst=1993 #start of SPI calibration (year)
monnum=1 
for i in range(0,1): # 12- all months
	im=i+monnum
	nc_file = xr.open_dataset(fnamebase+'SPI_rean_%s.nc' %im )
	bTl = nc_file['SPI']
	for i in range(0, bTl.shape[0]):
		iy=yst+i
		bT=bTl[i, :,:]
		bT = bT.rio.set_spatial_dims('lon', 'lat')
		bT.rio.crs
		bT.rio.set_crs("epsg:4326")
		if (im<10):
			fname=(fnamebase+"SPI_%s_0%s.tiff" %(iy,im))
		else:
			fname=(fnamebase+"SPI_%s_%s.tiff" %(iy,im))
		bT.rio.to_raster(fname)