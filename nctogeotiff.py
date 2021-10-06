import xarray as xr 
import rioxarray as rio 
yst=1993 #start of SPI calibration (year)
im=4
nc_file = xr.open_dataset('SPI_rean_%s.nc' %im )
bTl = nc_file['SPI']
print (bTl.shape)
for i in range(0, bTl.shape[0]):
	iy=yst+i
	bT=bTl[i, :,:]
	bT = bT.rio.set_spatial_dims('lon', 'lat')
	bT.rio.crs
	bT.rio.set_crs("epsg:4326")
	bT.rio.to_raster("SPI_%s_0%s.tiff" %(iy,im))