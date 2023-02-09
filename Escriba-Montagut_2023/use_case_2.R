# Use case 2: Exposome info from geo data

library(dsBase)
library(dsBaseClient)
library(resourcer)
library(DSLite)
library(dsExposomeClient)
library(dsExposome)

nc_data_pm25 <- newResource(
  name = "nc_data_pm25",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/V4NA03_PM25_NA_201801_201812-RH35.nc",
  format = "NetCDF"
)
nc_data_ss <- newResource(
  name = "nc_data_ss",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/GWRwSPEC_SSp_NA_201401_201412-wrtSPECtotal.nc",
  format = "NetCDF"
)
nc_data_so4 <- newResource(
  name = "nc_data_so4",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/GWRwSPEC_SO4p_NA_201501_201512-wrtSPECtotal.nc",
  format = "NetCDF"
)
nc_data_no3 <- newResource(
  name = "nc_data_no3",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/GWRwSPEC_NITp_NA_201201_201212-wrtSPECtotal.nc",
  format = "NetCDF"
)
nc_data_nh4 <- newResource(
  name = "nc_data_nh4",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/GWRwSPEC_NH4p_NA_201601_201612-wrtSPECtotal.nc",
  format = "NetCDF"
)
nc_data_om <- newResource(
  name = "nc_data_om",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/GWRwSPEC_OMp_NA_201301_201312-wrtSPECtotal.nc",
  format = "NetCDF"
)
nc_data_bc <- newResource(
  name = "nc_data_bc",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/GWRwSPEC_BCp_NA_201401_201412-wrtSPECtotal.nc",
  format = "NetCDF"
)
nc_data_soil <- newResource(
  name = "nc_data_soil",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/GWRwSPEC_SOILp_NA_201301_201312-wrtSPECtotal.nc",
  format = "NetCDF"
)
clinical <- newResource(
  name = "clinical",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/data4Xavi.Rdata",
  format = "data.frame"
)

dslite.server <- newDSLiteServer(resources=list(
  nc_data_pm25 = nc_data_pm25,
  nc_data_ss = nc_data_ss,
  nc_data_so4 = nc_data_so4,
  nc_data_no3 = nc_data_no3,
  nc_data_nh4 = nc_data_nh4,
  nc_data_om = nc_data_om,
  nc_data_bc = nc_data_bc,
  nc_data_soil = nc_data_soil,
  clinical = clinical
),
config = DSLite::defaultDSConfiguration(include=c("dsBase", 
                                                  "resourcer", 
                                                  "dsExposome")))
builder <- DSI::newDSLoginBuilder()
builder$append(server = "server1", url = "dslite.server", driver = "DSLiteDriver")
logindata.dslite.cnsim <- builder$build()
conns <- datashield.login(logindata.dslite.cnsim, assign=F)

datashield.assign.resource(conns = conns, symbol = "nc_data_pm25", resource = "nc_data_pm25")
datashield.assign.resource(conns = conns, symbol = "nc_data_ss", resource = "nc_data_ss")
datashield.assign.resource(conns = conns, symbol = "nc_data_so4", resource = "nc_data_so4")
datashield.assign.resource(conns = conns, symbol = "nc_data_no3", resource = "nc_data_no3")
datashield.assign.resource(conns = conns, symbol = "nc_data_nh4", resource = "nc_data_nh4")
datashield.assign.resource(conns = conns, symbol = "nc_data_om", resource = "nc_data_om")
datashield.assign.resource(conns = conns, symbol = "nc_data_bc", resource = "nc_data_bc")
datashield.assign.resource(conns = conns, symbol = "nc_data_soil", resource = "nc_data_soil")
datashield.assign.resource(conns = conns, symbol = "clinical", resource = "clinical")

datashield.assign.expr(conns = conns, symbol = "nc_data_pm25", expr = quote(as.resource.object(nc_data_pm25)))
datashield.assign.expr(conns = conns, symbol = "nc_data_ss", expr = quote(as.resource.object(nc_data_ss)))
datashield.assign.expr(conns = conns, symbol = "nc_data_so4", expr = quote(as.resource.object(nc_data_so4)))
datashield.assign.expr(conns = conns, symbol = "nc_data_no3", expr = quote(as.resource.object(nc_data_no3)))
datashield.assign.expr(conns = conns, symbol = "nc_data_nh4", expr = quote(as.resource.object(nc_data_nh4)))
datashield.assign.expr(conns = conns, symbol = "nc_data_om", expr = quote(as.resource.object(nc_data_om)))
datashield.assign.expr(conns = conns, symbol = "nc_data_bc", expr = quote(as.resource.object(nc_data_bc)))
datashield.assign.expr(conns = conns, symbol = "nc_data_soil", expr = quote(as.resource.object(nc_data_soil)))
datashield.assign.expr(conns = conns, symbol = "clinical", expr = quote(as.resource.object(clinical)))

nc_data_pm25_var <- ds.netcdf_vars("nc_data_pm25")[[1]]
nc_data_ss_var <- ds.netcdf_vars("nc_data_ss")[[1]]
nc_data_so4_var <- ds.netcdf_vars("nc_data_so4")[[1]]
nc_data_no3_var <- ds.netcdf_vars("nc_data_no3")[[1]]
nc_data_nh4_var <- ds.netcdf_vars("nc_data_nh4")[[1]]
nc_data_om_var <- ds.netcdf_vars("nc_data_om")[[1]]
nc_data_bc_var <- ds.netcdf_vars("nc_data_bc")[[1]]
nc_data_soil_var <- ds.netcdf_vars("nc_data_soil")[[1]]

ds.ncvar_get("nc_data_pm25", "LON", "lon_pm25")
ds.ncvar_get("nc_data_ss", "LON", "lon_ss")
ds.ncvar_get("nc_data_so4", "LON", "lon_so4")
ds.ncvar_get("nc_data_no3", "LON", "lon_no3")
ds.ncvar_get("nc_data_nh4", "LON", "lon_nh4")
ds.ncvar_get("nc_data_om", "LON", "lon_om")
ds.ncvar_get("nc_data_bc", "LON", "lon_bc")
ds.ncvar_get("nc_data_soil", "LON", "lon_soil")

ds.ncvar_get("nc_data_pm25", "LAT", "lat_pm25")
ds.ncvar_get("nc_data_ss", "LAT", "lat_ss")
ds.ncvar_get("nc_data_so4", "LAT", "lat_so4")
ds.ncvar_get("nc_data_no3", "LAT", "lat_no3")
ds.ncvar_get("nc_data_nh4", "LAT", "lat_nh4")
ds.ncvar_get("nc_data_om", "LAT", "lat_om")
ds.ncvar_get("nc_data_bc", "LAT", "lat_bc")
ds.ncvar_get("nc_data_soil", "LAT", "lat_soil")

ds.ncvar_get("nc_data_pm25", nc_data_pm25_var, nc_data_pm25_var)
ds.ncvar_get("nc_data_ss", nc_data_ss_var, nc_data_ss_var)
ds.ncvar_get("nc_data_so4", nc_data_so4_var, nc_data_so4_var)
ds.ncvar_get("nc_data_no3", nc_data_no3_var, nc_data_no3_var)
ds.ncvar_get("nc_data_nh4", nc_data_nh4_var, nc_data_nh4_var)
ds.ncvar_get("nc_data_om", nc_data_om_var, nc_data_om_var)
ds.ncvar_get("nc_data_bc", nc_data_bc_var, nc_data_bc_var)
ds.ncvar_get("nc_data_soil", nc_data_soil_var, nc_data_soil_var)

ds.ncatt_get("nc_data_pm25", nc_data_pm25_var, "_FillValue", "fillvalue_pm25")
ds.ncatt_get("nc_data_ss", nc_data_ss_var, "_FillValue", "fillvalue_ss")
ds.ncatt_get("nc_data_so4", nc_data_so4_var, "_FillValue", "fillvalue_so4")
ds.ncatt_get("nc_data_no3", nc_data_no3_var, "_FillValue", "fillvalue_no3")
ds.ncatt_get("nc_data_nh4", nc_data_nh4_var, "_FillValue", "fillvalue_nh4")
ds.ncatt_get("nc_data_om", nc_data_om_var, "_FillValue", "fillvalue_om")
ds.ncatt_get("nc_data_bc", nc_data_bc_var, "_FillValue", "fillvalue_bc")
ds.ncatt_get("nc_data_soil", nc_data_soil_var, "_FillValue", "fillvalue_soil")

ds.NetCDF_fillvalue_matrix(nc_data_pm25_var, "fillvalue_pm25", nc_data_pm25_var)
ds.NetCDF_fillvalue_matrix(nc_data_ss_var, "fillvalue_ss", nc_data_ss_var)
ds.NetCDF_fillvalue_matrix(nc_data_so4_var, "fillvalue_so4", nc_data_so4_var)
ds.NetCDF_fillvalue_matrix(nc_data_no3_var, "fillvalue_no3", nc_data_no3_var)
ds.NetCDF_fillvalue_matrix(nc_data_nh4_var, "fillvalue_nh4", nc_data_nh4_var)
ds.NetCDF_fillvalue_matrix(nc_data_om_var, "fillvalue_om", nc_data_om_var)
ds.NetCDF_fillvalue_matrix(nc_data_bc_var, "fillvalue_bc", nc_data_bc_var)
ds.NetCDF_fillvalue_matrix(nc_data_soil_var, "fillvalue_soil", nc_data_soil_var)

ds.get_exposure_from_geo("lat_pm25",
                         "lon_pm25", 
                         nc_data_pm25_var, 
                         nc_data_pm25_var, 
                         "clinical", "lat", "lon", "id",
                         new.obj = "exposures_new")
ds.get_exposure_from_geo("lat_ss", 
                         "lon_ss", 
                         nc_data_ss_var, 
                         nc_data_ss_var, 
                         "clinical", "lat", "lon", "id", 
                         "exposures_new", "id", "exposures_new")
ds.get_exposure_from_geo("lat_so4", 
                         "lon_so4", 
                         nc_data_so4_var, 
                         nc_data_so4_var, 
                         "clinical", "lat", "lon", "id", 
                         "exposures_new", "id", "exposures_new")
ds.get_exposure_from_geo("lat_no3", 
                         "lon_no3", 
                         nc_data_no3_var, 
                         nc_data_no3_var, 
                         "clinical", "lat", "lon", "id", 
                         "exposures_new", "id", "exposures_new")
ds.get_exposure_from_geo("lat_nh4", 
                         "lon_nh4", 
                         nc_data_nh4_var, 
                         nc_data_nh4_var, 
                         "clinical", "lat", "lon", "id", 
                         "exposures_new", "id", "exposures_new")
ds.get_exposure_from_geo("lat_om",
                         "lon_om", 
                         nc_data_om_var, 
                         nc_data_om_var, 
                         "clinical", "lat", "lon", "id", 
                         "exposures_new", "id", "exposures_new")
ds.get_exposure_from_geo("lat_bc", 
                         "lon_bc", 
                         nc_data_bc_var, 
                         nc_data_bc_var, 
                         "clinical", "lat", "lon", "id", 
                         "exposures_new", "id", "exposures_new")
ds.get_exposure_from_geo("lat_soil",
                         "lon_soil", 
                         nc_data_soil_var, 
                         nc_data_soil_var, 
                         "clinical", "lat", "lon", "id", 
                         "exposures_new", "id", "exposures_new")

ds.make(toAssign = "c('Expos', 'Expos', 'Expos', 'Expos', 'Expos', 'Expos', 'Expos', 'Expos')", newobj = "fams")
ds.make(toAssign = paste0("c('",
                          paste(c(nc_data_pm25_var,
                          nc_data_ss_var,
                          nc_data_so4_var,
                          nc_data_no3_var,
                          nc_data_nh4_var,
                          nc_data_om_var,
                          nc_data_bc_var,
                          nc_data_soil_var), collapse = "','"), "')"), 
        newobj = "expos")

ds.dataFrame(x = c("fams", "expos"), newobj = "desc")

dsExposomeClient::ds.loadExposome(
  exposures = "exposures_new",
  phenotypes = "clinical",
  exposures.idcol = "id",
  phenotypes.idcol = "id",
  object_name = "expo_set2",
  description = "desc",
  description.expCol = "expos",
  description.famCol = "fams",
  exposures.asFactor = 2
)

exwas_results <- ds.exwas("ihd ~ age + smoke", Set = "expo_set2", family = "binomial", exposures_family = "Expos", type = "meta")
ds.plotExwas(exwas_results, thld_pvalue = 0.05)
