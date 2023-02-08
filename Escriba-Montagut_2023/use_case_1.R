# Use case 1: Multi-cohort ExWAS

library(dsBaseClient)
library(DSLite)
library(resourcer)
library(dsExposomeClient)

expo_set1 <- newResource(
  name = "expo_set1",
  url = "file:C:/Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/exposome1.Rdata",
  format = "ExposomeSet"
)
expo_set2 <- newResource(
  name = "expo_set2",
  url = "file:C:/Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/exposome2.Rdata",
  format = "ExposomeSet"
)

dslite.server <- newDSLiteServer(resources=list(expo_set1 = expo_set1, expo_set2 = expo_set2),
                                 config = DSLite::defaultDSConfiguration(
                                   include=c("dsBase", 
                                             "resourcer", 
                                             "dsExposome")
                                 ))
builder <- DSI::newDSLoginBuilder()
builder$append(server = "server1", url = "dslite.server", resource = "expo_set1", driver = "DSLiteDriver")
builder$append(server = "server2", url = "dslite.server", resource = "expo_set2", driver = "DSLiteDriver")
logindata.dslite.cnsim <- builder$build()
conns <- datashield.login(logindata.dslite.cnsim, assign=T, symbol = "expos")
datashield.assign.expr(conns = conns, symbol = "exposome_set", expr = quote(as.resource.object(expos)))

missings <- ds.tableMissings("exposome_set")
ds.plotMissings(missings)

ds.familyNames("exposome_set")
ds.exposome_variables("exposome_set", target = "phenotypes")

res <- ds.exwas(
  model = "flu ~ 1", 
  Set = "exposome_set", 
  family = "binomial", 
  type = "pooled", 
  exposures_family = c("Air Pollutants", "Metals"),
  tef = FALSE
)

a <- ds.plotExwas(res, thld_pvalue = 0.05)
a
