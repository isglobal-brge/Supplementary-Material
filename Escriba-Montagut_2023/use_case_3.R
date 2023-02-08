library(dsBase)
library(dsBaseClient)
library(resourcer)
library(DSLite)
library(dsExposomeClient)
library(dsExposome)

exposome_data <- newResource(
  name = "exposome_data",
  url = "file:C://Users/Xavier/OneDrive/ISGlobal/PhD/Paper-dsExposome/use_cases/data/exposome.RData",
  format = "data.frame"
)


dslite.server <- newDSLiteServer(resources=list(
  exposome_data = exposome_data
),
config = DSLite::defaultDSConfiguration(include=c("dsBase", 
                                                  "resourcer", 
                                                  "dsExposome")))
builder <- DSI::newDSLoginBuilder()
builder$append(server = "server1", url = "dslite.server", driver = "DSLiteDriver")
logindata.dslite.cnsim <- builder$build()
conns <- datashield.login(logindata.dslite.cnsim, assign=F)

datashield.assign.resource(conns, "res", "exposome_data")
datashield.assign.expr(conns, "test", quote(as.resource.data.frame(res, 2)))
