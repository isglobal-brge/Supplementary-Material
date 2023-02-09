require(DSOpal)
require(DSI)
require(dsBaseClient)
require(dsExposomeClient)
require(dsMTLClient)

# Loading additional required packages
require(dplyr)
require(ggplot2)
require(ggrepel)
require(EnhancedVolcano)
require(tidyverse)
require(ggrepel)
require(reshape2)
require(RColorBrewer)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "BIB", url = "https://opal.isglobal.org/repo",
               user =  "invited", password = "Invited88_@",
               profile = "rock-inma")
builder$append(server = "EDEN", url = "https://opal.isglobal.org/repo",
               user =  "invited", password = "Invited88_@",
               profile = "rock-inma")
builder$append(server = "KANC", url = "https://opal.isglobal.org/repo",
               user =  "invited", password = "Invited88_@",
               profile = "rock-inma")
builder$append(server = "MoBA", url = "https://opal.isglobal.org/repo",
               user =  "invited", password = "Invited88_@",
               profile = "rock-inma")
builder$append(server = "Rhea", url = "https://opal.isglobal.org/repo",
               user =  "invited", password = "Invited88_@",
               profile = "rock-inma")
builder$append(server = "INMASAB", url = "https://opal.isglobal.org/repo",
               user =  "invited", password = "Invited88_@",
               profile = "rock-inma")
logindata <- builder$build()
conns <- DSI::datashield.login(logins = logindata)

# We assign post-natal data from all cohorts to an object called resource_pos
DSI::datashield.assign.resource(conns[1], "resource_pos", "HELIX.postnatal_BIB")
DSI::datashield.assign.resource(conns[2], "resource_pos", "HELIX.postnatal_EDE")
DSI::datashield.assign.resource(conns[3], "resource_pos", "HELIX.postnatal_KAN")
DSI::datashield.assign.resource(conns[4], "resource_pos", "HELIX.postnatal_MOB")
DSI::datashield.assign.resource(conns[5], "resource_pos", "HELIX.postnatal_RHE")
DSI::datashield.assign.resource(conns[6], "resource_pos", "HELIX.postnatal_SAB")

# We resolve the resource for post-natal data
DSI::datashield.assign.expr(conns = conns, symbol = "exposome_set_pos",
                            expr = as.symbol("as.resource.object(resource_pos)"))

DSI::datashield.assign.resource(conns[1], "pheno", "HELIX.subclinical_cardio_BIB")
DSI::datashield.assign.resource(conns[2], "pheno", "HELIX.subclinical_cardio_EDE")
DSI::datashield.assign.resource(conns[3], "pheno", "HELIX.subclinical_cardio_KAN")
DSI::datashield.assign.resource(conns[4], "pheno", "HELIX.subclinical_cardio_MOB")
DSI::datashield.assign.resource(conns[5], "pheno", "HELIX.subclinical_cardio_RHE")
DSI::datashield.assign.resource(conns[6], "pheno", "HELIX.subclinical_cardio_SAB")

DSI::datashield.assign.expr(conns = conns, symbol = "pheno_dt",
                            expr = as.symbol("as.resource.data.frame(pheno, strict = TRUE)"))

# Adding phenotype data to post-natal datasets.
ds.addPhenoData2ExposomeSet("exposome_set_pos", "pheno_dt", 
                            identifier_ExposomeSet = "HelixID", 
                            identifier_new_phenotypes = "HelixID")

# Adjust by cohort and other confounders. 
# Estimated execution time: approx. 6.830064 mins
res.pooled_adjusted_pos <- ds.exwas("hs_bp_sys ~ hs_child_age_days_None + hs_c_height_None + e3_sex_None + h_age_None + h_mbmi_None", type = "pooled", 
                                    adjust.by.study = TRUE,
                                    exposures_family = "OCs", 
                                    Set = "exposome_set_pos", family = "gaussian",
                                    tef = FALSE)
res.pooled_adjusted_pos2 <- ds.exwas("hs_bp_sys ~ hs_child_age_days_None + hs_c_height_None + e3_sex_None + h_age_None + h_mbmi_None", type = "pooled", 
                                    adjust.by.study = TRUE,
                                    exposures_family = "Indoor air", 
                                    Set = "exposome_set_pos", family = "gaussian",
                                    tef = FALSE)
res.pooled_adjusted_pos3 <- ds.exwas("hs_bp_sys ~ hs_child_age_days_None + hs_c_height_None + e3_sex_None + h_age_None + h_mbmi_None", type = "pooled", 
                                    adjust.by.study = TRUE,
                                    exposures_family = "Phthalates", 
                                    Set = "exposome_set_pos", family = "gaussian",
                                    tef = FALSE)
res.pooled_adjusted_pos4 <- ds.exwas("hs_bp_sys ~ hs_child_age_days_None + hs_c_height_None + e3_sex_None + h_age_None + h_mbmi_None", type = "pooled", 
                                    adjust.by.study = TRUE,
                                    exposures_family = "PBDEs", 
                                    Set = "exposome_set_pos", family = "gaussian",
                                    tef = FALSE)
