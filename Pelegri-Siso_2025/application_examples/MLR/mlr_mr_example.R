library(BigDataStatMeth)

setwd("/Users/mailos/DOCTORAT_Local/BigDataStatMeth/Colesterol_MLR")

# Define url
modelfile <- "https://raw.githubusercontent.com/isglobal-brge/Supplementary-Material/master/Pelegri-Siso_2021/application_examples/MLR/data/Cholesterol_model.zip"
Yfile <- "https://raw.githubusercontent.com/isglobal-brge/Supplementary-Material/master/Pelegri-Siso_2021/application_examples/MLR/data/Cholesterol_Y.zip"


# Import data from url
bdImportData_hdf5( inFile = modelfile, 
                   destFile = "mlr_mr_example.hdf5", 
                   destGroup = "data", 
                   destDataset = "X",
                   header = TRUE,
                   rownames = FALSE, 
                   overwrite = TRUE)
bdImportData_hdf5( inFile = Yfile,
                   destFile = "mlr_mr_example.hdf5",
                   destGroup = "data",
                   destDataset = "Y",
                   header = TRUE, 
                   rownames = FALSE, 
                   overwrite = TRUE)


# Set number of partitions
m <- 100


# 
# # Step 1 :
# Split datasets X abd Y by rows and store data to data file
bdSplit_matrix_hdf5( filename = "mlr_mr_example.hdf5",
                     group = "data",
                     dataset = "X",
                     outgroup = "Step1/Xrows",
                     nblocks = m, bycols = FALSE,
                     force = TRUE)
bdSplit_matrix_hdf5( filename = "mlr_mr_example.hdf5",
                     group = "data", 
                     dataset = "Y", 
                     outgroup = "Step1/Yrows", 
                     nblocks = m, 
                     bycols = FALSE, 
                     force = TRUE)

# Step 2 :
# Get splitted dataset names and apply QR to all blocks
x.blocks <- bdgetDatasetsList_hdf5("mlr_mr_example.hdf5", "Step1/Xrows")
bdapply_Function_hdf5( filename = "mlr_mr_example.hdf5", group = "Step1/Xrows", 
                       datasets = x.blocks, 
                       outgroup = "Step2/Xrows", 
                       func = "QR", 
                       force = TRUE )


#  Step 3 :
# Merge R in Rt
x.blocks.qr <- bdgetDatasetsList_hdf5("mlr_mr_example.hdf5", "Step2/Xrows")
bdBind_hdf5(filename = "mlr_mr_example.hdf5", 
            group = "Step2/Xrows", 
            datasets = x.blocks.qr[which(x.blocks.qr %like% ".R")],
            outgroup = "Step3/merged", outdataset = "Rt", 
            func = "bindRows", force = TRUE )

bdapply_Function_hdf5( filename = "mlr_mr_example.hdf5", 
                       group = "Step3/merged", 
                       datasets = "Rt", 
                       outgroup = "Step3/Final_QR", 
                       func = "QR", 
                       force = TRUE )

# Step 4 :
bdSplit_matrix_hdf5("mlr_mr_example.hdf5", "Step3/Final_QR", "Rt.Q", 
                    outgroup = "Step4/splitted", 
                    nblocks = m, 
                    bycols = FALSE, force = TRUE )

# Step 5 :
# Get splitted matrices names
Rt.Q.divide <- bdgetDatasetsList_hdf5("mlr_mr_example.hdf5", "Step4/splitted")
# multiply previous splitted matrices with Q descomposed matrices from model (X)
bdapply_Function_hdf5(  filename = "mlr_mr_example.hdf5", group = "Step2/Xrows", 
                        datasets = x.blocks.qr[which(x.blocks.qr %like% ".Q")], 
                        outgroup = "Step5", func = "blockmult",
                        b_group = "Step4/splitted", b_datasets = Rt.Q.divide,
                        force = TRUE )

# Step 6 :
QResults <- bdgetDatasetsList_hdf5("mlr_mr_example.hdf5", "Step5")
y.block <- bdgetDatasetsList_hdf5("mlr_mr_example.hdf5", "Step1/Yrows")

bdapply_Function_hdf5(  filename = "mlr_mr_example.hdf5", group = "Step5", 
                        datasets = QResults, 
                        outgroup = "Step6/V", func = "CrossProd",
                        b_group = "Step1/Yrows", b_datasets = y.block,
                        force = TRUE )
# Step 7 :
bdReduce_matrix_hdf5( filename ="mlr_mr_example.hdf5", 
                      group = "Step6/V", 
                      reducefunction = "+", 
                      outgroup = "Step7/Final_QR", 
                      outdataset = "V.sum", 
                      force = TRUE)

# Load results from hdf5 data file to get final Betas using the package rhdf5
h5f = H5Fopen("mlr_mr_example.hdf5")
    Rfinal <- h5f$Step3$Final_QR$Rt.R
    V <- h5f$Step7$Final_QR$V.sum
    vars <- h5f$data$.X_dimnames$`2`
h5closeAll()

beta <- bdSolve(A = Rfinal, B = diag(1, nrow(Rfinal), ncol(Rfinal))) %*% V
rownames(beta) <- vars$chr
beta
