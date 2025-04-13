require(dplyr)
require(arrow)

######## 0. Make sure fastbmdR is installed
if (!requireNamespace("fastbmdR", quietly = TRUE)) {
  
  # Check if remotes is installed, and install it if not
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  # Install fastbmdR from GitHub
  remotes::install_github("jessica-ewald/fastbmdR@v0.0.0.9000")
}
library(fastbmdR)


######## 1. Parse arguments and define parameters
args <- commandArgs(trailingOnly = TRUE)

input_path <- args[1]
output_path <- args[2]
num_sds <- args[3]

ctrl <- "DMSO"

######## 2. Calculate BMDs from global mahalanobis distances
dat <- read_parquet(input_path) %>% as.data.frame()

compounds <- unique(dat$Metadata_Compound)
compounds <- compounds[!grepl(ctrl, compounds)]

dat_dmso <- dat[dat$Metadata_Compound == ctrl, ]
dat_comp <- dat[dat$Metadata_Compound != ctrl, ]

feat_cols <- colnames(dat)[!grepl("Metadata", colnames(dat))]

bmd_res <- data.frame()
for (compound in compounds){
  comp_fit <- dat[dat$Metadata_Compound == compound, ]

  if (grepl("_ap", input_path)) {
    cmpd_ctrls <- grepl(compound, dat$Metadata_Compound) & grepl("DMSO", dat$Metadata_Compound)
    dmso_fit <- dat[cmpd_ctrls, ]
  } else {
    cmpd_plates <- comp_fit$Metadata_Plate
    dmso_fit <- dat_dmso[dat_dmso$Metadata_Plate %in% cmpd_plates, ]
  }

  dat_fit <- rbind(dmso_fit, comp_fit)
  dat_fit <- dat_fit[order(dat_fit$Metadata_Concentration), ]

  dat_mat <- t(dat_fit[,feat_cols])
  rownames(dat_mat) <- feat_cols

  dose <- c(dat_fit$Metadata_Log10Conc)

  pod <- scoresPOD(dat_mat, dose, log10.dose = TRUE, num.sds = num_sds,
                   filt.var = "SDres")
  if (!is.null(pod)) {
    pod$Metadata_Compound <- compound
    bmd_res <- rbind(bmd_res, pod)
  }
}

write_parquet(as.data.frame(bmd_res), output_path)