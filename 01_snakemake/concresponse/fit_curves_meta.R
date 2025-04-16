require(dplyr)
require(arrow)
require(parallel)

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

dat_path <- args[1]
output_path <- args[2]
num_sds <- args[3]
meta_nm <- args[4]

ctrl <- "DMSO"


######## 2. Calculate BMDs from cell counts
dat <- read_parquet(dat_path) %>% as.data.frame()
dat <- dat[dat$Metadata_well_type != "JUMP_control", ]

compounds <- unique(dat$Metadata_Compound)
compounds <- compounds[compounds != ctrl]

cc_dmso <- dat[dat$Metadata_Compound == ctrl, ]
cc_comp <- dat[dat$Metadata_Compound != ctrl, ]

# Define compound-specific fitting function
fit_compound <- function(compound) {
  cc_comp_fit <- cc_comp[cc_comp$Metadata_Compound == compound, ]
  cmpd_plates <- cc_comp_fit$Metadata_Plate
  cc_dmso_fit <- cc_dmso[cc_dmso$Metadata_Plate %in% cmpd_plates, ]
  cc_fit <- rbind(cc_dmso_fit, cc_comp_fit)
  cc_fit <- cc_fit[order(cc_fit$Metadata_Concentration), ]

  cc <- matrix(cc_fit[, meta_nm], nrow = 1)
  rownames(cc) <- "cc"

  dose <- cc_fit$Metadata_Log10Conc

  cc_pod <- scoresPOD(cc, dose, log10.dose = TRUE, num.sds = num_sds)
  if (!is.null(cc_pod)) {
    cc_pod$Metadata_Compound <- compound
    return(cc_pod)
  } else {
    return(NULL)
  }
}

# Run in parallel
bmd_list <- mclapply(compounds, fit_compound, mc.cores = detectCores())
bmd_res <- do.call(rbind, bmd_list)



######## 3. Merge OASIS IDs and write output
meta_info <- distinct(dat[, c("Metadata_Compound", "Metadata_OASIS_ID")])
bmd_res <- merge(as.data.frame(bmd_res), meta_info, by = "Metadata_Compound")

write_parquet(bmd_res, output_path)