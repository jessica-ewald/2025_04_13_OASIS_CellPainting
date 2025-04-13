require(dplyr)
require(arrow)
require(data.table)


######## 1. Parse arguments and define parameters
args <- commandArgs(trailingOnly = TRUE)

bmd_path <- args[1]
cc_pod_path <- args[2]
pod_path <- args[3]


######## 2. Select POD as minimum BMD across gmd and all categories
bmd <- read_parquet(bmd_path) %>% as.data.frame()

# This is what the SD of residuals must be less than
bmd$res.thresh <- 3 * bmd$SDctrl
bmd <- as.data.table(bmd)

bmd <- bmd[all.pass == TRUE]
bmd <- bmd[SDres < res.thresh] # Filter by the SD of the residual

min_bmd <- bmd[,
               .SD[bmd == min(bmd)],
               by = Metadata_Compound,
               .SDcols = c("gene.id", "mod.name", "bmd",
                           "b", "c", "d", "e", "f",
                           "bmdl", "bmdu", "bmr", "lof.p")]

min_bmd <- as.data.frame(min_bmd)


######## 3. Read in CC PODs
cc <- read_parquet(cc_pod_path) %>% as.data.frame()
cc$bmd[cc$all.pass != TRUE] <- 9999
cc <- cc[, c("Metadata_Compound", "bmd")]
colnames(cc)[2] <- "cc_POD"

min_bmd <- merge(min_bmd, cc, by = "Metadata_Compound")
min_bmd$PAC_below_cc_POD <- min_bmd$bmd < min_bmd$cc_POD


######## 4. Write out results
write_parquet(min_bmd, pod_path)
