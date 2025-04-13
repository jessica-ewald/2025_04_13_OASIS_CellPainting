require(dplyr)
require(arrow)
require(ggplot2)
require(ggforce)

if (!requireNamespace("fastbmdR", quietly = TRUE)) {
  
  # Check if remotes is installed, and install it if not
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  # Install fastbmdR from GitHub
  remotes::install_github("jessica-ewald/fastbmdR@v0.0.0.9000")
}
library(fastbmdR)

#### 0. Define helper functions for plot fits
Exp2 <- function(b,c,d,e,f,dose){
  return(e * exp(b * dose))
}

Exp3 <- function(b,c,d,e,f,dose){
  return(e * (exp(sign(b) * (abs(b) * dose)^d)))
}

Exp4 <- function(b,c,d,e,f,dose){
  return(e * (c - (c - 1) * exp((-1) * b * dose)))
}

Exp5 <- function(b,c,d,e,f,dose){
  return(e * (c - (c - 1) * exp((-1) * (b * dose)^d)))
}

Hill <- function(b,c,d,e,f,dose){
  return(c + (d - c)/(1 + (dose/e)^b))
}

Pow <- function(b,c,d,e,f,dose){
  return(e + b * (dose^c))
}

Poly2 <- function(b,c,d,e,f,dose){
  return(b + c * (dose) + d * (dose^2))
}

Lin <- function(b,c,d,e,f,dose){
  return(d + b * (dose))
}


#### 1. Read in data and set common parameters
args <- commandArgs(trailingOnly = TRUE)

cp_pod_path <- args[1]
cc_pod_path <- args[2]
dist_path <- args[3]
cp_plot_path <- args[4]

cp_pods <- read_parquet(cp_pod_path) %>% as.data.frame()
cc_pods <- read_parquet(cc_pod_path) %>% as.data.frame()
dat <- read_parquet(dist_path) %>% as.data.frame()

highest_dose <- max(unique(dat$Metadata_Log10Conc))
highest_dose <- round(highest_dose + (0.025 * highest_dose), 1)
dose <- seq(0, highest_dose, by = 0.1)

n_cols <- 5
n_rows <- 5
n_per_page <- n_cols * n_rows
pdf_w <- 12
pdf_h <- 10


#### 2. Make cell painting plots
cp_compounds <- unique(cp_pods$Metadata_Compound)
plot_results <- data.frame(Metadata_Log10Conc = dose)
fitted_feats <- data.frame()
for (compound in cp_compounds){
  temp_pod <- cp_pods[cp_pods$Metadata_Compound == compound, ] %>% distinct()
  feat_type <- temp_pod$gene.id
  model <- temp_pod$mod.name
  b <- temp_pod$b
  c <- temp_pod$c
  d <- temp_pod$d
  e <- temp_pod$e
  f <- temp_pod$f
  
  f_dose <- switch(model,
                   "Exp2" = Exp2(b, c, d, e, f, dose),
                   "Exp3" = Exp3(b, c, d, e, f, dose),
                   "Exp4" = Exp4(b, c, d, e, f, dose),
                   "Exp5" = Exp5(b, c, d, e, f, dose),
                   "Hill" = Hill(b, c, d, e, f, dose),
                   "Pow" = Pow(b, c, d, e, f, dose),
                   "Poly2" = Poly2(b, c, d, e, f, dose),
                   "Lin" = Lin(b, c, d, e, f, dose))
  
  plot_results[, compound] <- f_dose
  
  # get observations
  dat_temp <- dat[dat$Metadata_Compound == compound, ]
  if (grepl("_ap", cp_pod_path)) {
    cmpd_ctrls <- grepl(compound, dat$Metadata_Compound) & grepl("DMSO", dat$Metadata_Compound)
    temp_dmso <- dat[cmpd_ctrls, ]
    temp_dmso <- temp_dmso[sample(1:720, 20), c(feat_type, "Metadata_Log10Conc")]
  } else {
    temp_plates <- unique(dat_temp$Metadata_Plate)
    temp_dmso <- dat[(dat$Metadata_Compound == "DMSO") & (dat$Metadata_Plate %in% temp_plates), 
                    c(feat_type, "Metadata_Log10Conc")]
  }

  feats <- dat_temp[, c(feat_type, "Metadata_Log10Conc")]
  feats <- rbind(temp_dmso, feats)

  colnames(feats)[1] <- "Observations"
  feats$Metadata_Compound <- compound
  feats$Observation_Type <- feat_type
  fitted_feats <- rbind(fitted_feats, feats)
}

plot_results <- reshape2::melt(plot_results, id.vars = c("Metadata_Log10Conc"))
colnames(plot_results)[2:3] <- c("Metadata_Compound", "f_dose")

# merge dataframe together
fitted_feats$Metadata_Log10Conc <- round(fitted_feats$Metadata_Log10Conc, 1)
plot_results <- merge(plot_results, fitted_feats, by=c("Metadata_Compound", "Metadata_Log10Conc"), all.x = TRUE, all.y = FALSE)

# Add better plotting label
cmpd_type <- plot_results[,c("Metadata_Compound", "Observation_Type")] %>% na.omit() %>% distinct()
cmpd_type$Compound_Observation <- paste0(cmpd_type$Metadata_Compound, " (", cmpd_type$Observation_Type, ")")
cmpd_type <- cmpd_type[,c("Metadata_Compound", "Compound_Observation")]
plot_results <- merge(plot_results, cmpd_type, by="Metadata_Compound")

# Add pods
cp_pods <- cp_pods[, c("Metadata_Compound", "bmdl", "bmd", "bmdu", "cc_POD")]
cp_pods$cc_POD[cp_pods$cc_POD == 9999] <- NA
cp_pods$bmdl <- round(cp_pods$bmdl, 1)
cp_pods$bmd <- round(cp_pods$bmd, 1)
cp_pods$bmdu <- round(cp_pods$bmdu, 1)
cp_pods$cc_POD <- round(cp_pods$cc_POD, 1)
plot_results <- merge(plot_results, cp_pods, by = "Metadata_Compound", all.x = TRUE, all.y = FALSE)

# Remove duplicated bmd values
for (compound in cp_compounds){
  inds <- which(plot_results$Metadata_Compound == compound)
  plot_results$bmdl[inds[-1]] <- NA
  plot_results$bmd[inds[-1]] <- NA
  plot_results$bmdu[inds[-1]] <- NA
  plot_results$cc_POD[inds[-1]] <- NA
}

plot_compounds <- unique(plot_results$Metadata_Compound)
n_pages <- ceiling(length(plot_compounds) / n_per_page)
plot_results <- plot_results[order(plot_results$Metadata_Compound), ]
pdf(cp_plot_path, width = pdf_w, height = pdf_h)
for (i in 1:n_pages) {
tryCatch({
    p <- ggplot(plot_results, aes(x = Metadata_Log10Conc)) +

      geom_point(aes(y = Observations)) +
      geom_line(aes(y = f_dose)) +

      geom_vline(aes(xintercept = bmdl), linetype = "dashed", color = "red", na.rm = TRUE) +
      geom_vline(aes(xintercept = bmd), linetype = "solid", color = "red", na.rm = TRUE) +
      geom_vline(aes(xintercept = bmdu), linetype = "dashed", color = "red", na.rm = TRUE) +
      geom_vline(aes(xintercept = cc_POD), linetype = "solid", color = "blue", na.rm = TRUE) +

      xlim(0, highest_dose) +
      facet_wrap_paginate(~ Compound_Observation, ncol = n_cols, nrow = n_rows, page = i, scales = "free_y") +
      theme_bw() +
      theme(strip.text = element_text(size = 6))
    print(p)
  }, error = function(e) {
    # Handle error: print a message and skip this plot
    message(sprintf("Error in plotting page %d: %s", i, e$message))
  })
}
dev.off()