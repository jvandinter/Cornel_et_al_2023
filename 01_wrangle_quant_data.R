# Recreation of Figure for Annelisa
library(tidyverse)

project_dir = "/hpc/pmc_vanheesch/projects/Jip/counts"

# Read cohort metadata
meta <-
  read.delim(paste(project_dir, "cohort.txt", sep = "/"))

# Grab MR1 count info
mr1l <- pbapply::pbapply(meta, 1, function(x) {
  mr1 <-
    data.table::fread(paste(project_dir, x[1], sep = "/"),
                      sep = "\t",
                      skip = 3)[c(4094, 19125, 19222, 19223),]
  mr1$Biomaterial.Id <- x[3]
  return(mr1)
})

# Combine in single DF
mr1 <- do.call(plyr::rbind.fill, mr1l)

# Some files have a different annotation for MR1, find those, subset the metadata
# and rerun those with the correct location
table(mr1$GeneName %in% c("MR1", "HLA-A", "HLA-B", "HLA-C"))
to_rerun <-
  subset(meta, Biomaterial.Id %in% mr1[which(!(mr1$GeneName %in% c("MR1", "HLA-A", "HLA-B", "HLA-C"))),]$Biomaterial.Id)
mr1_correct <-
  mr1[which(mr1$GeneName %in% c("MR1", "HLA-A", "HLA-B", "HLA-C")),]

test <- data.table::fread(
  paste(project_dir,
        to_rerun$quant_files[1],
        sep = "/"),
  sep = "\t",
  skip = 3
)

# Rerun old version samples
mr1l <- pbapply::pbapply(to_rerun, 1, function(x) {
  mr1 <-
    data.table::fread(paste(project_dir, x[1], sep = "/"),
                      sep = "\t",
                      skip = 3)[c(4016, 18414, 18511, 18510),]
  mr1$Biomaterial.Id <- x[3]
  return(mr1)
})

# Make sure they have the same columns
mr1_rerun <- do.call(plyr::rbind.fill, mr1l)
mr1_rerun$`# ID` <- NA

mr1_full <- rbind(mr1_correct, mr1_rerun)

# Subset on interesting columns
mr1_full <-
  mr1_full[, c("Counts", "CPM", "FPKM", "GeneID", "GeneName", "Biomaterial.Id")]

# Connect to interesting metadata
plotdata <-
  dplyr::left_join(mr1_full, meta[, c("Biomaterial.Id", "V2")])
unique(plotdata$V2)

# Edit sample labels for plotting
plotdata$V2 <- ifelse(
  plotdata$V2 == "Hepatoblastoma",
  "HBL",
  ifelse(
    plotdata$V2 == "Medulloblastoma",
    "MB",
    ifelse(plotdata$V2 == "Ependymoma", "EPN", plotdata$V2)
  )
)
plotdata$embryonal <-
  ifelse(plotdata$V2 %in% c("NBL", "MRT", "HBL", "ATRT", "MB"),
         "Embryonal",
         "Not Embryonal")

plotdata_wide <-
  data.frame(Biomaterial.Id = unique(plotdata$Biomaterial.Id)) %>%
  dplyr::left_join(plotdata[, 6:8]) %>%
  dplyr::distinct(Biomaterial.Id, GeneName, .keep_all = TRUE) %>%
  dplyr::mutate(
    MR1 =
      plotdata[which(plotdata$GeneName == "MR1"), c("CPM")],
    HLA_A =
      plotdata[which(plotdata$GeneName == "HLA-A"), c("CPM")],
    HLA_B =
      plotdata[which(plotdata$GeneName == "HLA-B"), c("CPM")],
    HLA_C =
      plotdata[which(plotdata$GeneName == "HLA-C"), c("CPM")]
  )

# Create list of samples
write.table(
  as.data.frame(table(plotdata$V2)),
  file = paste(project_dir, "samples_per_type.txt", sep = "/"),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

write.table(
  plotdata,
  file = paste(project_dir, "MR1_plotdata.txt", sep = "/"),
  row.names = F,
  quote = F,
  sep = "\t"
)

write.table(
  plotdata_wide,
  file = paste(project_dir, "mr1_hla_plotdata.txt", sep = "/"),
  row.names = F,
  quote = F,
  sep = "\t"
)
