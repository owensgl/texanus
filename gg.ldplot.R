library(devtools)
install_github("timknut/gg.ldplot")
library(gg.ldplot)
ld <- read.table("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.gen1.geno.ld", header = TRUE, stringsAsFactors = F)
ld %>% filter(CHR == "HanXRQChr00c0041") -> tmp.ld
colnames(tmp.ld) <- c("CHR","POS1","POS2","N_CHR","R.2")
plot_LDmatrix(tmp.ld)
x = tmp.ld

assert_that(is.data.frame(x))
chrom <- unique(x$CHR)
new_df <- dplyr::data_frame(
  CHR = chrom,
  POS1 = unique(x$POS1),
  POS2 = unique(x$POS1),
  N_CHR = unique(x$N_CHR),
  R.2 = 1,
  D = NA,
  Dprime = NA
)
names(new_df)[5] <- "R.2"
names(x)[5] <- "R.2"
new_df <- dplyr::bind_rows(x, new_df)
new_df <- dplyr::arrange(new_df, POS1, POS2)
new_df