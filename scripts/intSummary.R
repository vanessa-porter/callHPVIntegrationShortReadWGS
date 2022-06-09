#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla
.libPaths("/projects/vporter_prj/R/x86_64-centos7-linux-gnu-library/4.0")

## ---------------------------------------------------------------------------
## HPV Integration in Illumina Data
## Vanessa Porter, June 2022
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))

## ---------------------------------------------------------------------------
## LOAD INPUT 
## ---------------------------------------------------------------------------

# Make help options
option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL,
              help="Filtered SV VCF file", metavar="character"),
  make_option(c("-d", "--depth"), type="character", default=NULL,
              help="Samtools depth file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = "integration_sites.txt",
              help="Output file name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir

## read in the files
#vcf <- read.delim("/projects/hpv_nanopore_prj/htmcp/illumina_call_integration/output/HTMCP-03-06-02058/int/hpvIntFilt.vcf", 
#                  comment.char = "#", header = F)
#depth <- read.delim("/projects/hpv_nanopore_prj/htmcp/illumina_call_integration/output/HTMCP-03-06-02058/int/hpvSitesDepth.txt", 
#                    header = F)
vcf <- read.delim(opt$vcf, comment.char = "#", header = F)
depth <- read.delim(opt$depth, header = F)

## ---------------------------------------------------------------------------
## Organize data 
## ---------------------------------------------------------------------------

# extract # of supporting reads from vcf info
info <- strsplit(vcf$V8, ";")
count <- lapply(info, function(x){grep("PAIR_COUNT", x, value = T)[2]})
count <- unlist(count)
count <- gsub("PAIR_COUNT=", "", count)
count <- as.numeric(count)

# clean up the HPV sites
hpv <- vcf$V5
hpv <- gsub("A|C|G|T", "", hpv)
hpv <- gsub('\\[|\\]', "", hpv)
hpvSites <- as.data.frame(matrix(data = unlist(strsplit(hpv, ":")), ncol = 2, byrow = T))

# bring together in a dataframe
df <- data.frame(chr=vcf$V1, pos=vcf$V2, hpv.chr=hpvSites$V1, hpv.pos=hpvSites$V2,
           nreads=count, depth=depth$V3)

df$VAF <- df$nreads / df$depth
bed <- df[,c(1,2)]
bed$end <- bed$pos + 1

## ---------------------------------------------------------------------------
## Save files
## ---------------------------------------------------------------------------

write.table(df, file = paste0(out, "/integration_summary.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(bed, file = paste0(out, "/integration_sites.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
