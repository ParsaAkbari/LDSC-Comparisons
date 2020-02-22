#!/software/R-3.3.0/bin/Rscript
source(paste0(Sys.getenv('SCRIPT_DIR'), "/general_functions.R"))
# list of all new variants in their clumps
OUT_DIR <- Sys.getenv('OUT_DIR')
new_var_clumped <- read.table(paste0(OUT_DIR, "/condout/ldclump/hg19_clump_table.tsv"), header=T, stringsAsFactors=F)
new_var_clumped$novel <- rep("NA", nrow(new_var_clumped))

# list of cell paper variants
old.vars <- read.table(paste0(OUT_DIR, "/condout/ldclump_compare/cell_paper_variants.tsv"), header=T, stringsAsFactors=F)
traits <- as.vector(read.table(Sys.getenv('PHENO_NAMES'), header=F, stringsAsFactors=F)$V1)
traits = filter_traits(traits)
print(length(traits))
## extract all vars which aren't condsigs
trait_hg19s = c()
for (trait in traits) {
  trait_hg19s = unique(c(trait_hg19s, get_trait_hg19s(trait)))
}
# convert old vars to hg19
old.vars$chr <- as.numeric(gsub("(.+):.+_.+_.+","\\1", old.vars$VARIANT))
old.vars$bp <- gsub(".+:(.+)_.+_.+","\\1", old.vars$VARIANT)
old.vars$hg19 <- paste0(old.vars$chr, ":", old.vars$bp)
# read plink.ld file and create hg19 ids for each pair
plink.ld <- read.table(paste0(OUT_DIR, "/condout/ldclump_compare/plink.ld"), header=T, stringsAsFactors=F)
plink.ld$hg19A <- paste0(plink.ld$CHR_A, ":", plink.ld$BP_A)
plink.ld$hg19B <- paste0(plink.ld$CHR_B, ":", plink.ld$BP_B)
plink.ld$R2 <- as.numeric(plink.ld$R2)
plink.ld = plink.ld[plink.ld$R2 > 0.8,]
# loop through all the new clumps, do any of them have a variant in ld with old variant?
novel <- c()
not.novel <- c()
for (clump in unique(new_var_clumped$clump_id)) {
  print(clump)
  clump_novel <- 1
  # get all the variants in this clump
  new.clump.vars <- new_var_clumped$hg19[new_var_clumped$clump_id == clump]
  # loop through all the new variants in the clump
  for (new.var in new.clump.vars) {
    # find the associated set of rows with pairwise ld for this SNP
    partnersA <- plink.ld[plink.ld$hg19A == new.var & plink.ld$R2 > 0.8,]$hg19B
    partnersB <- plink.ld[plink.ld$hg19B == new.var & plink.ld$R2 > 0.8,]$hg19A
    partners <- c(partnersA, partnersB, new.var)
     if (sum(partners %in% old.vars$hg19) > 0) {
      clump_novel <- 0
      not.novel <- c(not.novel, clump)
      new_var_clumped[new_var_clumped$clump_id == clump,]$novel <- 0
      break
    }
  }
  if (clump_novel == 1) {
    novel <- c(novel, clump)
    new_var_clumped[new_var_clumped$clump_id == clump,]$novel <- 1
  }
}
print(paste("There are", length(novel), "novel clumps"))
print(paste("And", length(not.novel), "not novel clumps"))
print(paste("Total of", length(novel) + length(not.novel), "which equals total (sanity check)", length(unique(new_var_clumped$clump_id))))

n <- paste("There are", length(novel), "novel clumps \n And", length(not.novel), "not novel clumps \n Total of", length(novel) + length(not.novel), "which equals total (sanity check)", length(unique(new_var_clumped$clump_id)))
write(n, paste0(OUT_DIR, "/condout/ldclump_compare/written_results.txt"))

write.table(new_var_clumped, paste0(OUT_DIR, "/condout/ldclump_compare/hg19_clump_table_novel.tsv"), quote=F, row.names=F)
