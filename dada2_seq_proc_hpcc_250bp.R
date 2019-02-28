# salloc -N 1 -c 20 --time=1:00:000 --mem-per-cpu=3000M

library(dada2);library(data.table)
set.seed(1390)
CORES = parallel::detectCores()
path <- "/mnt/research/germs"
#maximum_read_length = 150

fnFs <- sort(list.files(path, pattern="_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1_"), `[`, 1)
number_samples <- length(sample.names)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#trunclen is the number of bps to trim the R1 and R2 reads to.. typically based on the read quality report
out <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, truncLen=c(240,240), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=CORES) # On Windows set multithread=FALSE

derepFs <- derepFastq(filtFs, verbose=FALSE);names(derepFs) <- sample.names
derepRs <- derepFastq(filtRs, verbose=FALSE);names(derepRs) <- sample.names

errF <- learnErrors(filtFs, multithread=CORES, randomize = TRUE)
errR <- learnErrors(filtRs, multithread=CORES, randomize = TRUE)

dadaFs <- dada(derepFs, err=errF, multithread=CORES, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=CORES, pool=TRUE)
merged_reads <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)

seq_table <- makeSequenceTable(merged_reads)
seq_table <- removeBimeraDenovo(seq_table, method="consensus", multithread = CORES)
saveRDS(t(seq_table), paste0(path,"/../seq_table.RDS"))

tax_table <- assignTaxonomy(seq_table, "/mnt/research/germs/databases/dada2/rdp_train_set_16.fa.gz", multithread = CORES)
saveRDS(tax_table, paste0(path,"/../tax_table.RDS"))
