library(readr)
library(stringr)
library(dplyr)
library(edgeR)
library(tibble)


# Alignment --------------------------------------------------------------------

## fastQC quality check
# fastqc --extract --outdir=fastqc-reports *.fq

## alignment
# $ hisat2 -q -x /home/seq/resources/mm10/genome/genome -1 FCx1_5_1.fq -2 FCx1_5_2.fq -S FCx1_5.sam -p 8

## conversion sam to bam
# $ samtools view -bS FCx11_15.sam > FCx11_15.bam

## bam sorting
# samtools sort -m 2G -@ 8 <in.bam> <out.prefix>

## indexing
# samtools index <in_sorted.bam>

## counting reads in features
# ftp://ftp.ensembl.org/pub/release-90/gff3/mus_musculus/
# featureCounts -p -B -t gene -T 4 -a Mus_musculus.GRCm38.90.gff3 -o counts/FCx1_5_count.sub FCx11_15.sam >fc.out 2>&1


# NETKO Sequence ---------------------------------------------------------------

# 1. gtf reference for igv
# 2. knock out gene (fasta) -> hisat2-build -> alignment


# Analysis edgeR -------------------------------------------------------------

## reading the counts from a file
files <- str_subset(dir("data/counts"), 'FC.*_count.sub')
labels <- str_match(files, '(.*)_count.sub')[, 2]
# group - grouping factor, class: DGEList
FC_dge <- readDGE(files, "data/counts", group=c(2,1,2,1), labels=labels, columns=c(1, 7), header=T, skip=1)

## filtering out lowly expressed genes
FC_keep <- rowSums(cpm(FC_dge) > 1) >= 2
FC_dge <- FC_dge[FC_keep, , keep.lib.sizes=F]

## statystyka
FC_dge <- estimateDisp(FC_dge)
FC_lrt <- exactTest(FC_dge, pair=c(1,2))

gtf_annot <- read_tsv("data/counts/FCx1_5_count.sub", comment="#")[, 1:4]  # inne pliki takie same
wh_FC <- match(rownames(FC_lrt$table), gtf_annot$Geneid)
gtf_annot_FC <- gtf_annot[wh_FC, c("Chr", "Start", "End")]

FC <- FC_lrt$table %>% 
    as_tibble %>% 
    rownames_to_column(var="name") %>% 
    select(name, PValue) %>% 
    bind_cols(as_tibble(FC_dge$counts)[,c(2,1,4,3)]) %>% 
    add_column(chr=gtf_annot_FC$Chr, start=gtf_annot_FC$Start, end=gtf_annot_FC$End, .after="name") %>%
    {add_column(., fdr=p.adjust(unlist(.$PValue), method="fdr"), .after="PValue")} %>% 
    add_column(mean_wt=apply(FC_dge$counts, 1, function(x) mean(x[c("FCx1_5","FCx6_10")]))) %>% 
    add_column(mean_ko=apply(FC_dge$counts, 1, function(x) mean(x[c("FCx11_15","FCx16_20")]))) %>% 
    {add_column(., fold=apply(.[, c("mean_wt", "mean_ko")], 1, function(x) max(x[1]/x[2], x[2]/x[1])))} %>% 
    {add_column(., regulation=apply(.[, c("mean_wt", "mean_ko")], 1, function(x) if(diff(x) > 0) "up" else "down" ))} %>% 
    arrange(fdr)

write_excel_csv(FC, "FC.csv")

### LC
files <- str_subset(dir("data/counts"), 'LC.*_count.sub')
labels <- str_match(files, '(.*)_count.sub')[, 2]
LC_dge <- readDGE(files, "data/counts", group=c(2,1,2,1), labels=labels, columns=c(1, 7), header=T, skip=1)

LC_keep <- rowSums(cpm(LC_dge) > 1) >= 2
LC_dge <- LC_dge[LC_keep, , keep.lib.sizes=F]

LC_dge <- estimateDisp(LC_dge)
LC_lrt <- exactTest(LC_dge, pair=c(1,2))

wh_LC <- match(rownames(LC_lrt$table), gtf_annot$Geneid)
gtf_annot_LC <- gtf_annot[wh_LC, c("Chr", "Start", "End")]

LC <- LC_lrt$table %>% 
    as_tibble %>% 
    rownames_to_column(var="name") %>% 
    select(name, PValue) %>% 
    bind_cols(as_tibble(LC_dge$counts)[,c(2,1,4,3)]) %>% 
    add_column(chr=gtf_annot_LC$Chr, start=gtf_annot_LC$Start, end=gtf_annot_LC$End, .after="name") %>%
    {add_column(., fdr=p.adjust(unlist(.$PValue), method="fdr"), .after="PValue")} %>% 
    add_column(mean_wt=apply(LC_dge$counts, 1, function(x) mean(x[c("LC1_5","LC6_10")]))) %>% 
    add_column(mean_ko=apply(LC_dge$counts, 1, function(x) mean(x[c("LC11_15","LC16_20")]))) %>% 
    {add_column(., fold=apply(.[, c("mean_wt", "mean_ko")], 1, function (x) max(x[1]/x[2], x[2]/x[1])))} %>% 
    {add_column(., regulation=apply(.[, c("mean_wt", "mean_ko")], 1, function(x) if(diff(x) > 0) "up" else "down" ))} %>% 
    arrange(fdr)

write_excel_csv(LC, "LC.csv")


## LC, FC together
LCFC <- LC %>% bind_rows(FC)
LCFC <- LCFC[duplicated(LCFC$name) | duplicated(LCFC$name, fromLast=T),] %>% 
    arrange(name) %>% 
    select(1:10, 15:18, 11:14)
write_excel_csv(LCFC, "LCFC.csv")


loc <- str_match(chr8, "\t(\\d+)\t(\\d+)")
wh <- as.numeric(loc[, 2]) >= 92961047 & as.numeric(loc[, 3]) <= 93001667
str_subset(chr8[wh], "Slc")
