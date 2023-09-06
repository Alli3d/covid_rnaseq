# Load packages
library(tximport)
library(DESeq2)
library(plotly)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)

sample_table = unique(read.csv("SraRunTable.txt") %>% 
                        dplyr::select(Sample.Name, source_name, Treatment,
                                      Cell_Line, Cell_type, Time_point))




quantsf_files = paste0(pull(sample_table, Sample.Name), '/quant.sf')
names(quantsf_files) = pull(sample_table, Sample.Name)

gene_map = read.csv("gene_map.csv", col.names = c('enstid', 'ensgid'))




sample_files_a549 = quantsf_files[7:10]
count_data_a549 = tximport(files = sample_files_a549,
                           type = "salmon",
                           tx2gene = gene_map,
                           ignoreTxVersion = TRUE)
sample_table_a549 = sample_table[7:10,]
sample_table_a549$conditions = factor((c('mock', 'infected', 'infected', 'infected')),
                                      levels=c('mock', 'infected'))


dds_a549 = DESeqDataSetFromTximport(txi = count_data_a549,
                                    colData = sample_table_a549,
                                    design =~conditions)

dds_a549 = estimateSizeFactors(dds_a549)
dds_a549 = estimateDispersions(dds_a549)
plotDispEsts(dds_a549)

dds_a549 = nbinomWaldTest(dds_a549)

vst_a549 = varianceStabilizingTransformation(dds_a549)

plotPCA(vst_a549, intgroup='conditions')

a549_results = results(dds_a549)
summary(a549_results)

plotMA(a549_results)
vst_mat = assay(vst_a549)



a549_f1 = a549_results[complete.cases(a549_results),]
a549_f1 = as.data.frame(a549_f1)
a549_f1 = rownames_to_column(a549_f1, var='ensgene')

a549_f2 = a549_f1[a549_f1$padj < 0.05,]
a549_f3 = a549_f2[abs(a549_f2$log2FoldChange) > 1, ]

# annotation

ensembl99 = useEnsembl(biomart="ensembl", version=99)
ensembl99 = useDataset("hsapiens_gene_ensembl", mart=ensembl99)

a549_anno = getBM(attributes=c('ensembl_gene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                'strand',
                                'gene_biotype',
                                'external_gene_name',
                                'description'),
                   filters = c('ensembl_gene_id'),
                   values = a549_f1$ensgene,
                   mart = ensembl99)


annotated_a549df = right_join(a549_f1, a549_anno,
                          by=c('ensgene'='ensembl_gene_id'))

a549_anno2 = filter(annotated_a549df, padj < 0.05)
a549_anno3 = filter(a549_anno2, abs(log2FoldChange) > 1)

a549_degs = a549_anno3$ensgene
a549_hm = vst_mat[a549_degs,]
rownames(a549_hm) = a549_anno3$external_gene_name

heatmap(a549_hm)

pheatmap(a549_hm, fontsize_row=4, scale='row')

entrez_a549 = getBM(attributes=c('entrezgene_id'),
                    filters = c('ensembl_gene_id'),
                    values = a549_anno3$ensgene,
                    mart = ensembl99)
entrez_a549 = entrez_a549$entrezgene_id
entrez_a549 = as.character(entrez_a549)
a549_universe = getBM(attributes=c('entrezgene_id'),
                      filters = c('ensembl_gene_id'),
                      values = annotated_a549df$ensgene,
                      mart = ensembl99)
a549_universe = a549_universe$entrezgene_id
a549_universe = as.character(a549_universe)

ego_a549 = enrichGO(gene = entrez_a549,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    universe = a549_universe,
                    readable = TRUE)
barplot(ego_a549, showCategory=20)
dotplot(ego_a549, showCategory=20)

fold_changes = a549_anno3$log2FoldChange
names(fold_changes) = a549_anno3$external_gene_name

cnetplot(ego_a549, 
         showCategory=10, 
         foldChange = fold_changes)
goplot(ego_a549)
