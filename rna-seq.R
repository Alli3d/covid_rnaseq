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

count_data = tximport(files = quantsf_files,
         type = "salmon",
         tx2gene = gene_map,
         ignoreTxVersion = TRUE)


#formatting sample table to be useful
sample_table = as.data.frame(sample_table)
colnames(sample_table)[1] = "Sample"
conditions = c('mock_nhbe', 'infected_nhbe',
               'mock_a549', 'infected_a549')
conditions = rep(conditions, each=3)
conditions = conditions[0:10]
conditions = factor(conditions)
sample_table$conditions = conditions

deseq_dataset = DESeqDataSetFromTximport(txi = count_data,
                                         colData = sample_table,
                                         design =~conditions)

deseq_dataset = estimateSizeFactors(deseq_dataset)
counts(deseq_dataset, normalized=TRUE)[1:6, 1:3]
vst = varianceStabilizingTransformation(deseq_dataset)

plotPCA(vst, intgroup='conditions') +
  theme_bw()



# splitting data
dds1 = deseq_dataset[,1:6]
dds2 = deseq_dataset[,7:10]

dds1 = estimateSizeFactors(dds1)
dds2 = estimateSizeFactors(dds2)

vst1 = varianceStabilizingTransformation(dds1)
plotPCA(vst1, intgroup='conditions')+
  theme_bw()

vst2 = varianceStabilizingTransformation(dds2)
plotPCA(vst2, intgroup='conditions')+
  theme_bw()


sample_files_nhbe = quantsf_files[1:6]
count_data_nhbe = tximport(files = sample_files_nhbe,
                      type = "salmon",
                      tx2gene = gene_map,
                      ignoreTxVersion = TRUE)
sample_table_nhbe = sample_table[1:6,]
sample_table_nhbe$conditions = factor(rep(c('mock', 'infected'), each = 3),
                                 levels=c('mock', 'infected'))
dds_nhbe = DESeqDataSetFromTximport(txi = count_data_nhbe,
                                         colData = sample_table_nhbe,
                                         design =~conditions)

dds_nhbe = estimateSizeFactors(dds_nhbe)
dds_nhbe = estimateDispersions(dds_nhbe)
plotDispEsts(dds_nhbe)

dds_nhbe = nbinomWaldTest(dds_nhbe)

result_table = results(dds_nhbe)
summary(result_table)

result_df = as.data.frame(result_table)
View(result_df)

plotCounts(dds_nhbe, gene='ENSG00000265794', intgroup = 'conditions')

filter_df1 = result_df[complete.cases(result_df),]
filter_df1 = rownames_to_column(filter_df1, var='ensgene')


# filter results by p value
# padj < 0.05
# log2FC > 1 or < -1

filter_df2 = filter_df1[filter_df1$padj < 0.05, ]
filter_df3 = filter_df2[abs(filter_df2$log2FoldChange) > 1, ]

View(filter_df3)


# plotting!!
plotMA(result_table)

filter_df1$dif_expressed = filter_df1$padj < 0.05 & abs(filter_df1$log2FoldChange) > 1



ensembl99 = useEnsembl(biomart="ensembl", version=99)
View(listDatasets(ensembl99))
ensembl99 = useDataset("hsapiens_gene_ensembl", mart=ensembl99)
View(listAttributes(ensembl99))
View(listFilters(ensembl99))
getBM(attributes=c('ensembl_gene_id','ensembl_gene_id_version',
                   'ensembl_transcript_id','ensembl_transcript_id_version',
                   'external_gene_name'),
                 filters = c('ensembl_gene_id'),
                 values = filter_df1$ensgene[1:6],
                 mart = ensembl99)

annotation = getBM(attributes=c('ensembl_gene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                'strand',
                                'gene_biotype',
                                'external_gene_name',
                                'description'),
                   filters = c('ensembl_gene_id'),
                   values = filter_df1$ensgene,
                   mart = ensembl99)


annotated_df = right_join(filter_df1, annotation,
                         by=c('ensgene'='ensembl_gene_id'))
g = ggplot(annotated_df, aes(x=log2FoldChange,
                           y=-log10(padj),
                           name=external_gene_name)) +
  geom_point(aes(colour=dif_expressed)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_hline(yintercept = -log10(0.05))

ggplotly(g)

anno_df2 = annotated_df[annotated_df$padj < 0.05, ]
anno_df3 = anno_df2[abs(anno_df2$log2FoldChange) > 1, ]

degs = anno_df3$ensgene

vst_nhbe = varianceStabilizingTransformation(dds_nhbe)
vst_nhbe_matrix = assay(vst_nhbe)

data_for_hm = vst_nhbe_matrix[degs, ]
rownames(data_for_hm) = anno_df3$external_gene_name

pheatmap(data_for_hm, fontsize_row=4, cutree_cols=2, cutree_rows = 2)


## testing for go enrichment

ent_gene = getBM(attributes=c('entrezgene_id'),
                   filters = c('ensembl_gene_id'),
                   values = anno_df3$ensgene,
                   mart = ensembl99)
ent_gene = as.character(ent_gene$entrezgene_id)

ent_uni = getBM(attributes=c('entrezgene_id'),
                 filters = c('ensembl_gene_id'),
                 values = annotated_df$ensgene,
                 mart = ensembl99)
ent_uni = as.character(ent_uni$entrezgene_id)

ego = enrichGO(gene = ent_gene,
               OrgDb = org.Hs.eg.db ,
               ont = 'BP',
               universe = ent_uni)

barplot(ego)
dotplot(ego)
cnetplot(ego)

write_tsv(anno_df3, "filtered_nhbe_results.txt")

