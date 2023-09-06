library(tidyverse)
library(VennDiagram)
library(pheatmap)
source('rna-seq.R')
source('covid-a549-analysis.R')



nhbe_ensegene = annotated_df$ensgene
a549_ensegene = a549_f1$ensgene

# common genes between a549 and nhbe
common_ensgene = intersect(nhbe_ensegene,a549_ensegene)

nhbe_common = annotated_df[annotated_df$ensgene %in% common_ensgene,
             c('ensgene','log2FoldChange')]
a549_common = a549_f1[a549_f1$ensgene %in% common_ensgene,
                        c('ensgene', 'log2FoldChange')]

common_fold_changes = left_join(nhbe_common, a549_common, by=c('ensgene'='ensgene'))

ggplot(data=common_fold_changes,
       aes(x =log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point(alpha=0.3) +
  geom_abline(slope=1, intercept=0) +
  xlim(min=-5, max=5)+
  ylim(min=-5, max=5)

cor.test(common_fold_changes$log2FoldChange.x,
         common_fold_changes$log2FoldChange.y,
         method= 'spearman')



nhbe_diff = anno_df2$ensgene
a549_diff = a549_anno2$ensgene

union_diff = union(nhbe_diff, a549_diff)
union_fc_hm = filter(common_fold_changes,
                     ensgene %in% union_diff)


