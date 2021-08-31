# Column descriptors

## First column (unlabelled)
Gene symbol made unique per R. This will include the base gene symbol, but if it is repeated for different clusters it will add .X to it for each replicate. 

## Second column- p_val
P value based on Wilcoxon rank sum test

## Third column- avg_log2FC
The average fold change (log2 scale) between your comparators. 
For marker genes, this is the cells in the respective cluster vs all other cells. 

## Fourth column- pct.1
Percentage of cells in the first comparator (numerator) that contain the respective gene. For marker genes, this is the cells in cluster of note.  

## Fifth column- pct.2
Percentage of cells in the second comparator (denominator) that contain the respective gene. For marker genes, this is all cells not in the cluster of note

## Sixth column- p_val_adj
Bonferroni corrected p-value. When choosing a p-value, this is preferred as it reduces false positives. 

## Seventh column- cluster
The respective cluster for this row. 

## Last column- gene
The gene symbol for this row. 
