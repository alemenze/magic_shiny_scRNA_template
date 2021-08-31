# Re-clustering

The reclustering functions akin to the original clustering. It is highly recommended to not tweak the default parameters without understanding the underlying parameters. Pre-designed reclusterings can be included in the base image upon request. 

Upon reclustering the object will be added to the object list and can be accessed in the other tabs. 

## Advanced Options
- nFeatures upper limit: To remove putative doublets. This cannot detect doublets in the middle of the range, but can help remove the more blatant ones.
- nFeatures lower limit: To remove low quality GEMs. Generally with the current V3/NextGEM kits this should be cut at 1000 as the lower limit to remove noise GEMs. Certain cell types might require a slight lowering of this value. 
- Mitochondrial percent cutoff: Also to remove low quality cells. Mitochondrial reads are a proxy for dead/dying cells. 
- Ribosomal percent cutoff: Can be used to help remove in poor library preps. Generally this can vary by cell type, with T-cells being on the higher end. 
- UMAP Dimensions: The principal components to include in the UMAP clustering and identification of nearest neighbors. This can change depending upon the complexity of the sample composition. 
- Clustering Resolution: Define the clarity of de novo cluster identification. 

