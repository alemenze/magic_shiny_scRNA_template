# Set up
## Loading the libraries of use
```
library(Seurat);library(dplyr);library(hdf5r);library(cowplot);
```

## Load a metadata file noting the different samples and create a list
```
opt <- read.csv(file='metadata.csv')
samplelist <- list()
```

# Import each sample and process the same
```
cellcounts <- data.frame(matrix(ncol=2,nrow=0)); colnames(cellcounts) <- c('type','count')
for (i in 1:nrow(opt)){
    row <- opt[i,]
    filename <- row$samples
    samp <- row$name
    samp<-paste(samp)

    directory <- row$directory

    ###Import raw data
    data_in <- Read10X_h5(filename=paste(directory,filename,sep=''), use.names=TRUE, unique.features=TRUE);

    ###Create Seurat object
    s <- CreateSeuratObject(counts=data_in, project=,min.cells=3, min.features=200);

    ### Add metadata tags. You can add as many as needed from the metadata table
    s@meta.data$samp <- paste(samp);

    ### Mitochondrial
    s[['percent.mito']] <- PercentageFeatureSet(s, pattern='^mt-');
    ####Human:
    ####s[['percent.mito']] <- PercentageFeatureSet(s, pattern='^MT-');

    ### Ribosomal
    s[['percent.ribo']] <- PercentageFeatureSet(s,pattern='^Rp[sl][[:digit:]]');
    #####Human:
    ##### s[['percent.ribo']] <- PercentageFeatureSet(s,pattern='^RP[SL][[:digit:]]');  

    ## Visualize pre-filter

    qc_in_plot <- FeatureScatter(s, feature1='nCount_RNA',feature2='nFeature_RNA');
    assign(paste(samp,'_feature_in',sep=''),qc_in_plot);

    ### Violin/Dot plots

    qc_in_plot <- VlnPlot(s, features=c('nFeature_RNA','nCount_RNA','percent.mito','percent.ribo'),ncol=4);
    assign(paste(samp,'_vln_in',sep=''),qc_in_plot);


    ###Cell count pre-filter
    cellcounts[nrow(cellcounts)+1,] <- list(sample=paste(samp, '_prefilt',sep=''), count=length(WhichCells(s)))

    ## Filter
    s <- subset (s, subset= nFeature_RNA < 6000 & percent.mito < 10 & percent.ribo < 45 & nFeature_RNA > 1000 );

    ### Violin/Dot plots
    out_plot <- VlnPlot(s, features=c('nFeature_RNA','nCount_RNA','percent.mito','percent.ribo'),ncol=4);
    assign(paste(samp,'_vln_out',sep=''),out_plot);

    ### Scatter plots

    out_plot <- FeatureScatter(s, feature1='nCount_RNA',feature2='nFeature_RNA');

    ###Cell count post-filter
    cellcounts[nrow(cellcounts)+1,] <- list(sample=paste(samp, '_postfilt',sep=''), count=length(WhichCells(s)))

    ## Transform and regress the data
    s <- SCTransform(s, vars.to.regress=c('percent.mito','percent.ribo'));

    ## Rename temp files
    assign(paste(samp),s);

    samplelist <- c(samplelist, eval(parse(text=samp)))

}
```

# Perform the integration and clustering
```
combined.features <- SelectIntegrationFeatures(object.list=samplelist, nfeatures=3000);
combined.list <- PrepSCTIntegration(object.list=samplelist, anchor.features=combined.features, verbose=FALSE);
combined.anchors <- FindIntegrationAnchors(object.list=combined.list, normalization.method='SCT', anchor.features=combined.features, verbose=FALSE);
combined <- IntegrateData(anchorset=combined.anchors, normalization.method='SCT');

combined <- RunPCA(combined, verbose=FALSE, npcs=30);
combined <- RunUMAP(combined, reduction='pca', dims=1:15); 
combined <- FindNeighbors(combined, reduction='pca', dims=1:15); 
combined <- FindClusters(combined,resolution=0.5); 

```

# Find markers per cluster
```
DefaultAssay(combined) <- 'SCT'
clustnum <- nrow(table(Idents(combined),combined$samp))
for (i in 0:(clustnum-1)){
  run <- try({
    clust.markers <- FindMarkers(combined, ident.1=i, only.pos=TRUE, print.bar=FALSE)
    assign(paste('clust', i,'.markers', sep=''), clust.markers)
  })
  if(inherits(run, 'try-error')) {
    print(paste(i,' needs to be manually markered', sep=''))
  }
}

all.markers <- FindAllMarkers(combined, min.pct=0, logfc.threshold=0)
pos.markers <- FindAllMarkers(combined, only.pos=TRUE, min.pct=0.25)
```