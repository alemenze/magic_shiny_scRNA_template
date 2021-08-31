observe({
    updateSelectInput(session, "ISeuratObject", choices=c(names(data2$data)))
    updateSelectizeInput(session, "RClusterSet", choices=c(unique(levels(data2$data[[input$CCSeuratObject]]$seurat_clusters))), server=TRUE, options = list(maxOptions = 50))
})
#Idents
################################################################
output$reident <- renderRHandsontable({
    df=data.frame(ClusterNumber=unique(levels(data2$data[[input$ISeuratObject]]$seurat_clusters)), NewIdents='')

    table=rhandsontable(df, rowHeaderWidth=200) %>% hot_cols(colWidths=250) %>% hot_table(highlightCol = TRUE, highlightRow = TRUE, colHeaders = NULL)
    table
})

observeEvent(input$submitIdents, {
    newIdents <- hot_to_r(input$reident)[['NewIdents']]

    Idents(data2$data[[input$ISeuratObject]]) <- 'seurat_clusters'
    new.cluster.ids <- c(newIdents)
    names(new.cluster.ids) <- levels(data2$data[[input$ISeuratObject]])
    data2$data[[input$ISeuratObject]] <- RenameIdents(data2$data[[input$ISeuratObject]], new.cluster.ids)
    data2$data[[input$ISeuratObject]] <- StashIdent(data2$data[[input$ISeuratObject]], save.name=input$inputId)
    Idents(data2$data[[input$ISeuratObject]]) <- 'seurat_clusters'

    showNotification("New Identity added to the metadata.", type='message') 
})


# Reclustering
####################################################################
reclust_reactive <- reactive({
    withProgress(message='Calculating Re-clustering... this may take a bit',
        detail='Please stand by...',
        {
            tryCatch({
            shiny::setProgress(value=0.2, detail='Subsetting data...')
                reduced <- subset(data2$data[[input$ISeuratObject]], idents=c(input$RClusterSet))

            shiny::setProgress(value=0.3, detail='Digesting and splitting data...')
                DefaultAssay(reduced) <- 'RNA'
                reduced <- DietSeurat(reduced, assays='RNA')
                reduced.list <- SplitObject(reduced, split.by='samp')

            shiny::setProgress(value=0.5, detail='Re-normalizing...')
                for (i in 1:length(reduced.list)){
                    reduced.list[[i]] <- subset (reduced.list[[i]], subset= nFeature_RNA < input$RNFeatures & percent.mito < input$RPercentMito & percent.ribo < input$RPercentRibo & nFeature_RNA > input$RNFeaturesLow);
                    reduced.list[[i]] <- SCTransform(reduced.list[[i]], verbose=FALSE)
                }

            shiny::setProgress(value=0.7, detail='Integrating data...')
                reduced.features <- SelectIntegrationFeatures(object.list=reduced.list, nfeatures=3000);
                reduced.list <- PrepSCTIntegration(object.list=reduced.list, anchor.features=reduced.features, verbose=FALSE);
                reduced.anchors <- FindIntegrationAnchors(object.list=reduced.list, normalization.method='SCT', anchor.features=reduced.features, verbose=FALSE);
                reduced <- IntegrateData(anchorset=reduced.anchors, normalization.method='SCT');

            shiny::setProgress(value=0.9, detail='Clustering...')
                reduced <- RunPCA(reduced, verbose=FALSE, npcs=30);
                reduced <- RunUMAP(reduced, reduction='pca', dims=1:input$RUmapDim); 
                reduced <- FindNeighbors(reduced, reduction='pca', dims=1:input$RUmapDim); 
                reduced <- FindClusters(reduced,resolution=0.2);

            DefaultAssay(reduced) <- 'SCT'
            pos.markers <- FindAllMarkers(reduced, only.pos=TRUE, min.pct=0.25)
            dataM$markerslist[[paste(input$reclustId, 'Markers',sep='')]] <- pos.markers

            data2$data[[input$reclustId]] <- reduced

            rm(reduced)
            rm(pos.markers)
            rm(reduced.anchors)
            rm(reduced.list)
            rm(reduced.features)

            showNotification("Reclustering Completed- please see the associated tabs and select this object.", type='message', duration=NULL) 
                          
            },error=function(e)
            {
                status=paste("Re-clustering Error: ", e$message)
                showNotification(id="errorNotify", status, type='error', duration=NULL)
            })               
        }
    )

})

observeEvent(input$submitReclust, {
    reclust_reactive()
})