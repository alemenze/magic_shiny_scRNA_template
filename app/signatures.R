observe({
    updateSelectInput(session, "SSeuratObject", choices=c(names(data2$data)))
})

# Cell Cycle
#################################################################
observe({
    updateSelectInput(session, "CycleGroupBy", choices=c(colnames(data2$data[[input$SSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
    updateSelectInput(session, "CycleSplitBy", choices=c(colnames(data2$data[[input$SSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')   
})

cell_cycler <- reactive({
    withProgress(message='Performing signature enrichment',
        detail='Please stand by...',
        {
            shiny::setProgress(value=0.2, detail='Converting gene names')
            # If Mouse:
            #######################################
            convertHumanGeneList <- function(x){

            require("biomaRt")
            human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
            mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

            genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

            humanx <- unique(genesV2[, 2])

            # Print the first 6 genes found to the screen
            print(head(humanx))
            return(humanx)
            }

            shiny::setProgress(value=0.6, detail='Scoring the matrix')
            m.s.genes <- convertHumanGeneList(cc.genes$s.genes)
            m.g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)

            score <- CellCycleScoring(data2$data[[input$SSeuratObject]],s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

            shiny::setProgress(value=0.8, detail='Creating the plot')

            if(input$CyclePlotType=='Vln'){
                cyclerplot <- VlnPlot(score, features=c(input$CycleScore),
                    pt.size=input$CycleVlnPointSize,
                    group.by=input$CycleGroupBy,
                    split.by=input$CycleSplitBy
                )
            }
            if(input$CyclePlotType=='Feature'){
                cyclerplot <- DimPlot(score, pt.size=input$CycleFeaturePointSize)
            }
            return(cyclerplot)

        }
    )
})

observe({
    output$cycleplot <- renderPlot({
        cell_cycler()
    }, height=input$CHeight, width=input$CWidth)
})

output$DownloadCycle <- downloadHandler(
    filename=function(){
        paste('Cellcycle',input$DownCycleFormat,sep='.')
    },
    content=function(file){   
        if(input$DownCycleFormat=='jpeg'){
            jpeg(file, height=input$HHeight, width=input$HWidth)
            print(cell_cycler())
            dev.off()
        }
        if(input$DownCycleFormat=='png'){
            png(file, height=input$HHeight, width=input$HWidth)
            print(cell_cycler())
            dev.off()
        }
        if(input$DownCycleFormat=='tiff'){
            tiff(file, height=input$HHeight, width=input$HWidth)
            print(cell_cycler())
            dev.off()
        }
    }
)

# Gene Signature
#################################################################
observe({
    updateSelectInput(session, "SignGroupBy", choices=c(colnames(data2$data[[input$SSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
    updateSelectInput(session, "SignSplitBy", choices=c(colnames(data2$data[[input$SSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')

    temp <- msigdbr(species='mouse', category=input$SignSource)
    updateSelectInput(session, "SignPath", choices=unique(temp$gs_name))   
})


signature_plotter <- reactive({
    withProgress(message='Performing gene set signature enrichment',
        detail='Please stand by...',
        {
            shiny::setProgress(value=0.4, detail='Scoring the matrix')
            temp <- msigdbr(species='mouse', category=input$SignSource)
            genenames <- subset(temp, gs_name==input$SignPath)

            tryCatch({
                score <- AddModuleScore(data2$data[[input$SSeuratObject]],features=c(genenames$gene_symbol), name='Module', ctrl=10)

                shiny::setProgress(value=0.8, detail='Creating the plot')

                signplot <- VlnPlot(score, features='Module1',
                    pt.size=input$SignVlnPointSize,
                    group.by=input$SignGroupBy,
                    split.by=input$SignSplitBy
                )
                return(signplot)  
            },error=function(e)
            {
                status=paste("Module Score Error: ", e$message)
                showNotification(id="errorNotify", status, type='error', duration=NULL)
            })               
        }
    )
})

observe({
    output$signatureplot <- renderPlot({
        signature_plotter()
    }, height=input$CHeight, width=input$CWidth)
})

output$DownloadSignature <- downloadHandler(
    filename=function(){
        paste('Signature',input$DownSignatureFormat,sep='.')
    },
    content=function(file){   
        if(input$DownSignatureFormat=='jpeg'){
            jpeg(file, height=input$CHeight, width=input$CWidth)
            print(signature_plotter())
            dev.off()
        }
        if(input$DownSignatureFormat=='png'){
            png(file, height=input$CHeight, width=input$CWidth)
            print(signature_plotter())
            dev.off()
        }
        if(input$DownSignatureFormat=='tiff'){
            tiff(file, height=input$CHeight, width=input$CWidth)
            print(signature_plotter())
            dev.off()
        }
    }
)

# Pseudotime
#################################################################

observe({
    updateSelectInput(session, "MonocleRoot", choices=c(unique(levels(data2$data[[input$MSeuratObject]]$seurat_clusters))))
})

monocle_partition <- reactive({
    withProgress(message='Calculating Monocle Pseudotime',
        detail='Please stand by...',
        {
            shiny::setProgress(value=0.4, detail='Calculating...')

            
            tryCatch({
                get_earliest_principal_node <- function(cds, origin=input$MonocleRoot){
                    cell_ids <- which(colData(cds)[, "seurat_clusters"] == origin)
                    
                    closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
                    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
                    root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
                    
                    root_pr_nodes
                }

                mon <- as.cell_data_set(data2$data[[input$SSeuratObject]])
                mon <- cluster_cells(mon)
                mon <- learn_graph(mon)
                mon <- order_cells(mon, root_pr_nodes=get_earliest_principal_node(mon))

            monplot <- plot_cells(
                cds = mon,
                color_cells_by = "pseudotime",
                show_trajectory_graph = TRUE
                )
            return(monplot)
                
            },error=function(e)
            {
                status=paste("Monocle Error: ", e$message)
                showNotification(id="errorNotify", status, type='error', duration=NULL)
            })               
        }
    )
})

observe({
    output$monocleplot <- renderPlot({
        monocle_partition()
    }, height=input$CHeight, width=input$CWidth)
})

output$DownloadMonocle <- downloadHandler(
    filename=function(){
        paste('Pseudotime',input$DownMonocleFormat,sep='.')
    },
    content=function(file){   
        if(input$DownMonocleFormat=='jpeg'){
            jpeg(file, height=input$CHeight, width=input$CWidth)
            print(monocle_partition())
            dev.off()
        }
        if(input$DownMonocleFormat=='png'){
            png(file, height=input$CHeight, width=input$CWidth)
            print(monocle_partition())
            dev.off()
        }
        if(input$DownMonocleFormat=='tiff'){
            tiff(file, height=input$CHeight, width=input$CWidth)
            print(monocle_partition())
            dev.off()
        }
    }
)