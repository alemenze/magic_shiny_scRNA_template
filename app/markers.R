observe({
    updateSelectInput(session, "MSeuratObject", choices=c(names(data2$data)))
})
observe({
    updateSelectInput(session, "MTSeuratObject", choices=c(names(dataM$markerslist)))
})

# Marker Tables
#################################################################

observe({
    updateSelectInput(session, "ClusterMarkers", choices=c(unique(levels(data2$data[[input$MSeuratObject]]$seurat_clusters)),'all'), selected='all')
})

output$markers <- renderDataTable({
    if(input$ClusterMarkers=='all'){
        df=dataM$markerslist[[input$MTSeuratObject]]
    } else {
        df= subset(dataM$markerslist[[input$MTSeuratObject]], cluster == input$ClusterMarkers)
    }

    DT::datatable(df, style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

output$DownloadMarkers <- downloadHandler(
    filename=function(){
        paste('topmarkers.csv',sep='')
    },
    content = function(file){
        if(input$ClusterMarkers=='all'){
            df=dataM$markerslist[[input$MTSeuratObject]]
        } else {
            df= subset(dataM$markerslist[[input$MTSeuratObject]], cluster == input$ClusterMarkers)
        }
        write.csv(df, file)
    }
)

# Heatmaps
#################################################################
observe({
    GenesList=unique(all.markers$gene)
    updateSelectizeInput(session, "HSelectedGenes", choices=GenesList, server=TRUE, options = list(maxOptions = 50)) 
    updateSelectInput(session, "HGroupBy", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')  
})

heatmap_plotter <- reactive({
    if(input$HeatTop=='Tops'){
        dataM$markerslist[[input$MTSeuratObject]] %>% group_by(cluster) %>% top_n(n=input$HeatMapTops, wt=avg_log2FC) -> top

        mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

        heatmap_plot <- DoHeatmap(data2$data[[input$MSeuratObject]], features=top$gene, group.by=input$HGroupBy) + NoLegend()+scale_fill_gradientn(colours = rev(mapal)) +theme(axis.text.y = element_text(size = 0))

    }
    if(input$HeatTop=='Select'){
        validate(need(length(input$HSelectedGenes)>0, message = "Please choose at least 1 gene."))
        mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

        heatmap_plot <- DoHeatmap(data2$data[[input$MSeuratObject]], features=c(input$HSelectedGenes, group.by=input$HGroupBy)) + NoLegend()+scale_fill_gradientn(colours = rev(mapal)) +theme(axis.text.y = element_text(size = 0))
    }

    return(heatmap_plot)
})

observe({
    output$heatmapplot <- renderPlot({
        heatmap_plotter()
    }, height=input$HHeight, width=input$HWidth)
})

output$DownloadHeat <- downloadHandler(
    filename=function(){
        paste('Heatmap',input$DownHeatFormat,sep='.')
    },
    content=function(file){   
        if(input$DownHeatFormat=='jpeg'){
            jpeg(file, height=input$HHeight, width=input$HWidth)
            print(heatmap_plotter())
            dev.off()
        }
        if(input$DownHeatFormat=='png'){
            png(file, height=input$HHeight, width=input$HWidth)
            print(heatmap_plotter())
            dev.off()
        }
        if(input$DownHeatFormat=='tiff'){
            tiff(file, height=input$HHeight, width=input$HWidth)
            print(heatmap_plotter())
            dev.off()
        }
    }
)

# Dot Plots
#################################################################
observe({
    GenesList=unique(all.markers$gene)
    updateSelectizeInput(session, "DSelectedGenes", choices=GenesList, server=TRUE, options = list(maxOptions = 50))  
    updateSelectInput(session, "DGroupBy", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')   
})

dotplot_plotter <- reactive({
    if(input$DotTop=='Tops'){
        dataM$markerslist[[input$MTSeuratObject]] %>% group_by(cluster) %>% top_n(n=input$DotTops, wt=avg_log2FC) -> tops

        dotplot_plot <- DotPlot(data2$data[[input$MSeuratObject]], features=c(tops$gene), group.by=input$DGroupBy)

    }
    if(input$DotTop=='Select'){
        validate(need(length(input$DSelectedGenes)>0, message = "Please choose at least 1 gene."))
        dotplot_plot <- DotPlot(data2$data[[input$MSeuratObject]], features=c(input$DSelectedGenes),group.by=input$DGroupBy)
    }

    return(dotplot_plot)
})

observe({
    output$dotplotplot <- renderPlot({
        dotplot_plotter()
    }, height=input$HHeight, width=input$HWidth)
})

output$DownloadDot <- downloadHandler(
    filename=function(){
        paste('Dotplot',input$DownDotFormat,sep='.')
    },
    content=function(file){   
        if(input$DownDotFormat=='jpeg'){
            jpeg(file, height=input$HHeight, width=input$HWidth)
            print(dotplot_plotter())
            dev.off()
        }
        if(input$DownDotFormat=='png'){
            png(file, height=input$HHeight, width=input$HWidth)
            print(dotplot_plotter())
            dev.off()
        }
        if(input$DownDotFormat=='tiff'){
            tiff(file, height=input$HHeight, width=input$HWidth)
            print(dotplot_plotter())
            dev.off()
        }
    }
)