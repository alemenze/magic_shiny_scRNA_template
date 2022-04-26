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
    GenesList=unique(rownames(data2$data[[input$MSeuratObject]]))
    updateSelectizeInput(session, "HSelectedGenes", choices=GenesList, server=TRUE, options = list(maxOptions = 50))  
})
observe({
    updateSelectInput(session, "HGroupBy", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters') 
})

heatmap_plotter <- reactive({
    if(input$HeatTop=='Tops'){
        dataM$markerslist[[input$MTSeuratObject]] %>% group_by(cluster) %>% top_n(n=input$HeatMapTops, wt=avg_log2FC) -> top

        mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

        heatmap_plot <- DoHeatmap(data2$data[[input$MSeuratObject]], features=top$gene, group.by=input$HGroupBy) + NoLegend()+scale_fill_gradientn(colours = rev(mapal)) +theme(axis.text.y = element_text(size = input$HLabSize))+theme(
            axis.text.x = element_text(size=as.numeric(input$MGAxisSize)),
            plot.title=element_text(size=as.numeric(input$MGTitleSize)),
            legend.key.size = unit(as.numeric(input$MGLegendKeySize), 'cm'),
            legend.title = element_text(size=as.numeric(input$MGLegendTitleSize)),
            legend.text = element_text(size=as.numeric(input$MGLegendFontSize))
        )

    }
    if(input$HeatTop=='Select'){
        validate(need(length(input$HSelectedGenes)>0, message = "Please choose at least 1 gene."))
        mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

        heatmap_plot <- DoHeatmap(data2$data[[input$MSeuratObject]], features=c(input$HSelectedGenes, group.by=input$HGroupBy)) + NoLegend()+scale_fill_gradientn(colours = rev(mapal)) +theme(axis.text.y = element_text(size = input$HLabSize)) +theme(
            axis.text.x = element_text(size=as.numeric(input$MGAxisSize)),
            plot.title=element_text(size=as.numeric(input$MGTitleSize)),
            legend.key.size = unit(as.numeric(input$MGLegendKeySize), 'cm'),
            legend.text = element_text(size=as.numeric(input$MGLegendFontSize))
        )
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
    GenesList=unique(rownames(data2$data[[input$MSeuratObject]]))
    updateSelectizeInput(session, "DSelectedGenes", choices=GenesList, server=TRUE, options = list(maxOptions = 50))  
})
observe({
    updateSelectInput(session, "DGroupBy", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')   
})

dotplot_plotter <- reactive({
    if(input$DotTop=='Tops'){
        dataM$markerslist[[input$MTSeuratObject]] %>% group_by(cluster) %>% top_n(n=input$DotTops, wt=avg_log2FC) -> tops

        dotplot_plot <- DotPlot(data2$data[[input$MSeuratObject]], features=c(tops$gene), group.by=input$DGroupBy)+theme(
            axis.text.x = element_text(size=as.numeric(input$MGAxisSize)),
            axis.text.y = element_text(size=as.numeric(input$MGAxisSize)),
            plot.title=element_text(size=as.numeric(input$MGTitleSize)),
            legend.key.size = unit(as.numeric(input$MGLegendKeySize), 'cm'),
            legend.text = element_text(size=as.numeric(input$MGLegendFontSize))
        )

    }
    if(input$DotTop=='Select'){
        validate(need(length(input$DSelectedGenes)>0, message = "Please choose at least 1 gene."))
        dotplot_plot <- DotPlot(data2$data[[input$MSeuratObject]], features=c(input$DSelectedGenes),group.by=input$DGroupBy)+theme(
            axis.text.x = element_text(size=as.numeric(input$MGAxisSize)),
            axis.text.y = element_text(size=as.numeric(input$MGAxisSize)),
            plot.title=element_text(size=as.numeric(input$MGTitleSize)),
            legend.key.size = unit(as.numeric(input$MGLegendKeySize), 'cm'),
            legend.text = element_text(size=as.numeric(input$MGLegendFontSize))
        )
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

# Avg Expression Tables
#################################################################
observe({
    updateSelectInput(session, "EGroupBy", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')   
})

exp_table_setup <- reactive({
    temp <- DotPlot(data2$data[[input$MSeuratObject]], features=unique(rownames(data2$data[[input$MSeuratObject]])), group.by=input$EGroupBy, assay='SCT')
    temp <- temp$data[,colnames(temp$data) %in% c('avg.exp','pct.exp','features.plot','id')]
    return(temp)
})

output$expressions <- renderDataTable({
    DT::datatable(exp_table_setup(), style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

output$DownloadAvgExp <- downloadHandler(
    filename=function(){
        paste('avg_expression_table_',input$MSeuratObject,'.csv',sep='')
    },
    content = function(file){
        write.csv(exp_table_setup(), file)
    }
)

# Avg Expression Heatmap
#################################################################
observe({
    GenesList=unique(rownames(data2$data[[input$MSeuratObject]]))
    updateSelectizeInput(session, "HEgene_select", choices=GenesList, server=TRUE, options = list(maxOptions = 50))   
})

exp_hm <- reactive({
    cluster_r=TRUE
    if (input$HERowClust=='FALSE'){
        cluster_r=FALSE
    }
    cluster_c=TRUE
    if (input$HEColClust=='FALSE'){
        cluster_c=FALSE
    }

    shownames=FALSE
    if (input$HEShowNames=='TRUE'){
        shownames=TRUE
    }

    my_colors = brewer.pal(n = 11, name = "RdBu")
    my_colors = colorRampPalette(my_colors)(50)
    my_colors = rev(my_colors)

    validate(need(length(input$HEgene_select)>1, message = "Please choose at least 2 genes."))
    chosen_genes <- c(input$HEgene_select)

    temp <- DotPlot(data2$data[[input$MSeuratObject]], features=unique(rownames(data2$data[[input$MSeuratObject]])), group.by=input$EGroupBy, assay='SCT')
    temp <- temp$data %>% select('features.plot','avg.exp','id') %>% group_by(features.plot) %>% pivot_wider(names_from=id, values_from=avg.exp) %>% drop_na(features.plot) %>% unnest(cols = everything()) 
    temp <- as.data.frame(temp)
    temp[is.na(temp)] = 0
    temp <- temp %>% tibble::column_to_rownames(var='features.plot')
    temp <- temp[rowSums(temp[])>0,]
    
    DataSet <- temp[row.names(temp) %in% chosen_genes,]

    plot <- pheatmap(DataSet, cluster_rows=cluster_r, cluster_cols=cluster_c, color=my_colors,
        scale=input$HEScale, show_rownames=shownames, fontsize_col=input$HEXsize, angle_col=input$HEang)
    
    return(plot)

})

observe({
    output$heatmap_avg_exp <- renderPlot({
        exp_hm()
    },height=input$HHeight, width=input$HWidth)
})

output$DownloadHE <- downloadHandler(
    filename=function(){
        paste('avg_expression_heatmap',input$DownHEFormat,sep='.')
    },
    content=function(file){   
        if(input$DownHEFormat=='jpeg'){
            jpeg(file, height=input$HHeight, width=input$HWidth)
            print(exp_hm())
            dev.off()
        }
        if(input$DownHEFormat=='png'){
            png(file, height=input$HHeight, width=input$HWidth)
            print(exp_hm())
            dev.off()
        }
        if(input$DownHEFormat=='tiff'){
            tiff(file, height=input$HHeight, width=input$HWidth)
            print(exp_hm())
            dev.off()
        }
    }
)