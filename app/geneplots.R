observe({
    GenesList=unique(rownames(data2$data[[input$GSeuratObject]]))
    updateSelectizeInput(session, "gene_select", choices=GenesList, server=TRUE, options = list(maxOptions = 50))   
})
observe({
    updateSelectInput(session, "GSeuratObject", choices=c(names(data2$data)))
})

# Feature Plots
#################################################################
feature_plotter <- reactive({
    validate(need(length(input$gene_select)>0, message = "Please choose at least 1 gene."))

    feature_plot <- FeaturePlot(data2$data[[input$GSeuratObject]], features=c(input$gene_select), 
        cols=c(input$NullColor, input$PosColor),
        pt.size=input$FeaturePointSize) +theme(
            axis.text.x = element_text(size=as.numeric(input$GPAxisSize)),
            axis.text.y = element_text(size=as.numeric(input$GPAxisSize)),
            plot.title=element_text(size=as.numeric(input$GPTitleSize)),
            legend.key.size = unit(as.numeric(input$GPLegendKeySize), 'cm'),
            legend.text = element_text(size=as.numeric(input$GPLegendFontSize))
        )
    
    return(feature_plot)
})

observe({
    output$featureplot <- renderPlot({
        feature_plotter()
    }, height=input$GHeight, width=input$GWidth)
})

output$DownloadFeature <- downloadHandler(
    filename=function(){
        paste('Feature',input$DownFeatureFormat,sep='.')
    },
    content=function(file){   
        if(input$DownFeatureFormat=='jpeg'){
            jpeg(file, height=input$GHeight, width=input$GWidth)
            print(feature_plotter())
            dev.off()
        }
        if(input$DownFeatureFormat=='png'){
            png(file, height=input$GHeight, width=input$GWidth)
            print(feature_plotter())
            dev.off()
        }
        if(input$DownFeatureFormat=='tiff'){
            tiff(file, height=input$GHeight, width=input$GWidth)
            print(feature_plotter())
            dev.off()
        }
    }
)

# Vln Plots
#################################################################
observe({
    updateSelectInput(session, "GroupBy", choices=c(colnames(data2$data[[input$GSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
    updateSelectInput(session, "SplitBy", choices=c(colnames(data2$data[[input$GSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')   
})

vln_plotter <- reactive({
    validate(need(length(input$gene_select)>0, message = "Please choose at least 1 gene."))

    vln_plot <- VlnPlot(data2$data[[input$GSeuratObject]], features=c(input$gene_select), 
        pt.size=input$VlnPointSize,
        group.by=input$GroupBy,
        split.by=input$SplitBy
        ) +theme(
            axis.text.x = element_text(size=as.numeric(input$GPAxisSize)),
            axis.text.y = element_text(size=as.numeric(input$GPAxisSize)),
            plot.title=element_text(size=as.numeric(input$GPTitleSize)),
            legend.key.size = unit(as.numeric(input$GPLegendKeySize), 'cm'),
            legend.text = element_text(size=as.numeric(input$GPLegendFontSize))
        )
    
    return(vln_plot)
})

observe({
    output$vlnplot <- renderPlot({
        vln_plotter()
    }, height=input$GHeight, width=input$GWidth)
})

output$DownloadVln <- downloadHandler(
    filename=function(){
        paste('Vln',input$DownVlnFormat,sep='.')
    },
    content=function(file){   
        if(input$DownVlnFormat=='jpeg'){
            jpeg(file, height=input$GHeight, width=input$GWidth)
            print(vln_plotter())
            dev.off()
        }
        if(input$DownVlnFormat=='png'){
            png(file, height=input$GHeight, width=input$GWidth)
            print(vln_plotter())
            dev.off()
        }
        if(input$DownVlnFormat=='tiff'){
            tiff(file, height=input$GHeight, width=input$GWidth)
            print(vln_plotter())
            dev.off()
        }
    }
)

# Ridge Plots
#################################################################
observe({
    updateSelectInput(session, "RGroupBy", choices=c(colnames(data2$data[[input$GSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
})

ridge_plotter <- reactive({
    validate(need(length(input$gene_select)>0, message = "Please choose at least 1 gene."))

    ridge_plot <- RidgePlot(data2$data[[input$GSeuratObject]], features=c(input$gene_select), 
        group.by=input$RGroupBy
        ) +theme(
            axis.text.x = element_text(size=as.numeric(input$GPAxisSize)),
            axis.text.y = element_text(size=as.numeric(input$GPAxisSize)),
            plot.title=element_text(size=as.numeric(input$GPTitleSize)),
            legend.key.size = unit(as.numeric(input$GPLegendKeySize), 'cm'),
            legend.text = element_text(size=as.numeric(input$GPLegendFontSize))
        )
    
    return(ridge_plot)
})

observe({
    output$ridgeplot <- renderPlot({
        ridge_plotter()
    }, height=input$GHeight, width=input$GWidth)
})

output$DownloadRidge <- downloadHandler(
    filename=function(){
        paste('Ridge',input$DownRidgeFormat,sep='.')
    },
    content=function(file){   
        if(input$DownRidgeFormat=='jpeg'){
            jpeg(file, height=input$GHeight, width=input$GWidth)
            print(ridge_plotter())
            dev.off()
        }
        if(input$DownRidgeFormat=='png'){
            png(file, height=input$GHeight, width=input$GWidth)
            print(ridge_plotter())
            dev.off()
        }
        if(input$DownRidgeFormat=='tiff'){
            tiff(file, height=input$GHeight, width=input$GWidth)
            print(ridge_plotter())
            dev.off()
        }
    }
)


# QC Plots
#################################################################
observe({
    updateSelectInput(session, "QCFeature", choices=c(colnames(data2$data[[input$GSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
    updateSelectInput(session, "QCGroupBy", choices=c(colnames(data2$data[[input$GSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
    updateSelectInput(session, "QCSplitBy", choices=c(colnames(data2$data[[input$GSeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')   
})

vln_plotter_qc <- reactive({
    validate(need(length(input$GroupBy)>0, message = "Please choose at least 1 feature."))

    vln_plot <- VlnPlot(data2$data[[input$GSeuratObject]], features=c(input$QCFeature), 
        pt.size=input$QCVlnPointSize,
        group.by=input$QCGroupBy,
        split.by=input$QCSplitBy
        ) +theme(
            axis.text.x = element_text(size=as.numeric(input$GPAxisSize)),
            axis.text.y = element_text(size=as.numeric(input$GPAxisSize)),
            plot.title=element_text(size=as.numeric(input$GPTitleSize)),
            legend.key.size = unit(as.numeric(input$GPLegendKeySize), 'cm'),
            legend.text = element_text(size=as.numeric(input$GPLegendFontSize))
        )
    
    return(vln_plot)
})

observe({
    output$vlnplot_qc <- renderPlot({
        vln_plotter_qc()
    }, height=input$GHeight, width=input$GWidth)
})

output$DownloadVlnQC <- downloadHandler(
    filename=function(){
        paste('QCVln',input$DownVlnQCFormat,sep='.')
    },
    content=function(file){   
        if(input$DownVlnQCFormat=='jpeg'){
            jpeg(file, height=input$GHeight, width=input$GWidth)
            print(vln_plotter_qc())
            dev.off()
        }
        if(input$DownVlnQCFormat=='png'){
            png(file, height=input$GHeight, width=input$GWidth)
            print(vln_plotter_qc())
            dev.off()
        }
        if(input$DownVlnQCFormat=='tiff'){
            tiff(file, height=input$GHeight, width=input$GWidth)
            print(vln_plotter_qc())
            dev.off()
        }
    }
)