# Clustering UMAP
#################################################################

observe({
    updateSelectInput(session, "SeuratObject", choices=c(names(data2$data)))
})
observe({
    updateSelectInput(session, "SeuratMeta", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
})

umapplotter <- reactive({
    umapplot <- DimPlot(data2$data[[input$SeuratObject]], reduction='umap', label=input$UMAPLabel, group.by=input$SeuratMeta,
        pt.size=input$UMAPPointSize, label.size=input$UMAPLabelSize) +theme(
            axis.text.x = element_text(size=as.numeric(input$UMAPAxisSize)),
            axis.text.y = element_text(size=as.numeric(input$UMAPAxisSize)),
            plot.title=element_text(size=as.numeric(input$UMAPTitleSize)),
            legend.key.size = unit(as.numeric(input$UMAPLegendKeySize), 'cm'),
            legend.text = element_text(size=as.numeric(input$UMAPLegendFontSize))
        )
    
    return(umapplot)
})


observe({
    output$umapplot <- renderPlot({
        umapplotter()
    }, height=input$CHeight, width=input$CWidth)
})

output$DownloadUMAP <- downloadHandler(
    filename=function(){
        paste('UMAP',input$DownUMAPFormat,sep='.')
    },
    content=function(file){   
        if(input$DownUMAPFormat=='jpeg'){
            jpeg(file, height=input$CHeight, width=input$CWidth)
            print(umapplotter())
            dev.off()
        }
        if(input$DownUMAPFormat=='png'){
            png(file, height=input$CHeight, width=input$CWidth)
            print(umapplotter())
            dev.off()
        }
        if(input$DownUMAPFormat=='tiff'){
            tiff(file, height=input$CHeight, width=input$CWidth)
            print(umapplotter())
            dev.off()
        }
    }
)

# Counts
#################################################################
observe({
    updateSelectInput(session, "SeuratSplit", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
})
counts_setup <- reactive({
    Idents(data2$data[[input$SeuratObject]]) <- as.character(input$SeuratMeta)
    ctable <- table(Idents(data2$data[[input$SeuratObject]]),data2$data[[input$SeuratObject]][[]][[input$SeuratSplit]])
    return(as.data.frame.matrix(ctable))
})

output$umapcounts <- renderDataTable({
    DT::datatable(counts_setup(), style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

output$DownloadCountTable <- downloadHandler(
    filename=function(){
        paste('cellcount_table_',input$SeuratObject,'.csv',sep='')
    },
    content = function(file){
        write.csv(counts_setup(), file)
    }
)