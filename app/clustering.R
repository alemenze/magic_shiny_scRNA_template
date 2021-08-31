observe({
    updateSelectInput(session, "SeuratObject", choices=c(names(data2$data)))
})
observe({
    updateSelectInput(session, "SeuratMeta", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
})

umapplotter <- reactive({
    umapplot <- DimPlot(data2$data[[input$SeuratObject]], reduction='umap', label=input$UMAPLabel, group.by=input$SeuratMeta,
        pt.size=input$UMAPPointSize, label.size=input$UMAPLabelSize); 
    
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