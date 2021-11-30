observe({
    updateSelectInput(session, "CCSeuratObject", choices=c(names(data2$data)))
})
observe({
    updateSelectInput(session, "CCIdent", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
})
observe({
    updateSelectizeInput(session, "CCCluster1", choices=c(unique(levels(data2$data[[input$CCSeuratObject]][[input$CCIdent]][[input$CCIdent]]))), server=TRUE, options = list(maxOptions = 50))   
})
observe({
    updateSelectizeInput(session, "CCCluster2", choices=c(unique(levels(data2$data[[input$CCSeuratObject]][[input$CCIdent]][[input$CCIdent]]))), server=TRUE, options = list(maxOptions = 50)) 
})
observe({
    updateSelectInput(session, "CCSubIdent", choices=c(colnames(data2$data[[input$SeuratObject]][[]]),'seurat_clusters'), selected='seurat_clusters')
})
observe({
    updateSelectInput(session, "CCSubIdentChoice", choices=c(unique(levels(data2$data[[input$CCSeuratObject]][[input$CCSubIdent]][[input$CCSubIdent]])),'None'), selected='None')
})
# Comparisons
####################################################################

crosser <- reactive({
    if(input$CCOptions=='FALSE'){
        validate(need(length(input$CCCluster1)>0,message = "Please choose at least 1 comparison for each."),
            need(length(input$CCCluster2)>0,message = "Please choose at least 1 comparison for each.")
        )
            
        withProgress(message='Performing cross comparison',
            detail='Please stand by...',
            {
                shiny::setProgress(value=0.6, detail='Performing cross comparison')
                comparison <- FindMarkers(data2$data[[input$CCSeuratObject]], 
                    ident.1=c(input$CCCluster1), 
                    ident.2=c(input$CCCluster2),
                    group.by=input$CCIdent,
                    logfc.threshold = input$CCLog,
                    min.pct=input$CCPct,
                    only.pos=input$CCPos
                )
            }
        )
    }
    if(input$CCOptions=='TRUE'){
       validate(need(length(input$CCCluster1)>0,message = "Please choose at least 1 comparison for each."),
            need(length(input$CCCluster2)>0,message = "Please choose at least 1 comparison for each."),
            need(input$CCSubIdentChoice != 'None', message='Please choose a subset group')
        )
            
        withProgress(message='Performing cross comparison',
            detail='Please stand by...',
            {
                shiny::setProgress(value=0.6, detail='Performing cross comparison')
                Idents(data2$data[[input$CCSeuratObject]]) <- input$CCSubIdent
                comparison <- FindMarkers(data2$data[[input$CCSeuratObject]], 
                    ident.1=c(input$CCCluster1), 
                    ident.2=c(input$CCCluster2),
                    group.by=input$CCIdent,
                    logfc.threshold = input$CCLog,
                    min.pct=input$CCPct,
                    only.pos=input$CCPos,
                    subset.ident=input$CCSubIdentChoice
                )
            }
        ) 
    }
    return(comparison)
    
})

output$crosscomparison <- renderDataTable({
    DT::datatable(crosser(), style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

output$DownloadCross <- downloadHandler(
    filename=function(){
        paste('cross_comparison_',input$CCCluster1,'_vs_',input$CCCluster2,'.csv',sep='')
    },
    content = function(file){
        write.csv(crosser(), file)
    }
)