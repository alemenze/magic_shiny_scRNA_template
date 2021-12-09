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

#########################Enrichments################################
####################################################################

organism='org.Mm.eg.db'
korg='mmu'
msig_org='Mus musculus'

# ORA
############################################################################
ORAReactive <- reactive({
    if(input$CCOptions=='FALSE'){
        validate(need(length(input$CCCluster1)>0,message = "Please choose at least 1 comparison for each."),
            need(length(input$CCCluster2)>0,message = "Please choose at least 1 comparison for each.")
        )
    }
    if(input$CCOptions=='TRUE'){
       validate(need(length(input$CCCluster1)>0,message = "Please choose at least 1 comparison for each."),
            need(length(input$CCCluster2)>0,message = "Please choose at least 1 comparison for each."),
            need(input$CCSubIdentChoice != 'None', message='Please choose a subset group')
        )
    }

    withProgress(message='Running Over Representation',
        detail='Please stand by...',
        {
            shiny::setProgress(value = 0.15, detail = "Loading Reference Database")

            library(organism, character.only=TRUE)

            shiny::setProgress(value = 0.45, detail = "Pulling comparison data")

            DataSet=crosser()

            DataSet <- subset(DataSet, p_val_adj < 0.05)
            

            if(input$Enrichpval =='padj'){
                DataSet = subset(DataSet, p_val_adj < 0.05)
            }
            else if (input$Enrichpval=='pvalue'){
                DataSet = subset(DataSet, pvalue < 0.05)
            }

            DataSetList <- rownames(DataSet)
            DataSetList <- bitr(DataSetList, fromType='SYMBOL',toType='ENTREZID', OrgDb=organism)
            
            shiny::setProgress(value = 0.65, detail = "Executing Over Representation Analysis")
            if(input$ORASet=='GO'){
                ora <- enrichGO(gene=DataSetList$ENTREZID, 
                    ont =input$ORAOnt, 
                    keyType = "ENTREZID", 
                    pvalueCutoff = 0.05,
                    OrgDb = organism, 
                    pAdjustMethod = input$ORAPVal)

                return(ora)
            }
            else if(input$ORASet=='kegg'){
                ora <- enrichKEGG(gene=DataSetList$ENTREZID,
                    organism=korg,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = input$ORAPVal
                )

                return(ora)
            }
            else if(input$ORASet=='msigdb'){
                gene_sets=msigdbr(species=msig_org, category=input$ORAMSig)
                gene_sets$gs_name = gsub("_", " ", gene_sets$gs_name, fixed = TRUE)
                msigdbr_t2g = gene_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

                ora <- enricher(DataSetList$ENTREZID, 
                    TERM2GENE=msigdbr_t2g, 
                    pvalueCutoff = 0.05,
                    pAdjustMethod = input$ORAPVal
                )

                return(ora)
            }
        }
    )

})


output$ora_table <- renderDataTable({
    DataIn <- ORAReactive()
    DT::datatable(as.data.frame(DataIn), style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

### Plot part
observe({
    output$oraplot <- renderPlot({
        DataIn <- ORAReactive()

        if(input$ORASet=='kegg'){
            updateSelectizeInput(session, 'ORAKeggChoice', choices=DataIn@result$Description)
        }
        if(input$ORAPlotType=='dotplot'){
            plot <- dotplot(DataIn, showCategory=as.numeric(input$ORANCategories), font.size=as.numeric(input$ORAFontSize)) + scale_y_discrete(labels=function(x) str_wrap(x, width=as.numeric(input$ORAYWidth)))
            plot
        }
        else if(input$ORAPlotType=='barchart'){
            plot <- barplot(DataIn, showCategory=as.numeric(input$ORANCategories),font.size=as.numeric(input$ORAFontSize))
            plot
        }
        else if(input$ORAPlotType=='emap'){
            plot <- emapplot(DataIn, showCategory = as.numeric(input$ORANCategories), font.size=as.numeric(input$ORAFontSize))
            plot
        }
    }, height=input$EnrichHeight, width=input$EnrichWidth)
})

output$DownloadORATable <- downloadHandler(
    filename=function(){
        paste(input$EnrichComp,'_',input$ORASet,'_ORA','.csv',sep='')
    },
    content = function(file){
        write.csv(ORAReactive(), file)
    }
)

output$DownloadORAPlot <- downloadHandler(
    filename=function(){
        paste(input$EnrichComp,'_',input$ORASet,'_ORA.',input$DownORAFormat,sep='')
    },
    content=function(file){   
        DataIn <- ORAReactive()
        if(input$ORAPlotType=='dotplot'){
            plot <- dotplot(DataIn, showCategory=as.numeric(input$ORANCategories), font.size=as.numeric(input$ORAFontSize)) + scale_y_discrete(labels=function(x) str_wrap(x, width=as.numeric(input$ORAYWidth)))
        }
        else if(input$ORAPlotType=='barchart'){
            plot <- barplot(DataIn, showCategory=as.numeric(input$ORANCategories),font.size=as.numeric(input$ORAFontSize))
        }
        else if(input$ORAPlotType=='emap'){
            plot <- emapplot(DataIn, showCategory = as.numeric(input$ORANCategories), font.size=as.numeric(input$ORAFontSize))
        }
        if(input$DownORAFormat=='jpeg'){
            jpeg(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
        if(input$DownORAFormat=='png'){
            png(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
        if(input$DownORAFormat=='tiff'){
            tiff(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
    }
)
