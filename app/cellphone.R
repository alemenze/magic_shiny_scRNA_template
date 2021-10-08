# Chord Diagrams
#################################################################

observe({
    updateSelectInput(session, "CellphoneChordSource", choices=c(names(cellphone_chords)))
})

chord_plotter <- reactive({
    col_fun <- colorRampPalette(c('red','green','blue'))(length(unique(cellphone_chords[[input$CellphoneChordSource]]$SOURCE)))
    chord <- chordDiagram(cellphone_chords[[input$CellphoneChordSource]], col=col_fun, grid.col=col_fun)
    
    return(chord)
})

observe({
    output$chordplot <- renderPlot({
        chord_plotter()
    }, height=input$CPHeight, width=input$CPWidth)
})

output$DownloadCPChord <- downloadHandler(
    filename=function(){
        paste('chord',input$DownCPChordFormat,sep='.')
    },
    content=function(file){   
        if(input$DownCPChordFormat=='jpeg'){
            jpeg(file, height=input$CPHeight, width=input$CPWidth)
            print(chord_plotter())
            dev.off()
        }
        if(input$DownCPChordFormat=='png'){
            png(file, height=input$CPHeight, width=input$CPWidth)
            print(chord_plotter())
            dev.off()
        }
        if(input$DownCPChordFormat=='tiff'){
            tiff(file, height=input$CPHeight, width=input$CPWidth)
            print(chord_plotter())
            dev.off()
        }
    }
)

# DotPlot Diagrams
#################################################################
observe({
    updateSelectInput(session, "CellphoneDotSource", choices=c(names(cellphone_dots)))
})

dot_setup <- reactive({
    selected_rows = NULL
    selected_columns = NULL

    all_pval = cellphone_dots[[input$CellphoneDotSource]][[1]]
    all_means = cellphone_dots[[input$CellphoneDotSource]][[2]]

    intr_pairs = all_pval$interacting_pair
    all_pval = all_pval[,-c(1:11)]
    all_means = all_means[,-c(1:11)]

    if(is.null(selected_rows)){
    selected_rows = intr_pairs
    }

    if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
    }

    sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
    sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

    df_names = expand.grid(selected_rows, selected_columns)
    pval = unlist(sel_pval)
    pval[pval==0] = 0.0009
    plot.data = cbind(df_names,pval)
    pr = unlist(as.data.frame(sel_means))
    pr[pr==0] = 1
    plot.data = cbind(plot.data,log2(pr))
    colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

    return(plot.data)
})

observe({
    plot_start <- dot_setup()
    interactions_genes <- unique(plot_start$pair)
    interactions_clustes <- unique(plot_start$clusters)

    updateSelectizeInput(session, "cpgenes_select", choices=interactions_genes, server=TRUE, options = list(maxOptions = 50))
    updateSelectizeInput(session, "cpclusters_select", choices=interactions_clusters, server=TRUE, options = list(maxOptions = 50)) 
})

cpdot_plotter <- reactive({
    validate(need(length(input$cpgenes_select)>0, message = "Please choose at least 1 gene interaction pair."))
    validate(need(length(input$cpclusters_select)>0, message = "Please choose at least 1 cluster interaction pair."))
    
    selected_rows = input$cpgenes_select
    selected_columns = input$cpclusters_select

    all_pval = cellphone_dots[[input$CellphoneDotSource]][[1]]
    all_means = cellphone_dots[[input$CellphoneDotSource]][[2]]

    intr_pairs = all_pval$interacting_pair
    all_pval = all_pval[,-c(1:11)]
    all_means = all_means[,-c(1:11)]

    if(is.null(selected_rows)){
    selected_rows = intr_pairs
    }

    if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
    }

    sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
    sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

    df_names = expand.grid(selected_rows, selected_columns)
    pval = unlist(sel_pval)
    pval[pval==0] = 0.0009
    plot.data = cbind(df_names,pval)
    pr = unlist(as.data.frame(sel_means))
    pr[pr==0] = 1
    plot.data = cbind(plot.data,log2(pr))

    colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

    my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

    plot <- ggplot(plot.data,aes(x=clusters,y=pair)) +
        geom_point(aes(size=-log10(pvalue),color=mean)) +
        scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.text=element_text(size=14, colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text.y = element_text(size=12, colour = "black"),
                axis.title=element_blank(),
                panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
})

observe({
    output$cpdotplot <- renderPlot({
        cpdot_plotter()
    }, height=input$CPHeight, width=input$CPWidth)
})

output$DownloadCPDot <- downloadHandler(
    filename=function(){
        paste('cellphone_dotplot',input$DownCPDotFormat,sep='.')
    },
    content=function(file){   
        if(input$DownCPDotFormat=='jpeg'){
            jpeg(file, height=input$CPHeight, width=input$CPWidth)
            print(cpdot_plotter())
            dev.off()
        }
        if(input$DownCPDotFormat=='png'){
            png(file, height=input$CPHeight, width=input$CPWidth)
            print(cpdot_plotter())
            dev.off()
        }
        if(input$DownCPDotFormat=='tiff'){
            tiff(file, height=input$CPHeight, width=input$CPWidth)
            print(cpdot_plotter())
            dev.off()
        }
    }
)

output$cptable <- renderDataTable({
    DT::datatable(dot_setup(), style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

output$DownloadCPTable <- downloadHandler(
    filename=function(){
        paste('cellphonetable.csv',sep='')
    },
    content = function(file){
        write.csv(dot_setup(), file)
    }
)