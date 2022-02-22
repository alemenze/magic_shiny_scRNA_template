# Chord Diagrams
#################################################################

observe({
    updateSelectInput(session, "CellphoneChordSource", choices=c(names(cellphone_chords)))
})

chord_plotter <- reactive({
    col_fun <- colorRampPalette(c(input$chordcol1,input$chordcol2,input$chordcol3))(length(unique(cellphone_chords[[input$CellphoneChordSource]]$SOURCE)))
    par(cex=input$ChordFontSize, mar=c(0,0,0,0))
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
    validate(need(input$CellphoneDotSource, message = "Please choose data source."))

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
    interactions_clusters <- unique(plot_start$clusters)

    updateSelectizeInput(session, "cpgenes_select", choices=c(interactions_genes), server=TRUE)
    updateSelectizeInput(session, "cpclusters_select", choices=(interactions_clusters), server=TRUE) 
})

cpdot_plotter <- reactive({
    validate(need(length(input$cpgenes_select)>0, message = "Please choose at least 1 gene interaction pair."))
    validate(need(length(input$cpclusters_select)>0, message = "Please choose at least 1 cluster interaction pair."))
    
    selected_rows = c(input$cpgenes_select)
    selected_columns = c(input$cpclusters_select)

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

    my_palette <- colorRampPalette(c(input$dotcol1, input$dotcol2, input$dotcol3, input$dotcol4), alpha=TRUE)(n=399)

    plot <- ggplot(plot.data,aes(x=clusters,y=pair)) +
        geom_point(aes(size=-log10(pvalue),color=mean)) +
        scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.text=element_text(size=as.numeric(input$cpdotX), colour = "black"),
                axis.text.x = element_text(angle = as.numeric(input$cpdotang), hjust = 1),
                axis.text.y = element_text(size=as.numeric(input$cpdotY), colour = "black"),
                axis.title=element_blank(),
                panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

    return(plot)
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
# Cellphone table
#################################################################

observe({
    updateSelectInput(session, "CellphoneTableSource", choices=c(names(cellphone_dots)))
})

cptable_setup <- reactive({
    validate(need(input$CellphoneTableSource, message = "Please choose data source."))

    selected_rows = NULL
    selected_columns = NULL

    all_pval = cellphone_dots[[input$CellphoneTableSource]][[1]]
    all_means = cellphone_dots[[input$CellphoneTableSource]][[2]]

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

output$cptable <- renderDataTable({
    DT::datatable(cptable_setup(), style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

output$DownloadCPTable <- downloadHandler(
    filename=function(){
        paste('cellphonetable_',input$CellphoneTableSource,'.csv',sep='')
    },
    content = function(file){
        write.csv(cptable_setup(), file)
    }
)

# Heatmaps
#################################################################

observe({
    updateSelectInput(session, "CellphoneHMSource", choices=c(names(cellphone_chords)))
})

cp_heatmap_plotter <- reactive({
    cluster_r=TRUE
    if (input$CPHRowClust=='FALSE'){
        cluster_r=FALSE
    }
    cluster_c=TRUE
    if (input$CPHColClust=='FALSE'){
        cluster_c=FALSE
    }

    if(input$CPHcolors=='Default'){
        my_colors = brewer.pal(n = 11, name = "RdBu")
        my_colors = colorRampPalette(my_colors)(50)
        my_colors = rev(my_colors)
    } else {
        my_colors <- colorRampPalette(c(input$cphmcol1, input$cphmcol2, input$cphmcol3), alpha=TRUE)(n=399)
    }
    

    pos_pairs <- expand.grid(cellphone_chords[[input$CellphoneHMSource]]$SOURCE, cellphone_chords[[input$CellphoneHMSource]]$TARGET)
    pos_pairs$pairs <- paste(pos_pairs$Var1, pos_pairs$Var2, sep='|')
    pairs <- unique(pos_pairs$pairs)

    df <- data.frame(matrix(NA, nrow = length(unique(cellphone_chords[[input$CellphoneHMSource]]$SOURCE)), ncol = (length(unique(cellphone_chords[[input$CellphoneHMSource]]$TARGET)))+1))

    df <- data.frame(Clusters=c(unique(cellphone_chords[[input$CellphoneHMSource]]$SOURCE)))
    for(f in unique(cellphone_chords[[input$CellphoneHMSource]]$SOURCE)){
        counts <- c()
        
        for(i in pairs){
            var1 <- strsplit(i, '|')[[1]][1]
            var2 <- strsplit(i, '|')[[1]][3]
            if(var1==f){
                row <- cellphone_chords[[input$CellphoneHMSource]][which(cellphone_chords[[input$CellphoneHMSource]]$SOURCE == var1 & cellphone_chords[[input$CellphoneHMSource]]$TARGET == var2 ),]
                counts <- c(counts, row$count)
            }
            
        }
        df[[paste(f)]] <- counts
    }
    rownames(df) <- df$Clusters
    df <- df[order(df$Clusters),]
    df <- df[,-1]
    df <- df %>% dplyr::select(order(colnames(df)))

    p <- pheatmap(df,  cluster_rows=cluster_r, cluster_cols=cluster_c, color=my_colors,
            scale=input$CPHScale, fontsize_col=input$CPHXsize, fontsize_row=input$CPHXsize, angle_col=input$CPHang)

    return(p)
})

observe({
    output$cpheatmap <- renderPlot({
        cp_heatmap_plotter()
    }, height=input$CPHeight, width=input$CPWidth)
})

output$DownloadCPHM <- downloadHandler(
    filename=function(){
        paste('heatmap',input$DownCPHMFormat,sep='.')
    },
    content=function(file){   
        if(input$DownCPHMFormat=='jpeg'){
            jpeg(file, height=input$CPHeight, width=input$CPWidth)
            print(cp_heatmap_plotter)
            dev.off()
        }
        if(input$DownCPHMFormat=='png'){
            png(file, height=input$CPHeight, width=input$CPWidth)
            print(cp_heatmap_plotter)
            dev.off()
        }
        if(input$DownCPHMFormat=='tiff'){
            tiff(file, height=input$CPHeight, width=input$CPWidth)
            print(cp_heatmap_plotter)
            dev.off()
        }
    }
)