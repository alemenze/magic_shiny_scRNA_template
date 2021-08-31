# Seurat processing
########################################################################################################
observeEvent(input$seurat_processing, {
    showModal(modalDialog(
        column(12,includeMarkdown("docs/seurat.md"), align='center', hr()),
        easyClose = TRUE,
        footer = modalButton("Close"),
        size='l'
    ))    
})

# Column headers
########################################################################################################
observeEvent(input$marker_columns, {
    showModal(modalDialog(
        column(12,includeMarkdown("docs/columns.md"), align='center', hr()),
        easyClose = TRUE,
        footer = modalButton("Close"),
        size='l'
    ))    
})
