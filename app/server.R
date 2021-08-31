function(input, output, session) {
    load('data.RData')
    data2 <- reactiveValues()
    dataM <- reactiveValues()
    data2$data <- data
    dataM$markerslist <- markerslist

    options(shiny.maxRequestSize=30*1024^2)
    source('ui.R',local=TRUE)
    source('modals.R',local=TRUE)
    source('clustering.R',local=TRUE)
    source('geneplots.R',local=TRUE)
    source('markers.R',local=TRUE)
    source('signatures.R',local=TRUE)
    source('cross.R', local=TRUE)
    source('idents.R',local=TRUE)
}