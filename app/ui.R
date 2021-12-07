library(shiny)
require(shinyjs)
library(shinythemes)
require(shinycssloaders)

library(tidyr)
library(tidyverse)
library(dplyr)
library(DT)
library(rhandsontable)
library(colourpicker)
library(RColorBrewer)
library(ggplot2)
library(circlize)
library(pheatmap)

library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(msigdbr)
library(stringr)

library(Seurat)
library(hdf5r)
library(cowplot)
library(SingleR)
#library(scRNAseq)
#library(slingshot)
library(SeuratWrappers)
library(monocle3)
library(patchwork)
library(magrittr)

tagList(
    tags$head(
        #includeHTML(("www/GA.html")), If including google analytics
        tags$style(type = 'text/css','.navbar-brand{display:none;}')
    ),
    fluidPage(theme = shinytheme('yeti'),
            windowTitle = "scRNA",
            titlePanel(
                fluidRow(
                column(2, tags$a(href='http://www.bioinformagic.io/', tags$img(height =75 , src = "MaGIC_Icon_0f344c.svg")), align = 'center'), 
                column(10, fluidRow(
                  column(10, h1(strong('scRNA Interactive Visualization Tool'), align = 'center')),
                  column(10, h2(strong('Project'), align = 'center'))))
                ),
                windowTitle = "scRNA" ),
                tags$style(type='text/css', '.navbar{font-size:20px;}'),
                tags$style(type='text/css', '.nav-tabs{padding-bottom:20px;}'),
                tags$head(tags$style(".modal-dialog{ width:1300px}")),

        navbarPage(title ="", id='NAVTABS',

        ## Intro Page
##########################################################################################################################################################
            tabPanel('Introduction',
               fluidRow(
                  column(2,
                  ),
                  column(8,
                      column(12, align = "center", 
                          style="margin-bottom:25px;",
                          h2("Introduction")),
                      hr(),
                      column(12, align="center",
                          style="margin-bottom:50px;", 
                          column(2,),
                          column(8,markdown("Welcome to the scRNAseq data explorer by [the Molecular and Genomics Informatics Core](http://www.bioinformagic.io).
                              This tool is intended to take an initially processed scRNA data object from the Core's pipeline and allow you to visualize and dive deeper in to the data. 
                              Each launch is customizable, if you wish to had additional modules added to your specific project please contact the Core. 
                          ")),
                          column(2,)       
                      )
                  ),
                  column(2,
                  )
                ),

                fluidRow(
                  column(2,
                  ),
                  column(8,
                      column(12, align = "center", 
                          style="margin-bottom:25px;",
                          h2("Primary processing")),
                      hr(),
                      column(12, align="center",
                          style="margin-bottom:50px;", 
                          column(2,),
                          column(8,markdown("Initial primary processing is performed via [10X's Cellranger tool.](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

                            This tool wraps the initial basecalling and demultiplexing with the reference genome alignment and UMI/cell barcode deconvolution. An example command:
                            ```
                            cellranger count --id=SampleName \
                              --transcriptome=/path/to/references/refdata-gex-GRCh38-2020-A \
                              --fastqs=/path/to/fastqs \
                              --sample=samplename_in_fastqdir
                            ``` 
                          ")),
                          column(2,)       
                      )
                  ),
                  column(2,
                  )
                ),

                fluidRow(
                  column(2,
                  ),
                  column(8,
                      column(12, align = "center", 
                          style="margin-bottom:25px;",
                          h2("Secondary processing")),
                      hr(),
                      column(12, align="center",
                          style="margin-bottom:50px;", 
                          column(2,),
                          column(8,markdown("Secondary processing is performed via [Seurat](https://satijalab.org/seurat/) in R. 

                            These steps include the initial filtering of low quality cells, regression of data, integration (if applicable), and UMAP clustering by nearest neighbors
                            
                            "),
                            actionButton('seurat_processing','Seurat Example',class='btn btn-info',style="margin-top:15px;")
                            ),
                          column(2,)       
                      )
                  ),

                  column(2,
                  )
                )
               
            ),

        ## UMAP Page
##########################################################################################################################################################
            tabPanel('UMAP Clustering',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('UMAP Plots', align='center'),
                            selectInput("SeuratObject", label="Select Object Source", choices=NULL),
                            selectInput("SeuratMeta", label="Select Grouping", choices=NULL),
                            radioButtons("UMAPLabel", label="Label Clusters", inline=TRUE, 
                                        choices=c("True"='TRUE',"False"='FALSE'), selected="TRUE"),
                            sliderInput("UMAPPointSize","Point Size: ", min=0.1, max=10, step=0.1, value=1),
                            sliderInput("UMAPLabelSize","Label Size: ", min=1, max=20, step=1, value=4),
                            sliderInput("UMAPAxisSize", "Axis Label Size: ", min=1, max=50, step=1, value=10),
                            sliderInput("UMAPTitleSize", "Title Size: ", min=1, max=30, step=1, value=10),
                            sliderInput("UMAPLegendKeySize", "Legend Key Size: ", min=0.25, max=5, step=0.25, value=1),
                            sliderInput("UMAPLegendFontSize", "Legend Font Size: ", min=1, max=30, step=1, value=10),
                            sliderInput('CHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                            sliderInput('CWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                        )
                    ),
                    column(9,
                        tabsetPanel(id='UMAP',
                            tabPanel(title='UMAP Plot', hr(),
                            withSpinner(type=6, color='#5bc0de',
                                        plotOutput("umapplot", height='100%')
                                ),
                            fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownUMAPFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadUMAP', 'Download the UMAP Plot'),style="margin-bottom:50px;")
                              )
                            )
                        )
                    )
                )
            ),

            ## Gene Plots Page
##########################################################################################################################################################
            tabPanel('Gene Plots',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('Gene Plots', align='center'),
                            selectInput("GSeuratObject", label="Select Object Source", choices=NULL),
                            conditionalPanel("input.GenePlots=='Feature Plot' || input.GenePlots=='Violin Plot' || input.GenePlots=='Ridge Plot'",
                              selectizeInput("gene_select", label="Name of gene",
                                    choices = NULL, multiple=FALSE, options=list(placeholder='Search'))
                            ),
                            conditionalPanel("input.GenePlots=='Feature Plot'",
                              sliderInput("FeaturePointSize","Point Size: ", min=0.1, max=10, step=0.1, value=1),
                              column(6, colourInput("NullColor", "Null Color", "#B3B3B3",palette = "limited")),
                              column(6, colourInput("PosColor", "Positive Color", "#0000ff",palette = "limited"))
                                ),
                            conditionalPanel("input.GenePlots=='Violin Plot'",
                              sliderInput("VlnPointSize","Point Size: ", min=0.1, max=10, step=0.1, value=1),
                              selectInput("GroupBy", label="Group By", choices=NULL),
                              selectInput("SplitBy", label="Split By", choices=NULL)
                              #idents
                                ),
                            conditionalPanel("input.GenePlots=='Ridge Plot'",
                              selectInput("RGroupBy", label="Group By", choices=NULL)
                                ),
                            conditionalPanel("input.GenePlots=='QC Plot'",
                              selectInput("QCFeature", label="Feature to use:", choices=NULL),
                              sliderInput("QCVlnPointSize","Point Size: ", min=0.1, max=10, step=0.1, value=1),
                              selectInput("QCGroupBy", label="Group By", choices=NULL),
                              selectInput("QCSplitBy", label="Split By", choices=NULL)
                                ),
                            sliderInput("GPAxisSize", "Axis Label Size: ", min=1, max=50, step=1, value=10),
                            sliderInput("GPTitleSize", "Title Size: ", min=1, max=30, step=1, value=10),
                            sliderInput("GPLegendKeySize", "Legend Key Size: ", min=0.25, max=5, step=0.25, value=1),
                            sliderInput("GPLegendFontSize", "Legend Font Size: ", min=1, max=30, step=1, value=10),
                            sliderInput('GHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                            sliderInput('GWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                        )
                    ),
                    column(9,
                        tabsetPanel(id='GenePlots',
                            tabPanel(title='Feature Plot', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                        plotOutput("featureplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownFeatureFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadFeature', 'Download the Feature Plot'),style="margin-bottom:50px;")
                              )
                            ),

                            tabPanel(title='Violin Plot', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                        plotOutput("vlnplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownVlnFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadVln', 'Download the Vln Plot'),style="margin-bottom:50px;")
                              )
                            ),
                            tabPanel(title='Ridge Plot', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                        plotOutput("ridgeplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                  column(12, selectInput("DownRidgeFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                  column(12, downloadButton('DownloadRidge', 'Download the Ridge Plot'),style="margin-bottom:50px;")
                                )
                            ),
                            tabPanel(title='QC Plot', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                        plotOutput("vlnplot_qc", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                  column(12, selectInput("DownVlnQCFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                  column(12, downloadButton('DownloadVlnQC', 'Download the QC Plot'),style="margin-bottom:50px;")
                                )
                            )
                        )
                    )
                )
            ),

            ## Marker Genes Page
##########################################################################################################################################################
            tabPanel('Marker Genes',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('Marker Genes', align='center'),
                            selectInput("MSeuratObject", label="Select Object Source", choices=NULL),
                            selectInput("MTSeuratObject", label="Select List Source", choices=NULL),
                            conditionalPanel("input.MarkerGenes=='Tables'",
                              selectInput("ClusterMarkers", label="View Cluster: ", choices=NULL),
                              actionButton('marker_columns','Column Descriptors',class='btn btn-info',style="margin-top:15px;")
                            ),
                            conditionalPanel("input.MarkerGenes=='Heatmaps'",
                              radioButtons("HeatTop",label='Gene selection by:',inline=TRUE, choices=c('Top Genes per Cluster'='Tops','Chosen Genes'='Select'), selected='Tops'),
                              selectInput("HGroupBy", label="Group By", choices=NULL),
                              conditionalPanel("input.HeatTop=='Select'",
                                selectizeInput("HSelectedGenes", "Please list genes of choice", choices = NULL, multiple=TRUE, options=list(placeholder='Search'))
                              ),
                              conditionalPanel("input.HeatTop=='Tops'",
                                sliderInput("HeatMapTops","Number of markers per cluster (if available): ", min=5, max=50, step=5, value=10)
                              )
                            ),
                            conditionalPanel("input.MarkerGenes=='Dot Plots'",
                              radioButtons("DotTop",label='Gene selection by:',inline=TRUE, choices=c('Top Genes per Cluster'='Tops','Chosen Genes'='Select'), selected='Tops'),
                              selectInput("DGroupBy", label="Group By", choices=NULL),
                              conditionalPanel("input.DotTop=='Select'",
                                selectizeInput("DSelectedGenes", "Please list genes of choice", choices = NULL, multiple=TRUE, options=list(placeholder='Search'))
                              ),
                              conditionalPanel("input.DotTop=='Tops'",
                                sliderInput("DotTops","Number of markers per cluster (if available): ", min=5, max=50, step=5, value=10)
                              )
                            ),
                            conditionalPanel("input.MarkerGenes=='Heatmaps' || input.MarkerGenes=='Dot Plots'",
                              sliderInput("MGAxisSize", "Axis Label Size: ", min=1, max=50, step=1, value=10),
                              sliderInput("MGTitleSize", "Title Size: ", min=1, max=30, step=1, value=10),
                              sliderInput("MGLegendKeySize", "Legend Key Size: ", min=0.25, max=5, step=0.25, value=1),
                              sliderInput("MGLegendFontSize", "Legend Font Size: ", min=1, max=30, step=1, value=10),
                              sliderInput('HHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                              sliderInput('HWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='MarkerGenes',
                            tabPanel(title='Tables', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('markers')
                                ),
                              fluidRow(
                                    column(12, align='center',downloadButton('DownloadMarkers', 'Download the Markers Table'))
                                )
                            ),

                            tabPanel(title='Heatmaps', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                  plotOutput("heatmapplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownHeatFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadHeat', 'Download the Heatmap'),style="margin-bottom:50px;")
                              )
                            ),

                            tabPanel(title='Dot Plots', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                    plotOutput("dotplotplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownDotFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadDot', 'Download the Dot Plot'),style="margin-bottom:50px;")
                              )
                            )
                        )
                    )
                )
            ),

            ## Signatures Page
##########################################################################################################################################################
            tabPanel('Gene Signatures',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('Gene Signatures', align='center'),
                            selectInput("SSeuratObject", label="Select Object Source", choices=NULL),
                            conditionalPanel("input.GeneSignatures=='Cell Cycling'",
                              h2('Cell Cycle Signatures', align='center'),
                              radioButtons("CyclePlotType",label='Plot By:',inline=TRUE, choices=c('Violin Plot'='Vln','Feature Plot'='Feature'), selected='Vln'),
                              conditionalPanel("input.CyclePlotType=='Vln'",
                                radioButtons("CycleScore",label='Score By:',inline=TRUE, choices=c('S Score'='S.Score','G2M Score'='G2M.Score'), selected='S.Score'),
                                sliderInput("CycleVlnPointSize","Point Size: ", min=0.1, max=10, step=0.1, value=1),
                                selectInput("CycleGroupBy", label="Group By", choices=NULL),
                                selectInput("CycleSplitBy", label="Split By", choices=NULL)  
                              ),
                              conditionalPanel("input.CyclePlotType=='Feature'",
                                sliderInput("FeatureCyclePointSize","Point Size: ", min=0.1, max=10, step=0.1, value=1)
                              )
                            ),
                            conditionalPanel("input.GeneSignatures=='Gene Lists'",
                              h2('Pathway Signatures', align='center'),
                              markdown("Warning- These may fail regularly. This is due to the dropout nature of single cell and there being insufficient genes to calculate module scoring."),
                              #Add upload option?
                              radioButtons("SignSource",label='MSigDB Collection:',inline=TRUE, choices=c('Hallmark Gene Sets'='H',
                                'Positional Gene Sets'='C1',
                                'Curated Gene Sets'='C2',
                                'Regulatory Target Gene Sets'='C3',
                                'Computational Gene Sets'='C4',
                                'Ontology Gene Sets'='C5',
                                'Oncogenic Signature Gene Sets'='C6',
                                'Immunologic Signature Gene Sets'='C7',
                                'Cell Type Signature Gene Sets'='C8'), selected='H'),
                              selectInput("SignPath", label="Select Gene Set", choices=NULL),
                              sliderInput("SignVlnPointSize","Point Size: ", min=0.1, max=10, step=0.1, value=1),
                              selectInput("SignGroupBy", label="Group By", choices=NULL),
                              selectInput("SignSplitBy", label="Split By", choices=NULL)
                            ),
                            ## Add GSEA?
                            conditionalPanel("input.GeneSignatures=='Pseudotime'",
                              h2('Monocle Pseudotime',align='center'),
                              selectInput("MonocleRoot", label='Select Root Node', choices=NULL),
                              sliderInput("GSLegendTitleSize", "Legend Title Size: ", min=1, max=30, step=1, value=10)

                            ),
                            sliderInput("GSAxisSize", "Axis Label Size: ", min=1, max=50, step=1, value=10),
                            sliderInput("GSTitleSize", "Title Size: ", min=1, max=30, step=1, value=10),
                            sliderInput("GSLegendKeySize", "Legend Key Size: ", min=0.25, max=5, step=0.25, value=1),
                            sliderInput("GSLegendFontSize", "Legend Font Size: ", min=1, max=30, step=1, value=10),
                            sliderInput('CHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                            sliderInput('CWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                        )
                    ),
                    column(9,
                        tabsetPanel(id='GeneSignatures',
                            tabPanel(title='Cell Cycling', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                        plotOutput("cycleplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownCycleFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadCycle', 'Download the Cell Cycle Plot'),style="margin-bottom:50px;")
                              )
                            ),

                            tabPanel(title='Gene Lists', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                        plotOutput("signatureplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownSignatureFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadSignature', 'Download the Signature Plot'),style="margin-bottom:50px;")
                              )
                            
                            ),
                            tabPanel(title='Pseudotime', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                        plotOutput("monocleplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownMonocleFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadMonocle', 'Download the Pseudotime Plot'),style="margin-bottom:50px;")
                              )
                            
                            )
                        )
                    )
                )
            ),

            ## Cross Comparisons Page
##########################################################################################################################################################
            tabPanel('Cross Comparisons',
                fluidRow(
                    column(3,
                        wellPanel(
                            conditionalPanel("input.CrossComparisons=='Tables'",
                              h2('Cross Comparisons',align='center'),
                              selectInput("CCSeuratObject", label="Select Object Source", choices=NULL),
                              selectInput("CCIdent", label="Select Grouping", choices=NULL),
                              selectizeInput("CCCluster1", label="Comparison Numerator: ", multiple=TRUE, choices=NULL, options=list(placeholder='Select')),
                              selectizeInput("CCCluster2", label="Comparison Denominator: ", multiple=TRUE, choices=NULL, options=list(placeholder='Select')), 
                              radioButtons("CCPos",label='Only Positive Markers',inline=TRUE, choices=c('True'='TRUE','False'='FALSE'), selected='TRUE'),
                              sliderInput('CCLog', label='LogFC Threshold: ', min=0.05, max=1, step=0.05, value=0.25),
                              sliderInput('CCPct', label='Minimum Percentage: ', min=0.01, max=0.5, step=0.01, value=0.1),
                              radioButtons("CCOptions",label='Subset Grouping Options',inline=TRUE, choices=c('True'='TRUE','False'='FALSE'), selected='FALSE'),
                              conditionalPanel("input.CCOptions=='TRUE'",
                                selectInput("CCSubIdent", label="Select Subset Grouping", choices=NULL),
                                selectInput("CCSubIdentChoice", label="Subset Ident: ", choices=NULL)
                              )
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='CrossComparisons',
                            tabPanel(title='Tables', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('crosscomparison')
                                ),
                              fluidRow(
                                    column(12, align='center',downloadButton('DownloadCross', 'Download the Cross Comparison Table'))
                                )
                            
                            )
                        )
                    )
                )
            ),

            ## Idents Page
##########################################################################################################################################################
            tabPanel('Add Identities or Recluster',
                fluidRow(
                    column(3,
                        wellPanel(
                          conditionalPanel("input.IdentReclust=='Ident relabel'",
                              h2('Set new Identities', align='center'),
                              selectInput("ISeuratObject", label="Select Object Source", choices=NULL),
                              textInput('inputId', label='Name of the new Identities', value = "NewIdents", width = NULL, placeholder = 'NewIdents'),
                              actionButton('submitIdents',"Submit New Idents",class='btn btn-info')                            
                          ),
                          conditionalPanel("input.IdentReclust=='Recluster'",
                            h2('Reclustering', align='center'),
                            selectizeInput("RClusterSet", label="Select clusters: ", multiple=TRUE, choices=NULL, options=list(placeholder='Select')),
                            textInput('reclustId', label='Name of the new reclustering object:', value = "Reclust", width = NULL, placeholder = 'Reclust'),
                            radioButtons("AdvancedOptions",label='Advanced Options',inline=TRUE, choices=c('True'='TRUE','False'='FALSE'), selected='FALSE'),
                            conditionalPanel("input.AdvancedOptions=='TRUE'",
                              sliderInput('RNFeatures', label='nFeatures upper limit: ', min=500, max=10000, step=500, value=7000),
                              sliderInput('RNFeaturesLow', label='nFeatures lower limit: ', min=250, max=2000, step=50, value=1000),
                              sliderInput('RPercentMito', label='Mitochondrial percent cutoff: ', min=0, max=100, step=5, value=10),
                              sliderInput('RPercentRibo', label='Ribosomal percent cutoff: ', min=0, max=100, step=5, value=45),
                              sliderInput('RUmapDim', label='Umap Dimensions: ', min=2, max=30, step=1, value=15),
                              sliderInput('RResolution', label='Clustering Resolution: ', min=0.1, max=2, step=0.05, value=0.2)
                            ),
                            actionButton('submitReclust',"Submit Re-clustering",class='btn btn-info')  
                          )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='IdentReclust',
                          tabPanel(title='Ident relabel', hr(), 
                            withSpinner(type=6, color='#5bc0de',
                                      rHandsontableOutput('reident')
                            ),style="margin-bottom:25px;"),

                            tabPanel(title='Recluster', hr(), 
                              column(12,includeMarkdown("docs/recluster.md"), align='left', hr())
                            )
                        )
                    )
                )
            ),

            ## Cellphone DB
##########################################################################################################################################################
            tabPanel('CellphoneDB',
                fluidRow(
                    column(3,
                        wellPanel(
                            conditionalPanel("input.Cellphone=='Chord Diagrams'",
                              h2('CellphoneDB Chord Diagrams',align='center'),
                              selectInput("CellphoneChordSource", label="Select source dataset", choices=NULL),
                              p('Color will be built on the palette chosen here:'),
                              column(4, colourInput("chordcol1", "Color 1", "#FF0000",palette = "limited")),
                              column(4, colourInput("chordcol2", "Color 2", "#00FF00",palette = "limited")),
                              column(4, colourInput("chordcol3", "Color 3", "#0000FF",palette = "limited")),
                              sliderInput('ChordFontSize', label='Font size: ', min=1, max=20, step=1, value=2)
                            ),
                            conditionalPanel("input.Cellphone=='Heatmaps'",
                              h2('CellphoneDB Heatmaps',align='center'),
                              selectInput("CellphoneHMSource", label="Select source dataset", choices=NULL),
                              radioButtons('CPHRowClust', label='Cluster rows: ', inline=TRUE,
                                choices=c('True'='TRUE','False'='FALSE'), selected='TRUE'),
                              radioButtons('CPHColClust', label='Cluster columns: ', inline=TRUE,
                                choices=c('True'='TRUE','False'='FALSE'), selected='TRUE'),
                              radioButtons('CPHScale', label='Scale By: ', inline=TRUE,
                                choices=c('Row'='row','Column'='column','None'='none'), selected='none'),
                              sliderInput("CPHXsize", 'Label Size', min=1, max=30, step=1, value=10),
                              radioButtons("CPHang", label='X-Axis Angle', inline=TRUE,
                                choices=c('0'=0,'45'=45,'90'=90,'270'=270,'315'=315), selected=45),
                              radioButtons("CPHcolors", label='Heatmap Colors', inline=TRUE, 
                                choices=c('Default'='Default', 'Select'='Select'), selected='Default'),
                              conditionalPanel("input.CPHcolors=='Select'",
                                p('Color will be built on the palette chosen here:'),
                                column(4, colourInput("cphmcol1", "Color 1", "#FF0000",palette = "limited")),
                                column(4, colourInput("cphmcol2", "Color 2", "#00FF00",palette = "limited")),
                                column(4, colourInput("cphmcol3", "Color 3", "#0000FF",palette = "limited"))
                              )
                            ),
                            conditionalPanel("input.Cellphone=='Dot Plots'",
                              h2('CellphoneDB Dot Plots',align='center'),
                              selectInput("CellphoneDotSource", label="Select source dataset", choices=NULL),
                              selectizeInput("cpgenes_select", label="Select gene interaction pairs: ", multiple=TRUE, choices=NULL, options=list(placeholder='Select')),
                              selectizeInput("cpclusters_select", label="Select cluster interaction pairs: ", multiple=TRUE, choices=NULL, options=list(placeholder='Select')),
                              p('Color will be built on the palette chosen here:'),       
                              column(3, colourInput("dotcol1", "Color 1", "#000000",palette = "limited")),
                              column(3, colourInput("dotcol2", "Color 2", "#0000FF",palette = "limited")),
                              column(3, colourInput("dotcol3", "Color 3", "#FFFF00",palette = "limited")),
                              column(3, colourInput("dotcol4", "Color 4", "#FF0000",palette = "limited")),
                              sliderInput('cpdotX', label='X Axis Font size: ', min=1, max=30, step=1, value=5),
                              sliderInput('cpdotY', label='Y Axis Font size: ', min=1, max=30, step=1, value=5),
                              radioButtons("cpdotang", label='X-Axis Angle', inline=TRUE,
                                        choices=c('0'=0,'45'=45,'90'=90,'270'=270,'315'=315), selected=45)
                            ),                         

                            conditionalPanel("input.Cellphone=='Chord Diagrams' || input.Cellphone=='Dot Plots' || input.Cellphone=='Heatmaps'",
                              sliderInput('CPHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                              sliderInput('CPWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                            ),
                            conditionalPanel("input.Cellphone=='Cellphone Tables'",
                              h2('CellphoneDB Table Data',align='center'),
                              selectInput("CellphoneTableSource", label="Select source dataset", choices=NULL)
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='Cellphone',
                            tabPanel(title='Chord Diagrams', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                    plotOutput("chordplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownCPChordFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadCPChord', 'Download the Chord Diagram'),style="margin-bottom:50px;")
                              )                            
                            ),
                            tabPanel(title='Heatmaps', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                    plotOutput("cpheatmap", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownCPHMFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadCPHM', 'Download the CellphoneDB Heatmap'),style="margin-bottom:50px;")
                              )                            
                            ),
                            tabPanel(title='Dot Plots', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                    plotOutput("cpdotplot", height='100%')
                                ),
                              fluidRow(align='center',style="margin-top:25px;",
                                column(12, selectInput("DownCPDotFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                column(12, downloadButton('DownloadCPDot', 'Download the CellphoneDB Dot Plot'),style="margin-bottom:50px;")
                              )                            
                            ),
                            tabPanel(title='Cellphone Tables', hr(),
                              withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('cptable')
                                ),
                              fluidRow(
                                    column(12, align='center',downloadButton('DownloadCPTable', 'Download the CellphoneDB Table'))
                                )
                            )
                        )
                    )
                )
            ),
            
        ## Footer
##########################################################################################################################################################
            tags$footer(
                wellPanel(
                    fluidRow(
                        column(4, align='center',
                        tags$a(href="https://github.com/alemenze/magic_shiny_scRNA_template", icon("github", "fa-3x")),
                        tags$h4('GitHub to submit issues/requests')
                        ),
                        column(4, align='center',
                        tags$a(href="http://www.bioinformagic.io/", icon("magic", "fa-3x")),
                        tags$h4('Bioinfor-MaGIC Home Page')
                        ),
                        column(4, align='center',
                        tags$a(href="https://alemenze.github.io/", icon("address-card", "fa-3x")),
                        tags$h4("Developer's Page")
                        )
                    ),
                    fluidRow(
                        column(12, align='center',
                            HTML('<a href="https://www.youtube.com/watch?v=dQw4w9WgXcQ">
                            <p>&copy; 
                                <script language="javascript" type="text/javascript">
                                var today = new Date()
                                var year = today.getFullYear()
                                document.write(year)
                                </script>
                            </p>
                            </a>
                            ')
                        )
                    ) 
                )
            )
        )#Ends navbarPage,
    )#Ends fluidpage
)#Ends tagList
