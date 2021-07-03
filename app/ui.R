library(shiny)
require(shinyjs)
library(shinythemes)
require(shinycssloaders)
library(tidyr)
library(tidyverse)
library(dplyr)
library(DT)
library(colourpicker)
library(RColorBrewer)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(msigdbr)
library(stringr)
library(pathview)

library(Seurat)
library(hdf5r)
library(cowplot)
library(SingleR)
library(scRNAseq)
library(slingshot)
library(SeuratWrappers)
library(monocle3)
library(patchwork)
library(magrittr)

dataIn <- load('./data.RData')

tagList(
    tags$head(
        #includeHTML(("www/GA.html")), If including google analytics
        tags$style(type = 'text/css','.navbar-brand{display:none;}')
    ),
    fluidPage(theme = shinytheme('yeti'),
            windowTitle = "scRNA",
            titlePanel(
                fluidRow(
                #column(4, img(height =75 , src = "")), if I ever want a logo in there
                column(12, h1(strong('scRNA'), align = 'center')),
                column(12, h2(strong('Tool'), align = 'center'))
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
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='UMAP',
                            tabPanel(title='UMAP Plot', hr()
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
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='Gene Plots',
                            tabPanel(title='Feature Plot', hr()
                            
                            ),
                            tabPanel(title='Violin Plot', hr()
                            
                            ),
                            tabPanel(title='Ridge Plot', hr()
                            
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
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='Marker Genes',
                            tabPanel(title='Tables', hr()
                            
                            ),
                            tabPanel(title='Heatmaps', hr()
                            
                            ),
                            tabPanel(title='Dot Plot', hr()
                            
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
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='Gene Signatures',
                            tabPanel(title='Cell Cycling', hr()
                            
                            ),
                            tabPanel(title='Gene Lists', hr()
                            
                            ),
                            tabPanel(title='Pseudotime', hr()
                            
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
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='Cross Comparisons',
                            tabPanel(title='Tables', hr()
                            
                            ),
                            tabPanel(title='Heatmaps', hr()
                            
                            ),
                            tabPanel(title='Violin Plots', hr()
                            
                            )
                        )
                    )
                )
            ),

            ## Idents Page
##########################################################################################################################################################
            tabPanel('Idents',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='Idents',
                            tabPanel(title='Ident relabel', hr()
                            
                            )
                        )
                    )
                )
            ),

            ## Reclustering Page
##########################################################################################################################################################
            tabPanel('Reclustering',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='Reclustering',
                            tabPanel(title='Cluster selection', hr()
                            
                            )
                        )
                    )
                )
            ),

            ## Autolabel Page
##########################################################################################################################################################
            tabPanel('Automated Labeling',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='Auto label',
                            tabPanel(title='Label plot', hr()
                            
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
                        tags$a(href="https://github.com/alemenze", icon("github", "fa-3x")),
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
