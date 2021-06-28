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


tagList(
    tags$head(
        #includeHTML(("www/GA.html")), If including google analytics
        tags$style(type = 'text/css','.navbar-brand{display:none;}')
    ),
    fluidPage(theme = shinytheme('superhero'),
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
                            column(8,markdown("Text 
                            ")),
                            column(2,)
                            
                        ),
                    column(2,
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
