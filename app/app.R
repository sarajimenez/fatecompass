# libraries ---------------------------------------------------------------
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(DT)
library(shiny)
library(shinydashboard)
library(shinythemes)

source("global.R")

# Shiny app ---------------------------------------------------------------

ui <- fluidPage(theme = shinytheme("sandstone"),
  
  fluidRow(
    column(10,
           h1("FateCompass: for identifying lineage-specific TF activity dynamics"),
           "The FateCompass pipeline aims to estimate transcription factor (TF) activities in a dynamic fashion 
           using single-cell transcriptomics by modeling cell transitions in a probabilistic manner and the observed 
           gene expression in terms of computationally predicted regulatory sites.",
           strong("The novelty of FateCompass lies in its integrative approach that allows the identification of temporal-lineage-specific regulators."),
           "As proof of concept we used two well-characterized", strong("pancreatic islet cell formation"), "data sets from mouse",em("in vivo"),"and human", 
           em("in vitro."),
           hr(),
           p("For more info use the links below: "),
           hr()),
    
    column(2,img(src = "fatecompass.png",height = 180, width = 200))
    ),
  
  splitLayout(
    actionButton(inputId  = "paper", label = "Paper", icon(name = "scroll", lib = "font-awesome"), width = '100%', onclick ="window.open('https://doi.org/10.1101/2022.04.01.486696', '_blank')"),
    actionButton(inputId  = "code", label = "Code", icon(name = "github", lib = "font-awesome"), width = '100%', onclick ="window.open('https://github.com/sarajimenez/fatecompass', '_blank')"),
    actionButton(inputId  = "database", label = "Binding Sites DB", icon(name = "database", lib = "font-awesome"), width = '100%', onclick ="window.open('https://swissregulon.unibas.ch/sr/swissregulon', '_blank')")
  ),
  
  hr(),
  
  tabsetPanel(
  tabPanel("Mouse in vivo",
           
           hr(),
           "Data set of pancreatic endocrinogenesis at embryonic day 15.5 published by Bastidas-Ponce et al. 2019.
           Click", a(href="https://doi.org/10.1242/dev.173849","here",target="_blank"), "for the original paper",
           hr(),
           "Cell types and FateCompass predictions",
           hr(),
           
           fluidRow(
             column(4, plotOutput("umap_groups_m")),
             column(8, DT::dataTableOutput("dma_m"))
           ),
           
           hr(),
           "TF activity and gene expression profiles in UMAP plot and over differentiation trajectories",
           hr(),
           
           fluidRow(
             column(3, selectInput(inputId  = "motifid_m",
                                   label    = "Select motif of interest",
                                   choices  = motifs_m,
                                   selected = "Nkx6.1_Evx1_Hesx1")),
             column(3, br()),
             column(3, selectInput(inputId  = "geneid_m",
                                   label    = "Select a TF from the motif family",
                                   choices  = genes_m,
                                   selected = "Nkx6.1")),
             column(3, br())
           ),
           
           fluidRow(
             column(3, plotOutput("umap_act_m")),
             column(3, plotOutput("avg_act_m")),
             column(3, plotOutput("umap_exp_m")),
             column(3, textOutput("error_1")),
             column(3, plotOutput("avg_exp_m"))
           )
           ),# tabpanel - Mouse
  
  tabPanel("Human in vitro", 
           
           hr(),
           "Dataset of",em("in vitro"), "differentiation of human embryonic stem cells (hESCs) towards beta-like cells 
           published by Veres et al. 2019. Click", a(href="https://doi.org/10.1038/s41586-019-1168-5","here",target="_blank"),"for the original paper",
           hr(),
           "Cell types and FateCompass predictions",
           hr(),
                                              
           fluidRow(
             column(4, plotOutput("umap_groups_h")),
             column(8, DT::dataTableOutput("dma_h"))
           ),
           
           hr(),
           "TF activity and gene expression profiles in UMAP plot and over differentiation trajectories",
           hr(),
           
           fluidRow(
             column(3, selectInput(inputId  = "motifid_h",
                                   label    = "Select motif of interest",
                                   choices  = motifs_h,
                                   selected = "NOTO_VSX2_DLX2_DLX6_NKX6.1")),
             column(3, br()),
             column(3, selectInput(inputId  = "geneid_h",
                                   label    = "Select a TF from the motif family",
                                   choices  = genes_h,
                                   selected = "NKX6.1")),
             column(3, br())
           ),
           
           fluidRow(
             column(3, plotOutput("umap_act_h")),
             column(3, plotOutput("avg_act_h")),
             column(3, plotOutput("umap_exp_h")),
             column(3, textOutput("error_3")),
             column(3, plotOutput("avg_exp_h"))
           )
           )# tabpanel - human
  
)#tabsetpanel

)#fluidpage

server <- function(input, output) {
  
  # Mouse tab
    observe({
      tfs_m <- unlist(strsplit(as.character(input$motifid_m),"_"))
      updateSelectInput(session = getDefaultReactiveDomain(),
                        inputId = "geneid_m",
                        label = NULL,
                        choices = tfs_m,
                        selected = "Nkx6.1")
    })
    output$umap_groups_m <- renderPlot({
      cell_type <- factor(i_m$V1)
      ggplot(u2_m, aes(x=umap1, y=umap2, colour = cell_type, fill = cell_type)) +
        geom_point(shape = 21, alpha = 0.5, size = 1, aes(colour = cell_type)) + 
        scale_colour_manual(values = cols_m, aesthetics = c("colour", "fill")) +
        theme_void()
    })
  
    output$dma_m <- DT::renderDataTable(dma_m,
                                      options = list(scrollX = TRUE),
                                      rownames = FALSE)
  
    output$umap_act_m <- renderPlot({
      motif_color <- a_m[,input$motifid_m]
      general_title_act <- "activity"
      mid <- 0
      ggplot(u2_m, aes(x=umap1, y=umap2, color = motif_color)) +
        geom_point(size=1) + 
        scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = mid) +
        ggtitle(paste(input$motifid_m, general_title_act)) +
        theme_void()
    })
    output$avg_act_m <- renderPlot({
      mpg_act <- data.frame(cbind(seq(0,1,length.out=2000),act_traj_alpha_m[,input$motifid_m],act_traj_beta_m[,input$motifid_m]))
      colnames(mpg_act) <- c("Simulated_time", "Alpha", "Beta")
      ggplot() + 
        geom_smooth(data = mpg_act, aes(x = Simulated_time, y = Alpha), colour = "#3477B4") +
        geom_smooth(data = mpg_act, aes(x = Simulated_time, y = Beta), colour = "#B3E08B") +
        xlab("Simulated time") +
        ylab("Average motif activity") +
        ggtitle(input$motifid_m) +
        theme_classic()
    })
    output$umap_exp_m <- renderPlot({
      if(input$geneid_m%in%colnames(g_m)){
        gene_color <- g_m[,input$geneid_m]
        general_title_exp <- "gene expression"
        mid <- max(gene_color)/2
        ggplot(u2_m, aes(x=umap1, y=umap2, color = gene_color)) +
          geom_point(size=1) + 
          scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = mid) +
          ggtitle(paste(input$geneid_m, general_title_exp)) +
          theme_void()
      }
    })
    output$error_1 <- renderText({
      if(!(input$geneid_m%in%colnames(g_m))){
        error_msj <- paste("The gene",input$geneid_m,"is not expressed")
        print(error_msj)
      }
    })
    output$avg_exp_m <- renderPlot({
      if(input$geneid_m%in%colnames(g_m)){
        mpg_exp <- data.frame(cbind(seq(0,1,length.out=2000),exp_traj_alpha_m[,input$geneid_m],exp_traj_beta_m[,input$geneid_m]))
        colnames(mpg_exp) <- c("Simulated_time", "Alpha", "Beta")
        ggplot() + 
          geom_smooth(data = mpg_exp, aes(x = Simulated_time, y = Alpha), colour = "#3477B4") +
          geom_smooth(data = mpg_exp, aes(x = Simulated_time, y = Beta), colour = "#B3E08B") +
          xlab("Simulated time") +
          ylab("Average gene expression") +
          ggtitle(input$geneid_m) +
          theme_classic()
      }
    })
    
    # Human tab
    observe({
      tfs_h <- unlist(strsplit(as.character(input$motifid_h),"_"))
      updateSelectInput(session = getDefaultReactiveDomain(),
                        inputId = "geneid_h",
                        label = NULL,
                        choices = tfs_h,
                        selected = "NKX6.1")
    })
    output$umap_groups_h <- renderPlot({
      cell_type <- factor(i_h$Labels)
      ggplot(u2_h, aes(x=umap1, y=umap2, colour = cell_type, fill = cell_type)) +
        geom_point(shape = 21, alpha = 0.3, size = 1, aes(colour = cell_type)) + 
        scale_colour_manual(values = cols_h, aesthetics = c("colour", "fill")) +
        theme_void()
    })
    
    output$dma_h <- DT::renderDataTable(dma_h,
                                        options = list(scrollX = TRUE),
                                        rownames = FALSE)
    
    output$umap_act_h <- renderPlot({
      motif_color <- a_h[,input$motifid_h]
      general_title_act <- "activity"
      mid <- 0
      ggplot(u2_h, aes(x=umap1, y=umap2, color = motif_color)) +
        geom_point(size=1) + 
        scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = mid) +
        ggtitle(paste(input$motifid_h, general_title_act)) +
        theme_void()
    })
    output$avg_act_h <- renderPlot({
      mpg_act <- data.frame(cbind(seq(0,1,length.out=1000),act_traj_alpha_h[,input$motifid_h],
                                  act_traj_beta_h[,input$motifid_h],act_traj_ec_h[,input$motifid_h]))
      colnames(mpg_act) <- c("Simulated_time", "Alpha", "Beta","EC")
      ggplot() + 
        geom_smooth(data = mpg_act, aes(x = Simulated_time, y = Alpha), colour = "#3477B4") +
        geom_smooth(data = mpg_act, aes(x = Simulated_time, y = Beta), colour = "#B3E08B") +
        geom_smooth(data = mpg_act, aes(x = Simulated_time, y = EC), colour = "#EE847C") +
        xlab("Simulated time") +
        ylab("Average motif activity") +
        ggtitle(input$motifid_h) +
        theme_classic()
    })
    output$umap_exp_h <- renderPlot({
      if(input$geneid_h%in%colnames(g_h)){
        gene_color <- g_h[,input$geneid_h]
        general_title_exp <- "gene expression"
        mid <- max(gene_color)/2
        ggplot(u2_h, aes(x=umap1, y=umap2, color = gene_color)) +
          geom_point(size=1) + 
          scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = mid) +
          ggtitle(paste(input$geneid_h, general_title_exp)) +
          theme_void()
      }
    })
    output$error_3 <- renderText({
      if(!(input$geneid_h%in%colnames(g_h))){
        error_msj <- paste("The gene",input$geneid_h,"is not expressed")
        print(error_msj)
      }
    })
    output$avg_exp_h <- renderPlot({
      if(input$geneid_h%in%colnames(g_h)){
        mpg_exp <- data.frame(cbind(seq(0,1,length.out=1000),exp_traj_alpha_h[,input$geneid_h],
                                    exp_traj_beta_h[,input$geneid_h],exp_traj_ec_h[,input$geneid_h]))
        colnames(mpg_exp) <- c("Simulated_time", "Alpha", "Beta","EC")
        ggplot() + 
          geom_smooth(data = mpg_exp, aes(x = Simulated_time, y = Alpha), colour = "#3477B4") +
          geom_smooth(data = mpg_exp, aes(x = Simulated_time, y = Beta), colour = "#B3E08B") +
          geom_smooth(data = mpg_exp, aes(x = Simulated_time, y = EC), colour = "#EE847C") +
          xlab("Simulated time") +
          ylab("Average gene expression") +
          ggtitle(input$geneid_h) +
          theme_classic()
      }
    })
}

shinyApp(ui = ui, server = server)