
options("scipen"=1, "digits"=4)
library(shiny)
library(dplyr)
library(tidyr)
library(data.table)
library(plotly)
library(ggplot2)
library(d3heatmap)
library(ggfortify)

converted_dt <- read.table("data/converted.csv",
                           sep = "|",
                           header = TRUE,
                           stringsAsFactors = FALSE)


tables_info <- read.table("data/tables.txt",
                          sep = "|",
                          header = TRUE,
                          stringsAsFactors = FALSE)

agg_methods <- c(
  "Raw Data",
  "Mean (columnwise)",
#  "Median (columnwise)",
  "Max (rowwise)")

shinyUI(fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      selectInput("dt",
                  "Dataset:",
                  c("None", tables_info$names)),
      br(),
      
      selectInput("agg",
                         "Aggregation method:",
                         agg_methods),
      br(),
      helpText("Mean method outputs an average value", 
               "for each column for each unique gene name."),
      br(),
      helpText("Max method outputs the values of the row", 
               "with the maximum sum of values among the", 
               "rows with the same gene name.")
    ),

    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Info", includeMarkdown("info.Rmd")),
                  tabPanel("Table", DT::dataTableOutput("table")),
                  tabPanel("t-test", 
                           fluidRow(
                             column(3, 
                                    selectizeInput(
                                      'cancer_types', 'Choose two types of cell lines', 
                                      choices = unique(converted_dt$type),
                                      multiple = TRUE, options = list(maxItems = 2)
                           )),# 1. cancertype
                             column(3, 
                                    radioButtons("corr_type", 
                                          "Correction Type",
                                         choices = list("None" = "none", 
                                                        "Bonferroni" = "bonferroni",
                                                        "Benjamini-Hochberg (FDR)" = "BH"), 
                                         selected = "none")), # 2. correction
                             column(3,
                                    radioButtons("eq_var", 
                                                 "Equal Variance",
                                                 choices = list("False" = 1, 
                                                                "True" = 2), 
                                                 selected = 2)),
                           column(3, 
                                  sliderInput("thresh", "Significance Threshold", 
                                              min = 0, max = 1, value = 0.5, step = 0.01))),
                           DT::dataTableOutput("ttest")),
                  tabPanel("Exploratory", 
                           fluidRow(
                             column(3,
                                    selectInput("plot_type",
                                                "Plot Type",
                                                c("box", "histogram"))),
                             column(3, 
                                    selectizeInput(
                                      'cell_lines', 'Choose cell lines (max. 10)', 
                                      choices = sort(unique(converted_dt$mapping)),
                                      multiple = TRUE, options = list(maxItems = 10)
                                    ))),
                           plotlyOutput("plot_ly")),
                  tabPanel("Heatmap",
                           fluidRow(column(3,
                                           radioButtons("comp_method", 
                                                        "Compression Method",
                                                        choices = list("PCA" = 1,
                                                                       "Random Sampling" = 2), 
                                                        selected = 1)),
                                    column(3,
                                           sliderInput("hm_sample", "Number of genes (only for random sampling)", 
                                       min = 0, max = 200000, value = 400, step = 100))),
                           d3heatmapOutput("h_map")),
                           #plotOutput("h_map")),
                  #tabPanel("Alt. Clustering",
                  #         plotlyOutput("pca"))
                  tabPanel("Dendrogram",
                           radioButtons("comp_method2",
                                        "Compression Method",
                                        choices = list("PCA" = 1,
                                                       "None" = 2), 
                                                        selected = 1),#),
                                    #column(3,
                                    #       sliderInput("hm_sample2", "Number of genes", 
                                    #                   min = 0, max = 200000, value = 400, step = 100))),
                           plotOutput("dendr"))
                  
      )
    ))
  
))
