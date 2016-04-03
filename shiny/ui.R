
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
options("scipen"=1, "digits"=4)
library(shiny)
library(dplyr)
library(tidyr)
library(data.table)

#available_tables <- c("RNA AFFY HG U-95 (A-E)",
#                      "RNA AFFY HG U-133 (A-B)",
#                      "RNA AFFY HG U-133 Plus 2.0")

tables_info <- read.table("data/tables.txt",
                          sep = "|",
                          header = TRUE,
                          stringsAsFactors = FALSE)

agg_methods <- c(
  "Raw Data",
  "Mean (columnwise)",
#  "Median (columnwise)",
  "Max (rowwise)")


# Define the overall UI
shinyUI(
  fluidPage(
    titlePanel("Bio Data Tables"),
    
    # Create a new Row in the UI for selectInputs
    fluidRow(
      column(6,
             selectInput("dt",
                         "Choose Dataset:",
                         c("None", tables_info$names))
      ),
      column(6,
             selectInput("agg",
                         "Choose the aggregation method:",
                         agg_methods)
      )
    ),
    # Create a new row for the table.
    fluidRow(
      DT::dataTableOutput("table")
    )
  )
)
