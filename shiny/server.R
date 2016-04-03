
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

options("scipen"=1, "digits"=4)

library(shiny)
library(dplyr)
library(tidyr)
library(data.table)

converted_dt <- read.table("data/converted.csv",
                           sep = "|",
                           header = TRUE,
                           stringsAsFactors = FALSE)

tables_info <- read.table("data/tables.txt",
                          sep = "|",
                          header = TRUE,
                          stringsAsFactors = FALSE)

summary_df <- read.csv("data/summary_proteome.csv",
                       stringsAsFactors = FALSE)

blank_symbols_acc <- summary_df %>% filter(., Symbol == "") %>% .$Accession %>%
  unique

summary_df$CleansedSymbol <- lapply(summary_df$Symbol %>% 
                                      strsplit(" "), 
                                    function(x) x[1]) %>% 
  unlist

good_symbols_acc <- summary_df %>%
  filter(., !is.na(CleansedSymbol)) %>%
  select(., Accession, CleansedSymbol) %>%
  unique


# Define a server for the Shiny app
shinyServer(function(input, output) {
  
  # Filter data based on selections
  output$table <- DT::renderDataTable(DT::datatable({
    
    
    if (input$dt != "None") {
      data_type <- tables_info %>%
        filter(., names == input$dt) %>% .$data_type
      
      if (data_type == "gene") {
        df <- fread(tables_info %>%
                          filter(., names == input$dt) %>% .$files,
                  sep = "\t", header = TRUE,
                  select = c(2, 10:69),
                  stringsAsFactors = FALSE, na.strings = "-")
      
        setnames(df, colnames(df)[1], "gene_name")
      }else if (data_type == "proteomes") {
        
        df <- fread(tables_info %>%
                      filter(., names == input$dt) %>% .$files,
                    select = c(1, 28:86))
        
        df <- merge(df,
                    data.table(good_symbols_acc),
                    by = "Accession")[, !"Accession", with = FALSE]
        
        setcolorder(df, c("CleansedSymbol", setdiff(names(df), "CleansedSymbol")))
        setnames(df, colnames(df)[1], "gene_name")
        
      }else if (data_type == "deep_proteomes") {
        
        df <- fread(tables_info %>%
                      filter(., names == input$dt) %>% .$files,
                    select = c(1, 28:36))
        
        df <- merge(df,
                    data.table(good_symbols_acc),
                    by = "Accession")[, !"Accession", with = FALSE]
        
        setcolorder(df, c("CleansedSymbol", setdiff(names(df), "CleansedSymbol")))
        setnames(df, colnames(df)[1], "gene_name")
        
      }
      #df <- data.frame(df)
      
      type <- input$agg
      
      df <- df %>%
        filter(., gene_name != "", gene_name != " ", gene_name != "_")
      if (type != "Raw Data") {
        df <- aggregate_df(df, type)
      }
      
      format(df, digits = 4)
        
    }
    
    
    
    
  }, 
  extensions = c('Buttons', 'FixedColumns'),
  options = list(pageLength = 50, 
                 dom = 'Bfrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print',I('colvis')),
                 scrollX = TRUE,
                 fixedColumns = list(leftColumns = 1))#,
  #filter = 'top'
  ))
  
})


aggregate_df <- function(df, type) {
  
  if (length(grep("Mean", type)) == 1) {
    method <- "mean"
    df[order(gene_name)
       ][, lapply(.SD, eval(parse(text = method)), na.rm = TRUE), by = "gene_name"]
  } else if (length(grep("Median", type)) == 1) {
    method <- "median"
  } else {
    method <- "max"
    df$rs <- rowSums(df %>% select(., -gene_name), na.rm = TRUE)
    #df <- df[,.SD[which.max(rs)], by = gene_name]
    #df[, rs:= NULL]
    df[order(-rs)
       ][, .SD[1], by = gene_name
         ][order(gene_name)
           ][, !"rs", with = FALSE]
    
  }
  #df_by_gene <- df %>%
  #  group_by(., gene_name)
  
  #df_by_gene %>%
  #  summarise_each(., eval(parse(text = paste0("funs(", method, "(., na.rm = FALSE))")))) %>%
  #  ungroup()
  #df[order(gene_name)
  #   ][, lapply(.SD, eval(parse(text = method))), by = "gene_name"]
  
  
}
