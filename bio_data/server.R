

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

shinyServer(function(input, output) {
  
  df <- reactive({
    if (input$dt != "None") {
      data_type <- tables_info %>%
        filter(., names == input$dt) %>% .$data_type
      
      if (data_type == "gene") {
        df <- fread(tables_info %>%
                      filter(., names == input$dt) %>% .$files,
                    sep = "\t", header = TRUE,
                    select = c(2, 10:69),
                    nrows = -1,
                    stringsAsFactors = FALSE, na.strings = "-")
        
        setnames(df, colnames(df)[1], "gene_name")
        setnames(df, colnames(df)[-1], 
                 data.table(converted_dt)[match(colnames(df)[-1] %>% 
                                            strsplit(., ":") %>% 
                                            lapply(., function(x) x[2]) %>% 
                                            unlist,
                                            gene)]$mapping)
        
        setcolorder(df, c(colnames(df)[1], sort(colnames(df)[-1])))
        
        
      }else if (data_type == "proteomes") {
        
        df <- fread(tables_info %>%
                      filter(., names == input$dt) %>% .$files,
                    select = c(1, 7, 28:86), # 7 is a fasta header
                    nrows = -1)
        
        df <- merge(df,
                    data.table(good_symbols_acc),
                    by = "Accession")[, !"Accession", with = FALSE]
        
        fasta_headers <- df[grep("^[0-9]+-", CleansedSymbol)]$`Fasta headers`
        gene_symb_fasta <- gsub(".*Gene_Symbol=([[:alnum:]-]+) .*",
                                '\\1',
                                fasta_headers)
        df[grep("^[0-9]+-", CleansedSymbol)]$CleansedSymbol <- gene_symb_fasta
        
        df[, `Fasta headers`:= NULL]
        
        
        setcolorder(df, c("CleansedSymbol", setdiff(names(df), "CleansedSymbol")))
        setnames(df, colnames(df)[1], "gene_name")
        
        setnames(df, colnames(df)[-1], 
                 colnames(df)[-1] %>% gsub("LFQ_(.*)", '\\1', .))
        
        #setnames(df, colnames(df)[-1], 
        #         unique(data.table(converted_dt)[,.[match(colnames(df)[-1], unique(proteoms))]$mapping)
        setcolorder(df, c(colnames(df)[1], sort(colnames(df)[-1])))
        
        
      }
      data_type <- tables_info %>%
        filter(., names == input$dt) %>% .$data_type

      type <- input$agg
      
      # drop NA columns
      df <- df[,which(unlist(lapply(df, 
                                    function(x) !all(is.na(x))))),with=F]
      
      df <- df %>%
        filter(., gene_name != "", gene_name != " ", gene_name != "_")
      if (type != "Raw Data") {
        df <- aggregate_df(df, type)
      }
      df
    }
  })

  output$table <- DT::renderDataTable(DT::datatable({
    
    
    if (length(df()) != 0)  
      
      
      format(df(), digits = 4)
        
    
    
    
    
  }, 
  extensions = c('Buttons', 'FixedColumns'),
  options = list(pageLength = 50, 
                 dom = 'Bfrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print',I('colvis')),
                 scrollX = TRUE,
                 fixedColumns = list(leftColumns = 1))#,
  #filter = 'top'
  ))
  
  # T test
  output$ttest <- DT::renderDataTable(DT::datatable({
    
    
    if ((length(input$cancer_types) == 2) 
        & (input$dt != "None") 
        & (input$agg != "Raw Data")) {
      
      if (input$eq_var == 1)
        eq_var <- FALSE
      else
        eq_var <- TRUE
      
      df_t <- df() %>%
        select(., gene_name, 
               contains(input$cancer_types[1]),
               contains(input$cancer_types[2]))
      
      col_type1 <- grep(input$cancer_types[1], colnames(df_t))
      col_type2 <- grep(input$cancer_types[2], colnames(df_t))
      
      df_dt <- as.data.frame(df_t)
      rownames(df_dt) <- df_dt$gene_name
      df_dt$gene_name <- NULL
      
      pval_arr <- apply(df_dt, 1, function(x) {
        tryCatch(t.test(x[col_type1], x[col_type2], 
                        var.equal = eq_var)$p.value,
                 error = function(x) NA)
      })
      
      data.table("gene_name" = names(pval_arr),
                 "p.value" = p.adjust(pval_arr, method = input$corr_type)) %>%
        filter(., p.value < input$thresh)
      
      #df_t %>%
      #  select(., gene_name) %>%
      #  cbind(., p.value = apply(copy(df_t)[,gene_name:=NULL], 1, 
      #                           function(x) 
      #                             t.test(x[col_type1], x[col_type2])$p.value))
      #df_t[, t.test(.SD[, col_type1, with = F],
      #      .SD[, col_type2, with = F])$p.value, by = "gene_name"]
      #df_t %>%
      #  select(., col_type2, col_type1)
      #df_t[col_typ1]
        #mutate(., sum1 = sum(select(., col_type1)), 
        #       sum2 = sum(select(., col_type2)))
        #mutate(p.value = t.test(col_type1, col_type2)$p.value)
      
      #pvalues <- apply(df_t, 1, function(x) {
      #  t.test(x[col_type1], x[col_type2], na.action = "na.omit")
      #})
        
      #data.table(str(df_t))
      #format(df_t, digits = 4)
      
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
  
  output$conversion <- DT::renderDataTable(DT::datatable({
  
    read.table("data/converted_final.csv",
               sep = "|",
               quote = "",
               header = TRUE,
               stringsAsFactors = FALSE)
    },
    extensions = c('Buttons'),
    options = list(pageLength = 72, 
                   dom = 'Bfrtip2',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print',I('colvis')),
                   scrollX = TRUE)
    
    
  ))
  
  output$plot_ly <- renderPlotly({
    if ((input$dt != "None") & (length(input$cell_lines) >= 1)) {
      
      chosen_cell_lines <- match(input$cell_lines, colnames(df()))
      #chosen_cell_lines <- sample(seq_len(length(colnames(df())))[-1], 10)
      
      df_t <- df() %>%
        select(., chosen_cell_lines) %>%
        gather(., key = "cell_line", value = "gene_expression")
      
      
      
      df_t %>%
        filter(cell_line %in% chosen_cell_lines)
      
      a <- list(
        autotick = FALSE,
        #ticks = "outside",
        #tick0 = 0,
        #dtick = 0.25,
        #ticklen = 8,
        #tickwidth = 4,
        tickfont = list(size = 8)
      )
      
      m = list(
        l = 50,
        r = 50,
        b = 200,
        t = 60,
        pad = 4
      )
      
      if (input$plot_type == "box") {
      
      plot_ly(df_t, 
              y = gene_expression, 
              color = cell_line, 
              type = input$plot_type) %>%
        layout(autosize = F, 
               xaxis = a, 
               width = 800, height = 700, margin = m)
      }else if (input$plot_type == "histogram"){
        
        
        plot_ly(df_t, 
                x = gene_expression, 
                color = cell_line, 
                type = input$plot_type,
                opacity = 0.6) %>%
          layout(autosize = F, 
                 xaxis = a, 
                 width = 800, height = 700, margin = m,
                 barmode = "overlay")
        
      }
      

        #layout(autosize = F, width = 1000, height = 800)
      #p <- df_t %>%
      #  ggplot(., aes(cell_line, gene_expression)) +
      #  geom_boxplot() +
      #  theme(axis.text.x=element_text(angle = -90, hjust = 0))
      
      #p <- ggplotly(p)
      #p
              
      
    }else
      return();

    
  })
  
  output$h_map <- renderD3heatmap({
  #output$h_map <- renderPlot({  
    # size of the bins depend on the input 'bins'
    if (input$dt != "None") {
      
      # 1 PCA
      # 2 Random Sampling
      
      if (input$comp_method == 2) {
        
        df_t <- df() %>%
          sample_n(., input$hm_sample) %>%
          select(., -gene_name) %>%
          t()
        
      }else{
        
        df_t <- df() %>%
          select(., -gene_name) %>%
          t() %>%
          prcomp(.)
        
        df_t <- df_t$x
        
      }
      
      df_t <- t(df_t)
        #gather(., key = "cell_line", value = "gene_expression")
      dist2 <- function(x, ...)
        as.dist(t(1-cor(t(x), method="pearson")))
    
      hc <- hclust(as.dist(t(1 - cor(df_t))))
      #plot(hc)
    d3heatmap(df_t, scale = "column", colors = "Spectral",
              distfun = dist2,
              Rowv = NULL,
              labRow = NULL,
              width = 700, height = 100,
              xaxis_height = 150,
              show_grid = FALSE)
    }
  })
  
  #output$h_map <- renderD3heatmap({
  output$dendr <- renderPlot({  
    if (input$dt != "None") {
      
      # 1 PCA
      # 2 Random Sampling
      
      if (input$comp_method2 == 2) {
        
        df_t <- df() %>%
          #sample_n(., input$hm_sample2) %>%
          select(., -gene_name) %>%
          t()
        
      }else{
        
        df_t <- df() %>%
          select(., -gene_name) %>%
          t() %>%
          prcomp(.)
        
        df_t <- df_t$x
        
      }
      
      df_t <- t(df_t)
      #gather(., key = "cell_line", value = "gene_expression")
      dist2 <- function(x, ...)
        as.dist(t(1-cor(t(x), method="pearson")))
      
      hc <- hclust(as.dist(t(1 - cor(df_t))))
      plot(hc, ann = FALSE, axes = FALSE)
      
    }
  })  
  
  output$pca <- renderPlotly ({#renderPlot({  
    if (input$dt != "None") {
      
      # 1 PCA
      # 2 Random Sampling
      
      #if (input$comp_method == 2) {
      #  
      #  df_t <- df() %>%
      #    sample_n(., input$hm_sample) %>%
      #    select(., -gene_name) %>%
      #    t()
      #  
      #}else{
        
      #  df_t <- df() %>%
      #    select(., -gene_name) %>%
      #    t() %>%
      #    prcomp(.)
      #  
      #  df_t <- df_t$x
        
      #}
      
      df_m <- df() %>%
        select(-gene_name) %>%
        t()
      
      df_t <- as.data.frame(df_m)
      df_t$cell_line_full <- row.names(df_m)
      df_t$cell_line <- row.names(df_m) %>%
        strsplit("_") %>%
        lapply(., function(x) x[1]) %>% unlist()
      
      #print(df_t)
      #gather(., key = "cell_line", value = "gene_expr")
      #print (rownames(t(df_t)))
      p <- df_t %>%
        select(-cell_line, -cell_line_full) %>%
        prcomp(.) %>%
        autoplot(., data = df_t, colour = "cell_line")
      ggplotly(p)
    }
  })
  
  
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
