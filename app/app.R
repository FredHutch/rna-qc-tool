library(shiny)
library(Seurat)
library(ggplot2)
library(plotly)
library(bslib)
library(glue)
library(reactable)
# library(webshot2)
# library(bs4Dash)

options(shiny.maxRequestSize = 1*1024^3)  # max upload size 1 GB



MIN_GENES <- function(Sobj, N){
  gene_counts <- Matrix::rowSums(Sobj[["RNA"]]$counts > 0)
  Sobj2 <- subset(Sobj, features = names(gene_counts[gene_counts >= N]))
  return(Sobj2)
}

# Finding cumulative median
CUM_MEDIAN <- function(x) {        
  sapply(1:length(x), function(i) median(x[1:i]))
}


# source: 1. https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/handling-sparsity.html
#         2. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1926-6#Sec7
SPARSITY <- function(Sobj){                        # returns proportion of genes that are not expressed in >80% of cells
  counts_matrix <- Sobj[["RNA"]]$counts
  total_genes <- dim(counts_matrix)[1]
  total_cells <- dim(counts_matrix)[2]
  zero_rate <- Matrix::rowSums(counts_matrix == 0) / total_cells  # proportion of cells where the gene is not expressed
  high_zero_rate_genes <- names(zero_rate[zero_rate > 0.80]) # identify genes not expressed in >80% of cells
  return(length(high_zero_rate_genes)/total_genes)
}


ui <- fluidPage(
  tags$head(includeHTML("google-analytics.html")),
  titlePanel(strong("Single-Cell RNA Quality Control Tool")),
  h5(p("An application for a visual and interactive quality control of scRNA-seq data")),
  hr(),
  
  page_sidebar(
    sidebar = sidebar(
      fileInput("file", "Upload HDF5 File (.h5)", accept = c(".h5")),
      sliderInput("mt_threshold", "Mitochondrial Gene Threshold (%)", min = 0, max = 100, value = 10),
      sliderInput("rb_threshold", "Ribosomal Gene Threshold (%)", min = 0, max = 100, value = 45),
      card( card_header('Features/cell'),
            card_body( layout_column_wrap(width = 0.5, gap = "5px",
                                          numericInput("nFeature_min", "Minimum", value = 200, min = 0),
                                          numericInput("nFeature_max", "Maximum", value = 10e6, min = 0)
                                          )
                       ),
            ),
      
      card( card_header('UMIs/cell'),
            card_body( layout_column_wrap(width = 0.5, gap = "5px",
                                          numericInput("nCount_min", "Minimum", value = 500, min = 0, ),
                                          numericInput("nCount_max", "Maximum", value = 10e6, min = 0 ),
                                          )
            )
      ),
      # numericInput("min_gene_cnt", "Minimum number of cells where a gene should be present", value = 3, min = 0),
      checkboxInput("density_plot", "View density plot", value = FALSE),
      actionButton("applyFilter", "Apply Filter"),
      card(
        card_header('Download script'),
        card_body(downloadButton("downloadPython", "Python"),
                  downloadButton("downloadR", "R"))
      ),
      # downloadButton("saveReport", "Download Report")
    ),
    
    navset_card_underline(
      # title = "",
      
      nav_panel("Table", column(4,reactableOutput("seuratSummary"))),
      
      nav_panel("Plots",
                fluidRow(
                  column(6, plotlyOutput("dist_umi")),
                  column(6, plotlyOutput("dist_feature"))
                ),
                
                
                fluidRow( column(4, plotlyOutput("plot_mt")),
                          column(4, plotlyOutput("plot_rb")),
                          column(4, plotlyOutput("plot_count_vs_feature"))
                ),
                fluidRow(column(4, plotlyOutput("knee_plot")),
                         column(4, plotlyOutput("novelty_score")),
                         column(4, plotlyOutput("gene_sensitivity"))
                         
                )
                
      )
    ),
    # card(
    #   card_header('Cardhead')
    # )
    
    
  )
)


server <- function(input, output, session) {
  
  # Reactive value to store the original Seurat object
  original_seurat <- reactiveVal()
  
  # Observe file upload and load the Seurat object
  observeEvent(input$file, {
    req(input$file)
    infile <- input$file$datapath
    seurat_obj <- CreateSeuratObject(counts = Read10X_h5(infile))
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
    seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
    original_seurat(seurat_obj)
  })
  
  # Reactive expression to apply filtering based on user input
  filtered_seurat <- reactive({ req(original_seurat())
    seurat_copy <- original_seurat()
    seurat_copy <- subset(seurat_copy, subset =  nFeature_RNA > input$nFeature_min & nFeature_RNA < input$nFeature_max &
                                                nCount_RNA > input$nCount_min & nCount_RNA < input$nCount_max &
                                                percent.mt < input$mt_threshold &
                                                percent.rb < input$rb_threshold)
    # DefaultAssay(seurat_copy) <- "RNAsub"
    # seurat_copy2 <- MIN_GENES(seurat_copy, input$min_gene_cnt)
    return(seurat_copy)
  })
  
  # Observe the applyFilter button and generate the summary and plots
  observeEvent(input$applyFilter, {
    
    
    output$seuratSummary <- renderReactable({ req(filtered_seurat())
      seurat <- filtered_seurat()
      metadata <- seurat@meta.data
      
      umi_counts <- round(sum(metadata$nCount_RNA))
      median_umi_counts <- round(median(metadata$nCount_RNA))
      total_cells <- round(nrow(metadata))
      total_features <- round(nrow(seurat))
      sparsity <- round(SPARSITY(seurat),2)
      
      summary_table <- data.frame( Metric = c("Total UMI Counts", "Median UMI Counts", "Total Cells", "Total Features", "Sparsity"),
                                   Value = c(umi_counts, median_umi_counts, total_cells, total_features, sparsity) )
      summary_table <- reactable(summary_table)
      
      return(summary_table)
    })
    
    
    ## PLOTS
    output$plot_mt <- renderPlotly({ req(filtered_seurat())
      p <- ggplot(filtered_seurat()@meta.data, aes(x = '', y = percent.mt)) + geom_violin(show.legend = FALSE, fill = "orange") +
        labs(title = '% of MT genes in cells', x='', y = '% MT genes') + theme_bw()
      p1 <- ggplotly(p) %>% layout(showlegend = FALSE)
      p1
    })
    
    output$plot_rb <- renderPlotly({ req(filtered_seurat())
      p <- ggplot(filtered_seurat()@meta.data, aes(x = '', y = percent.rb)) + geom_violin(show.legend = FALSE, fill = "cyan") +
        labs(title = '% of RB genes in cells', x='', y = '% RB genes') + theme_bw()
      p1 <- ggplotly(p) %>% layout(showlegend = FALSE)
      p1
    })
    
    output$plot_count_vs_feature <- renderPlotly({ req(filtered_seurat())
      p <- ggplot(filtered_seurat()@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + geom_point(size=0.1) + theme_classic() + labs(color = '%MT genes')
      # geom_vline(xintercept=input$nCount_min) + geom_hline(yintercept=input$nFeature_min) 
      p1 <- ggplotly(p) %>% layout(showlegend = FALSE)
      p1
    })
    
    
    output$dist_umi <- renderPlotly({ req(filtered_seurat())
      p <- ggplot(filtered_seurat()@meta.data, aes(x = nCount_RNA)) + theme_minimal()
      if (input$density_plot){
        p <- p + geom_density(color='purple') + labs(title = "Density Distribution of UMI counts/cell", y='Density', x='log 10 (UMI counts per cell)') + scale_x_log10()
      }
      else {
        p <- p + geom_histogram(bins = 100, fill = "purple", color = "black") + labs(title = "Frequency Distribution of UMI counts/cell", x = "UMI counts per cell", y = "Frequency")
      }
      p
    })
    
    output$dist_feature <- renderPlotly({ req(filtered_seurat())
      p <- ggplot(filtered_seurat()@meta.data, aes(x = nFeature_RNA)) + theme_minimal()
      if (input$density_plot){
        p <- p + geom_density(color='lightgreen') + labs(title = "Density Distribution of Feature counts/cell", y='Density', x='log10(Feature counts per cell)') + scale_x_log10()
      }
      else{
        p <- p + geom_histogram(bins = 100, fill = "lightgreen", color = "black") + labs(title = "Frequency Distribution of Feature counts/cell", x = "Feature counts per cell", y = "Frequency")
      }
      p
    })
    
    
    output$knee_plot <- renderPlotly({ req(filtered_seurat())
      umi_counts <- sort(filtered_seurat()@meta.data$nCount_RNA, decreasing = TRUE)
      knee_data <- data.frame( Rank = 1:length(umi_counts),UMIs = umi_counts )
      ggplot(knee_data, aes(x = Rank, y = UMIs)) +
        geom_line(color = "blue") + scale_x_log10() + scale_y_log10() +
        labs(title = "Knee Plot", x = "Barcode Rank (log10)", y = "Number of UMIs (log10)") + theme_minimal()
    })
    
    
    # source: https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html
    # number of genes detected per UMI - gives us an idea of the complexity of out dataset (more genes detected per UMI => more complex out data)
    # sobj$log10GenesPerUMI <- log10(sobj$nFeature_RNA) / log10(sobj$nCount_RNA)
    
    plot_novelty <- reactive({
      req(filtered_seurat())
      ggplot(filtered_seurat()@meta.data, aes(x=log10GenesPerUMI)) + geom_density(alpha=0.2, fill='red') +
        labs(title = "Novelty Score (Complexity)") + theme_classic()
    })
    
    
    output$novelty_score <- renderPlotly({
      p1 <- ggplotly(plot_novelty()) %>% layout(showlegend = FALSE)
      p1
    })
    
    output$gene_sensitivity <- renderPlotly({ req(filtered_seurat())
      meta <- filtered_seurat()@meta.data
      meta <- meta[order(meta$nCount_RNA, decreasing = FALSE),]
      meta$cum_median_gene <- CUM_MEDIAN(meta$nFeature_RNA)
      p <- ggplot(meta, aes(x = nCount_RNA, y = cum_median_gene)) + geom_line() +
        labs(title = 'Gene Sensitivity', x = 'Sequencing Depth (Reads/Cell)', y = 'Median Genes per Cell') + theme_classic()
      p1 <- ggplotly(p) %>% layout(showlegend = FALSE)
      p1
    })
  })
  
  ## DOWNLOADS
  output$downloadR <- downloadHandler(
    filename = function() { paste("Seurat_scRNA_", Sys.Date(), ".R", sep="")},
    content = function(file){
      scriptTemplate <- reactive({
        glue("
                              library(Seurat)
                              
                              # Load input file
                              infile <- _________
                              seurat_obj <- CreateSeuratObject(counts = Read10X_h5(infile))
                              
                              seurat_obj[[\"percent.mt\"]] <- PercentageFeatureSet(seurat_obj, pattern = \"^MT-\")
                              seurat_obj[[\"percent.rb\"]] <- PercentageFeatureSet(seurat_obj, pattern = \"^RPS|^RPL\")
                              
                              seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > {input$nFeature_min} & nCount_RNA > {input$nCount_min} &
                                                                        percent.mt < {input$mt_threshold} & percent.rb < {input$rb_threshold} &
                                                                        nFeature_RNA < {input$nFeature_max} & nCount_RNA < {input$nCount_max}
                                                                        
                              ") })
      writeLines(scriptTemplate(), file)
    }
  )
  
  output$downloadPython <- downloadHandler(
    filename = function() { paste("Scanpy_scRNA", Sys.Date(), ".py", sep="")},
    content = function(file){
      scriptTemplate <- reactive({
        glue("
                              import scanpy
                              
                              # Load input file
                              infile = _________
                              adata = sc.read_10x_h5(infile)
                              adata.var_names_make_unique()
                              adata.obs_names_make_unique()
                              
                              adata.var['mt'] = data.var_names.str.startswith('MT-')
                              adata.var['rb'] = data.var_names.str.startswith(('RPS','RPL'))
                              
                              sc.pp.filter_cells(adata, min_counts={input$nCount_min})
                              sc.pp.filter_cells(adata, min_genes={input$nFeature_min})
                              
                              adata = adata[(adata.obs['pct_counts_mt'] < {input$mt_threshold}) & (adata.obs['pct_counts_rb'] < {input$rb_threshold})]
                              ") })
      writeLines(scriptTemplate(), file)
    }
  )
  
  output$saveReport <- downloadHandler(
    filename = function() { paste("Figure", Sys.Date(), '.png', sep="")},
    content = function(file){
      
      temp_html <- tempfile(fileext = ".html")
      plotly::plotly_save(ggplotly(plot_novelty()), file = temp_html)
      
      # Use webshot2 to convert HTML file to PNG
      webshot(temp_html, file = file, vwidth = 800, vheight = 600)
    },
    contentType = "image/png"
  )
  
}

options <- list()
if (!interactive()) {
  options$launch.browser <- FALSE
  options$host <- "0.0.0.0"
  options$port <- 3838
}
shinyApp(ui = ui, server = server, options = options)



# nCount_RNA - number of UMIs per cell
# nFeature_RNA - number of genes detected per cell
