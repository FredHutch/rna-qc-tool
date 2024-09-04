# rna-qc-tool

A Shiny app for interactive visualization of sinlge-cell RNA quality control statistics.

## How it works?

Click on `Browse` button in the left pane and upload a (.h5) file. Specify the input parameters (or leave default values). After uploading is complete, click `Apply Filter` button. The application automatically annotates mitochondrial and ribosomal genes. Next, it filters all the cells and features based on parameters selected in the left pane, in real-time. There parameters can be adjusted anytime.

In a few seconds, 2 tabs on the right will be populated:

-   `Table`: lists all the commonly used statistics for filtered dataset

-   `Plots`: generates all commonly used figures to visualize filtered dataset

## Features

1.  The plots are rendered live (changing parameters will change both tabs live)

2.  The plots are interactive (hover mouse over plots to zoom, see specific values)

3.  `Download Script` button allows you to download scripts that are pre-populated with selected parameters for quality control in Python (using scanpy) and R (seurat)
