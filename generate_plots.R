#!/usr/bin/env Rscript

suppressMessages({
    library(ggplot2)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    output_dir <- "."
} else {
    output_dir <- args[1]
}

plots_dir <- file.path(output_dir, "plots")

cat("Output directory:", output_dir, "\n")
cat("Plots directory:", plots_dir, "\n\n")

plot_enrichment <- function(summary_file, output_file, title, n_top = 20) {
    if (!file.exists(summary_file)) {
        cat("File not found:", summary_file, "\n")
        return(invisible(NULL))
    }
    
    tryCatch({
        data <- fread(summary_file)
        
        # Check if Competitive-P column exists
        if (!"Competitive-P" %in% colnames(data)) {
            cat("Warning: No Competitive-P column in", basename(summary_file), "\n")
            return(invisible(NULL))
        }
        
        # Sort and take top N
        data <- data[order(`Competitive-P`)][1:min(n_top, .N)]
        
        # Remove NA
        data <- data[!is.na(`Competitive-P`)]
        
        if (nrow(data) == 0) {
            cat("No valid data in", basename(summary_file), "\n")
            return(invisible(NULL))
        }
        
        # Create plot
        p <- ggplot(data, aes(x = reorder(Set, -log10(`Competitive-P`)), 
                              y = -log10(`Competitive-P`))) +
            geom_col(fill = "steelblue", alpha = 0.8) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                      color = "red", size = 1) +
            geom_hline(yintercept = -log10(0.01), linetype = "dashed", 
                      color = "darkred", size = 1) +
            coord_flip() +
            labs(title = title,
                 x = "Pathway / Gene Set",
                 y = "-log10(Competitive P-value)") +
            theme_minimal(base_size = 12) +
            theme(axis.text.y = element_text(size = 9),
                  plot.title = element_text(face = "bold", size = 14))
        
        ggsave(output_file, p, width = 12, height = 8, dpi = 300)
        cat("✅ Plot saved:", basename(output_file), "\n")
        
    }, error = function(e) {
        cat("Error plotting", basename(summary_file), ":", conditionMessage(e), "\n")
    })
}

# Generate plots for each analysis
cat("Generating enrichment plots...\n\n")

plot_enrichment(
    file.path(output_dir, "hallmark.summary"),
    file.path(plots_dir, "hallmark_enrichment_barplot.png"),
    "Hallmark Gene Sets - Sickle Cell Disease"
)

plot_enrichment(
    file.path(output_dir, "custom_sicklecell.summary"),
    file.path(plots_dir, "custom_sicklecell_enrichment_barplot.png"),
    "Custom Sickle Cell Gene Sets"
)

plot_enrichment(
    file.path(output_dir, "kegg.summary"),
    file.path(plots_dir, "kegg_enrichment_barplot.png"),
    "KEGG Pathways - Sickle Cell Disease"
)

plot_enrichment(
    file.path(output_dir, "reactome.summary"),
    file.path(plots_dir, "reactome_enrichment_barplot.png"),
    "Reactome Pathways - Sickle Cell Disease"
)

cat("\n✅ All plots generated successfully!\n")
