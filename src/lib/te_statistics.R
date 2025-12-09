#' Comprehensive TE Statistics Analysis
#' This script provides modular functions for analyzing transposable element statistics
#' that can be easily integrated into R packages or used standalone.
#' @author Robert Schwarz
#' @date 2025-09-26
#' 
library(dplyr)
library(ggplot2)
library(gt)
library(kableExtra)
library(scales)
library(treemap)
library(DT)
library(plotly)
library(RColorBrewer)

example_data <- read.csv(file = './foo/te_stats_example.csv')
example_data <- example_data %>% filter(!duplicated(te_id), order %in% c("LINE", "SINE", "LTR", "DNA"))


example_data_zwei <- data.frame(
    ensembl_gene_id = paste0("ENSG0000", rep(1:10, 1)),
    order = c(rep("SINE", 2), rep("LINE", 4), rep("LTR",3), "DNA"),
    super_family = c("B1", "B1", "L1", "L1", "L2", "L2", "ERVK", "ERVK_Mal", "ERVK", "DNA_1"),
    family = c("B1_fam1", "B1_fam2", "L1_3end", "L1_3end", "L2_4", "L2_1", "ERVK_n", "ERVK_Mal_m", "ERVK_n", "DNA_1_2")
)


#' Generate TE Color Palette
#'
#' Creates a standardized color palette for TE classes
#'
#' @param te_classes Character vector of TE class names
#' @param palette_name Name of RColorBrewer palette (default: "Set2")
#' @return Named vector of colors
#' @export
generate_te_colors <- function(te_classes, palette_name = "Set2") {
    
    # Standard TE classes with default colors
    standard_colors <- c(
        "LTR" = "#1b9e77",
        "LINE" = "#7570b3", 
        "SINE" = "#d95f02",
        "DNA" = "#e7298a",
        "RC" = "#66a61e",
        "Unknown" = "#999999",
        "Other" = "#cccccc"
    )
    
    # Use standard colors if available, otherwise generate from palette
    if (all(te_classes %in% names(standard_colors))) {
        return(standard_colors[te_classes])
    } else {
        n_colors <- length(unique(te_classes))
        if (n_colors <= 8) {
            palette_colors <- RColorBrewer::brewer.pal(max(3, n_colors), palette_name)
        } else {
            palette_colors <- rainbow(n_colors)
        }
        names(palette_colors) <- unique(te_classes)
        return(palette_colors)
    }
}




#' Calculate TE Hierarchy Statistics
#'
#' Computes counts at Family, Superfamily, and Class levels
#'
#' @param te_df Data frame with columns: te_class, te_superfamily, te_family, ensembl_gene_id
#' @param entity_col Name of the column containing entities to count (default: "ensembl_gene_id")
#' @return List with family_stats, superfamily_stats, class_stats, and combined_stats
#' @export
calculate_te_hierarchy <- function(te_df, entity_col = "ensembl_gene_id") {
    
    # Validate input
    required_cols <- c("order", "super_family", "family", entity_col)
    missing_cols <- setdiff(required_cols, colnames(te_df))
    if (length(missing_cols) > 0) {
        stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
    }
    
    # Count unique entities at each hierarchical level
    family_stats <- te_df %>%
        group_by(order, super_family, family) %>%
        summarise(n_family = n_distinct(.data[[entity_col]]), .groups = "drop")
    
    superfamily_stats <- te_df %>%
        group_by(order, super_family) %>%
        summarise(n_superfamily = n_distinct(.data[[entity_col]]), .groups = "drop")
    
    class_stats <- te_df %>%
        group_by(order) %>%
        summarise(n_order = n_distinct(.data[[entity_col]]), .groups = "drop")
    
    # Combine all levels
    combined_stats <- family_stats %>%
        left_join(superfamily_stats, by = c("order", "super_family")) %>%
        left_join(class_stats, by = "order") %>%
        mutate(
            pct_family_in_superfamily = round(100 * n_family / n_superfamily, 2),
            pct_family_in_order = round(100 * n_family / n_order, 2),
            pct_superfamily_in_class = round(100 * n_superfamily / n_order, 2)
        ) %>%
        arrange(order, super_family, desc(n_family))
    
    return(list(
        family_stats = family_stats,
        superfamily_stats = superfamily_stats,
        class_stats = class_stats,
        combined_stats = combined_stats
    ))
}

#' Create Publication-Ready TE Table
#'
#' Generates formatted tables suitable for publications
#'
#' @param combined_stats Output from calculate_te_hierarchy()
#' @param top_n Number of top families to show per class (default: 10)
#' @param format Output format: "gt" or "kable" (default: "gt")
#' @param colors Named vector of colors for TE classes
#' @return Formatted table object
#' @export
create_te_table <- function(combined_stats, top_n = 10, format = "gt", colors = NULL) {
    
    # Select top N families per class
    top_families <- combined_stats %>%
        group_by(order) %>%
        slice_max(order_by = n_family, n = top_n) %>%
        ungroup()
    
    if (is.null(colors)) {
        colors <- generate_te_colors(unique(top_families$order))
    }
    
    if (format == "gt") {
        # GT table with grouping and colors
        color_fn <- scales::col_factor(palette = colors, 
                                       domain = names(colors), 
                                       na.color = "#cccccc")
        
        gt_tbl <- top_families %>%
            gt(groupname_col = "order") %>%
            cols_label(
                super_family = "Superfamily",
                family = "Family",
                n_family = "Count",
                pct_family_in_order = "% of Class"
            ) %>%
            cols_hide(columns = c(n_superfamily, n_order, pct_family_in_superfamily, pct_superfamily_in_class)) %>%
            data_color(columns = order, fn = color_fn) %>%
            fmt_number(columns = n_family, decimals = 0) %>%
            fmt_number(columns = pct_family_in_order, decimals = 1) %>%
            tab_header(title = paste("Top", top_n, "TE Families per Class")) %>%
            opt_stylize(style = 6, color = "blue")
        
        return(gt_tbl)
        
    } else if (format == "kable") {
        # KableExtra version
        top_families_kbl <- top_families %>%
            mutate(te_class_col = cell_spec(order, format = "html",
                                            background = colors[order],
                                            color = "white")) %>%
            dplyr::select(te_class_col, super_family, family, n_family, pct_family_in_order)
        
        kbl_tbl <- top_families_kbl %>%
            kbl(col.names = c("TE Class", "Superfamily", "Family", "Count", "% of Class"),
                escape = FALSE, booktabs = TRUE) %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
            pack_rows(index = table(top_families$order))
        
        return(kbl_tbl)
    }
}

#' Create TE Visualization Plots
#'
#' Generates various plots for TE statistics
#'
#'TODO
#' - add super_fam_barplot
#' 
#' @param stats_list Output from calculate_te_hierarchy()
#' @param colors Named vector of colors for TE classes
#' @param plot_types Vector of plot types to generate: "barplot", "treemap", "heatmap", "pie"
#' @return List of ggplot objects
#' @export
create_te_plots <- function(stats_list, colors = NULL, 
                            plot_types = c("barplot", "treemap", "heatmap")) {
    
    combined_stats <- stats_list$combined_stats
    class_stats <- stats_list$class_stats
    
    if (is.null(colors)) {
        colors <- generate_te_colors(unique(combined_stats$order))
    }
    
    plot_list <- list()
    
    # 1. Class overview barplot
    if ("barplot" %in% plot_types) {
        p_class <- ggplot(class_stats, aes(x = reorder(order, n_order), 
                                           y = n_order, fill = order)) +
            geom_col() +
            scale_fill_manual(values = colors) +
            coord_flip() +
            labs(title = "TE Count by Class",
                 x = "TE Class", y = "Count", fill = "TE Class") +
            theme_minimal() +
            theme(legend.position = "none")
        
        plot_list$class_barplot <- p_class
        
        # Top families barplot
        top_families <- combined_stats %>%
            group_by(order) %>%
            slice_max(order_by = n_family, n = 5) %>%
            ungroup()
        
        p_families <- ggplot(top_families, aes(x = reorder(family, n_family), 
                                               y = n_family, fill = order)) +
            geom_col() +
            scale_fill_manual(values = colors) +
            facet_wrap(~order, scales = "free") +
            coord_flip() +
            labs(title = "Top 5 TE Families per Class",
                 x = "TE Family", y = "Count", fill = "TE Class") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")
        
        plot_list$family_barplot <- p_families
    }
    
    # 2. Treemap
    if ("treemap" %in% plot_types && requireNamespace("treemap", quietly = TRUE)) {
        # Prepare data for treemap
        treemap_data <- combined_stats %>%
            dplyr::select(order, super_family, family, n_family) %>%
            arrange(desc(n_family))
        
        # Create treemap
        treemap_plot <- treemap::treemap(
            treemap_data,
            index = c("order", "super_family", "family"),
            vSize = "n_family",
            palette = colors,
            title = "TE Hierarchy Treemap"
        )
        
        plot_list$treemap <- treemap_plot
    }
    
    # 3. Heatmap of superfamily vs class
    if ("heatmap" %in% plot_types) {
        superfamily_stats <- stats_list$superfamily_stats
        
        p_heatmap <- ggplot(superfamily_stats, 
                            aes(x = order, y = super_family, fill = n_superfamily)) +
            geom_tile() +
            scale_fill_gradient(low = "darkblue", high = "darkred") +
            labs(title = "TE Superfamily Distribution Heatmap",
                 x = "TE Class", y = "TE Superfamily", fill = "Count") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        plot_list$heatmap <- p_heatmap
    }
    
    # 4. Pie chart for classes
    if ("pie" %in% plot_types) {
        p_pie <- ggplot(class_stats, aes(x = "", y = n_order, fill = order)) +
            geom_col() +
            scale_fill_manual(values = colors) +
            coord_polar("y", start = 0) +
            labs(title = "TE Class Distribution", fill = "TE Class") +
            theme_void()
        
        plot_list$pie_chart <- p_pie
    }
    
    return(plot_list)
}


#' Create Interactive TE Table
#'
#' Generates interactive DT table for web/HTML output
#'
#' @param combined_stats Output from calculate_te_hierarchy()
#' @param colors Named vector of colors for TE classes
#' @return DT datatable object
#' @export
create_interactive_table <- function(combined_stats, colors = NULL) {
    
    if (is.null(colors)) {
        colors <- generate_te_colors(unique(combined_stats$order))
    }
    
    # Format for display
    display_data <- combined_stats %>%
        dplyr::select(order, super_family, family, n_family, 
               pct_family_in_order, n_superfamily, n_order) %>%
        dplyr::rename(
            "TE Class" = order,
            "Superfamily" = super_family,
            "Family" = family,
            "Family Count" = n_family,
            "% of Class" = pct_family_in_order,
            "Superfamily Count" = n_superfamily,
            "Class Count" = n_order
        )
    
    dt_table <- DT::datatable(
        display_data,
        options = list(
            pageLength = 25,
            scrollX = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel')
        ),
        rownames = FALSE,
        caption = "Interactive TE Statistics Table"
    ) %>%
        DT::formatStyle(
            "TE Class",
            backgroundColor = DT::styleEqual(names(colors), colors)
        ) %>%
        DT::formatRound(c("% of Class"), 2)
    
    return(dt_table)
}

#' Main TE Analysis Function
#'
#' Comprehensive analysis of TE statistics with multiple outputs
#'
#' @param te_df Data frame with TE information
#' @param entity_col Column name containing entities to count (default: "entity_id")
#' @param top_n Number of top families per class for tables (default: 10)
#' @param colors Named vector of colors for TE classes (optional)
#' @param output_formats Vector of table formats: "gt", "kable", "interactive"
#' @param plot_types Vector of plot types to generate
#' @param export_path Path to export results (optional)
#' @return List containing all analysis results
#' @export
analyze_te_statistics <- function(te_df, 
                                  entity_col = "ensembl_gene_id",
                                  top_n = 10,
                                  colors = NULL,
                                  output_formats = c("gt", "interactive"),
                                  plot_types = c("barplot", "treemap", "heatmap"),
                                  export_path = NULL) {
    
    cat("Starting comprehensive TE analysis...\n")
    
    # 1. Calculate hierarchical statistics
    cat("Calculating hierarchical statistics...\n")
    stats <- calculate_te_hierarchy(te_df, entity_col)
    
    # 2. Generate colors if not provided
    if (is.null(colors)) {
        colors <- generate_te_colors(unique(stats$combined_stats$order))
    }
    
    # 3. Create tables
    cat("Creating formatted tables...\n")
    tables <- list()
    if ("gt" %in% output_formats) {
        tables$gt_table <- create_te_table(stats$combined_stats, top_n, "gt", colors)
    }
    if ("kable" %in% output_formats) {
        tables$kable_table <- create_te_table(stats$combined_stats, top_n, "kable", colors)
    }
    if ("interactive" %in% output_formats) {
        tables$interactive_table <- create_interactive_table(stats$combined_stats, colors)
    }
    
    # 4. Create plots
    cat("Generating visualization plots...\n")
    plots <- create_te_plots(stats, colors, plot_types)
    
    # 5. Summary statistics
    summary_stats <- list(
        total_entities = length(unique(te_df[[entity_col]])),
        total_families = nrow(stats$family_stats),
        total_superfamilies = nrow(stats$superfamily_stats),
        total_classes = nrow(stats$class_stats),
        top_class = stats$class_stats$order[which.max(stats$class_stats$n_order)],
        top_family = stats$combined_stats$family[which.max(stats$combined_stats$n_family)]
    )
    
    # 6. Export results if path provided
    if (!is.null(export_path)) {
        cat("Exporting results to:", export_path, "\n")
        dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
        
        # Export data
        write.csv(stats$combined_stats, file.path(export_path, "te_combined_stats.csv"), row.names = FALSE)
        write.csv(stats$class_stats, file.path(export_path, "te_class_stats.csv"), row.names = FALSE)
        
        # Export plots
        if ("barplot" %in% names(plots)) {
            ggsave(file.path(export_path, "te_class_barplot.png"), plots$class_barplot, width = 10, height = 6)
        }
        if ("family_barplot" %in% names(plots)) {
            ggsave(file.path(export_path, "te_family_barplot.png"), plots$family_barplot, width = 12, height = 8)
        }
    }
    
    # 7. Compile results
    results <- list(
        statistics = stats,
        tables = tables,
        plots = plots,
        summary = summary_stats,
        colors = colors
    )
    
    cat("Analysis complete!\n")
    cat("Summary:\n")
    cat("- Total entities:", summary_stats$total_entities, "\n")
    cat("- Total TE classes:", summary_stats$total_classes, "\n")
    cat("- Total TE families:", summary_stats$total_families, "\n")
    cat("- Most abundant class:", summary_stats$top_class, "\n")
    cat("- Most abundant family:", summary_stats$top_family, "\n")
    
    return(results)
}

# Example usage function
#' Example TE Analysis
#'
#' Demonstrates usage with sample data
#'
#' @export
example_te_analysis <- function() {
    # Create sample TE data
    set.seed(123)
    n_entities <- 1000
    
    sample_data <- data.frame(
        entity_id = paste0("entity_", 1:n_entities),
        te_class = sample(c("LTR", "LINE", "SINE", "DNA"), n_entities, 
                          replace = TRUE, prob = c(0.3, 0.35, 0.25, 0.1)),
        te_superfamily = sample(c("ERV1", "ERV2", "L1", "L2", "Alu", "MIR", "hAT", "Tc1"), 
                                n_entities, replace = TRUE),
        te_family = sample(paste0("Family_", 1:20), n_entities, replace = TRUE),
        stringsAsFactors = FALSE
    )
    
    # Run analysis
    results <- analyze_te_statistics(
        sample_data,
        entity_col = "entity_id",
        top_n = 5,
        output_formats = c("gt", "interactive"),
        plot_types = c("barplot", "heatmap", "pie")
    )
    
    return(results)
}

cat("TE Statistics Analysis Functions Loaded Successfully!\n")
cat("Main function: analyze_te_statistics()\n")
cat("Example usage: example_te_analysis()\n")

