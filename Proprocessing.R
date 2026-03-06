# Function to compute descriptive stats for one variable

# Function to compute descriptive stats for one variable
describe_var <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  test_use <- NA
  # Default p-value
  p <- NA_real_
  
  # Handle too-small samples
  if (n < 3) {
    return(data.frame(
      n = n,
      mean = NA, median = NA, sd = NA,
      skewness = NA, kurtosis = NA, z_test = NA
    ))
  }
  
  # Choose appropriate normality test
  if (n >= 3 && n <= 5000) {
    p <- tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
    test_used <- "Shapiro–Wilk"
  } else if (n > 5000) {
    p <- tryCatch(nortest::ad.test(x)$p.value, error = function(e) NA_real_)
    test_used <- "Anderson–Darling"
  }
  
  # Safe descriptive stats
  data.frame(
    n = n,
    mean = tryCatch(mean(x), error = function(e) NA_real_),
    median = tryCatch(median(x), error = function(e) NA_real_),
    sd = tryCatch(sd(x), error = function(e) NA_real_),
    skewness = tryCatch(e1071::skewness(x, type = 2), error = function(e) NA_real_),
    kurtosis = tryCatch(e1071::kurtosis(x, type = 2), error = function(e) NA_real_),
    Normality_Test = test_used,
    p_value = as.numeric(p)
  )
}



# Function to create histogram for indicator
plot_hist <- function(data, indicators = names(data), bins = 30, ...) {
  plot_list <- list()
  for (ind in indicators) {
    if (is.numeric(data[[ind]])) {
      df <- data.frame(value = na.omit(data[[ind]]))
      
      p <- ggplot(df, aes(x = value)) +
        geom_histogram(aes(y = after_stat(density)),
                       bins = bins, fill = "skyblue", 
                       color = "white") + 
        geom_density(color = "darkblue", linewidth = 1) +
        labs(
          title = paste("Histogram of", ind),
          x = ind,
          y = "Density"
        ) +
        theme_minimal(base_size = 11)
    }
    plot_list <- append(plot_list, list(p))
  }
    gridExtra::grid.arrange(grobs = plot_list, ...)
}

# Function to create qq plot for indicators
plot_qq <- function(data, indicators = names(data), bins = 30, ...) {

  plot_list <- list()
  for (ind in indicators) {
    if (is.numeric(data[[ind]])) {
      df <- data.frame(value = na.omit(data[[ind]]))

      q <- ggplot(df, aes(sample = value)) +
        stat_qq(color = "darkorange") +
        stat_qq_line(color = "red", linetype = "dashed") +
        labs(title = paste(ind, "- Q–Q Plot"), x = "Theoretical", y = "Sample") +
        theme_minimal(base_size = 10)

      # Add to list
      plot_list <- append(plot_list, list(q))
    }
  }

  # Arrange all plots in one big grid (2 per indicator)
  nrow <- ceiling(length(plot_list) / ncol)

  gridExtra::grid.arrange(grobs = plot_list, ...)
}

# --- Helper function to plot domain-ordered heatmaps ---

plot_corr_heatmap <- function(df, area_name = "all") {

  if(area_name == "all"){
        df_area <- df %>%
    st_drop_geometry() %>%
    select(any_of(indicator_order)) %>%
    mutate(across(where(is.character), ~ as.numeric(str_remove(.x, "%")))) %>%
    drop_na()
  } else{
    df_area <- df %>%
    filter(Joint_Board_area == area_name) %>%
    st_drop_geometry() %>%
    select(any_of(indicator_order)) %>%
    mutate(across(where(is.character), ~ as.numeric(str_remove(.x, "%")))) %>%
    drop_na()
  }
  # Prepare data

  # Compute correlation
  cor_matrix <- cor(df_area, use = "pairwise.complete.obs")
  cor_melt <- melt(cor_matrix)

  # Enforce consistent indicator order
  cor_melt <- cor_melt %>%
    mutate(
      Var1 = factor(Var1, levels = indicator_order),
      Var2 = factor(Var2, levels = indicator_order)
    )

  # Plot
  ggplot(cor_melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation"
    ) +
    labs(
      title = paste("Correlation Heatmap of", area_name, "Indicators"),
      caption = "Source: Scottish Government (SIMD 2020v2)",
      x = "Indicator", y = "Indicator"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, size = 7, face = "bold"),
      axis.text.y = element_text(size = 7, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.caption = element_text(size = 7, hjust = 1)
    )
}

align_by_overlap <- function(ref, target) {
  # Convert to integers (drop unused levels)
  ref <- as.integer(as.factor(ref))
  target <- as.integer(as.factor(target))
  
  # Build contingency table
  tab <- table(ref, target)
  mat <- as.matrix(tab)
  
  # Make sure it's square for LSAP (pad with zeros if needed)
  n <- max(nrow(mat), ncol(mat))
  padded <- matrix(0, n, n)
  padded[1:nrow(mat), 1:ncol(mat)] <- mat
  
  # Solve the assignment problem (maximize total overlap)
  assignment <- solve_LSAP(padded, maximum = TRUE)
  
  # Construct mapping (target cluster → ref cluster)
  mapping <- seq_len(ncol(mat))
  mapping[seq_along(assignment)] <- assignment[seq_along(assignment)]
  
  # Apply mapping to target labels
  new_target <- sapply(target, function(x) {
    if (x <= length(mapping)) mapping[x] else NA
  })
  
  return(factor(new_target, levels = sort(unique(ref))))
}