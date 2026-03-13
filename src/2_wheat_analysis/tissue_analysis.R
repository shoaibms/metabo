# temporal_analysis_core.R
library(nparLD)
library(effsize)
library(boot)
library(data.table)
library(logger)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set file paths
data_path <- "C:/Users/USER/Desktop/r/chem_data/metabo2/Merged_VIP_hub_r_Path2.csv"
out_dir <- "C:/Users/USER/Desktop/r/chem_data/metabo2/result/section3"

# Initialize logging
log_dir <- file.path(out_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_appender(appender_file(file.path(log_dir, "temporal_analysis.log")))

calculate_tss <- function(data) {
  log_info("Starting TSS calculation")
  
  # Log data dimensions
  log_info(sprintf("Number of metabolites: %d", length(unique(data$Metabolite))))
  log_info(sprintf("Number of genotypes: %d", length(unique(data$Genotype))))
  
  # Split data by tissue
  root_data <- data[Tissue.type == "R"]
  leaf_data <- data[Tissue.type == "L"]
  
  log_info(sprintf("Root data rows: %d, Leaf data rows: %d", 
                   nrow(root_data), nrow(leaf_data)))
  
  results <- data.table()
  
  for(metabolite in unique(data$Metabolite)) {
    for(genotype in unique(data$Genotype)) {
      tryCatch({
        # Get metabolite data for each tissue
        root_met <- root_data[Metabolite == metabolite & Genotype == genotype]
        leaf_met <- leaf_data[Metabolite == metabolite & Genotype == genotype]
        
        if(nrow(root_met) > 0 && nrow(leaf_met) > 0) {
          # Response Timing Score
          root_peak <- root_met[which.max(abs(Metabolite_Value)), Day]
          leaf_peak <- leaf_met[which.max(abs(Metabolite_Value)), Day]
          timing_diff <- abs(root_peak - leaf_peak)
          rts <- 1 / (1 + timing_diff)
          
          # Basic Response Rate Score - using original calculation for stability
          root_rates <- diff(root_met[order(Day)]$Metabolite_Value)
          leaf_rates <- diff(leaf_met[order(Day)]$Metabolite_Value)
          rrs <- ifelse(length(root_rates) > 0 && length(leaf_rates) > 0,
                        (cor(root_rates, leaf_rates, method="spearman", 
                             use="complete.obs") + 1) / 2,
                        NA)
          
          # Basic Magnitude Alignment Score
          mas <- (cor(root_met$Metabolite_Value, leaf_met$Metabolite_Value, 
                      method="spearman", use="complete.obs") + 1) / 2
          
          # Calculate total TSS
          tss <- mean(c(rts, rrs, mas), na.rm = TRUE)
          
          if(!is.na(tss)) {  # Only add if TSS is valid
            results <- rbind(results, data.table(
              Metabolite = metabolite,
              Genotype = genotype,
              TSS = tss,
              RTS = rts,
              RRS = rrs,
              MAS = mas
            ))
          }
        }
      }, error = function(e) {
        log_error(sprintf("Error processing metabolite %s, genotype %s: %s", 
                          metabolite, genotype, e$message))
      })
    }
  }
  
  log_info(sprintf("TSS calculation completed. Results for %d combinations", 
                   nrow(results)))
  
  if(nrow(results) == 0) {
    log_error("No valid results generated!")
    return(NULL)
  }
  
  return(results)
}


#' Perform FDR correction on temporal data
#' @param data Data frame containing metabolomics data
#' @return List of FDR-corrected results
perform_fdr_correction <- function(data) {
  tryCatch({
    log_info("Starting FDR correction")
    results <- data[, {
      # Calculate p-values for temporal changes using Kruskal-Wallis
      pvals <- kruskal.test(Metabolite_Value ~ Day)$p.value
      list(p_value = pvals)
    }, by = .(Metabolite, Tissue.type, Genotype)]
    
    # Apply BH correction
    results[, fdr := p.adjust(p_value, method = "BH")]
    log_info("FDR correction completed")
    return(results)
  }, error = function(e) {
    log_error("FDR correction failed: {e$message}")
    return(NULL)
  })
}

# Replace the perform_temporal_analysis function with this version
perform_temporal_analysis <- function(data) {
  log_info("Starting robust non-parametric analysis")
  
  # Convert to data.table
  DT <- as.data.table(data)
  
  # Prepare data structure
  DT[, `:=`(
    Day = as.numeric(as.character(Day)),
    Genotype = factor(Genotype),
    Tissue.type = factor(Tissue.type),
    Subject = factor(paste(Genotype, Replication, Batch))
  )]
  
  # Analyze by group using data.table
  results <- DT[, {
    # Friedman test for time effect
    time_test <- tryCatch({
      friedman.test(
        y = Metabolite_Value,
        groups = Day,
        blocks = Subject
      )
    }, error = function(e) {
      list(statistic = NA_real_, p.value = NA_real_)
    })
    
    # Kruskal-Wallis for genotype effect
    genotype_test <- tryCatch({
      kruskal.test(
        Metabolite_Value ~ Genotype
      )
    }, error = function(e) {
      list(statistic = NA_real_, p.value = NA_real_)
    })
    
    # Interaction analysis
    interaction_p <- tryCatch({
      unique_days <- unique(Day)
      p_values <- numeric(length(unique_days))
      
      for(i in seq_along(unique_days)) {
        d <- unique_days[i]
        time_data <- .SD[Day == d]
        w <- wilcox.test(
          Metabolite_Value ~ Genotype,
          data = time_data,
          paired = FALSE
        )
        p_values[i] <- w$p.value
      }
      min(p_values)
    }, error = function(e) {
      NA_real_
    })
    
    # Effect sizes
    time_effect <- tryCatch({
      cliff.delta(Metabolite_Value, Day)$estimate
    }, error = function(e) {
      NA_real_
    })
    
    genotype_effect <- tryCatch({
      cliff.delta(Metabolite_Value, Genotype)$estimate
    }, error = function(e) {
      NA_real_
    })
    
    # Return results as a list
    list(
      time_stat = as.numeric(time_test$statistic),
      time_p = time_test$p.value,
      genotype_stat = as.numeric(genotype_test$statistic),
      genotype_p = genotype_test$p.value,
      interaction_p = interaction_p,
      interaction_stat = interaction_p * length(unique(Day)),
      time_effect = time_effect,
      genotype_effect = genotype_effect
    )
  }, by = .(Metabolite, Tissue.type)]
  
  # Add FDR corrections
  results[, `:=`(
    time_p_adj = p.adjust(time_p, method = "BH"),
    genotype_p_adj = p.adjust(genotype_p, method = "BH"),
    interaction_p_adj = p.adjust(interaction_p, method = "BH")
  )]
  
  log_info("Temporal analysis completed")
  return(results)
}

#' Calculate effect sizes with bootstrap CIs
#' @param data Data frame containing metabolomics data
#' @return Data table with effect sizes and CIs
calculate_effect_sizes <- function(data) {
  tryCatch({
    log_info("Starting effect size calculation")
    results <- data[, {
      first_day <- min(Day)
      last_day <- max(Day)
      
      first_values <- Metabolite_Value[Day == first_day]
      last_values <- Metabolite_Value[Day == last_day]
      
      cd <- cliff.delta(first_values, last_values)
      
      boot_data <- data.frame(
        value = c(first_values, last_values),
        group = rep(c("first", "last"), c(length(first_values), length(last_values)))
      )
      
      boot_fn <- function(d, i) {
        subset_data <- d[i,]
        cliff.delta(subset_data$value ~ subset_data$group)$estimate
      }
      
      boot_res <- boot(boot_data, boot_fn, R = 5000)
      
      list(
        effect_size = cd$estimate,
        ci_lower = ci$percent[4],
        ci_upper = ci$percent[5],
        magnitude = cd$magnitude
      )
    }, by = .(Metabolite, Tissue.type, Genotype)]
    
    log_info("Effect size calculation completed")
    return(results)
  }, error = function(e) {
    log_error("Effect size calculation failed: {e$message}")
    return(NULL)
  })
}

#' Calculate temporal autocorrelation
#' @param data Data frame containing metabolomics data
#' @return Data table with autocorrelation results
calculate_temporal_correlation <- function(data) {
  tryCatch({
    log_info("Starting temporal correlation analysis")
    results <- data[, {
      spearman <- cor.test(as.numeric(Day), Metabolite_Value, 
                           method = "spearman", exact = FALSE)
      kendall <- cor.test(as.numeric(Day), Metabolite_Value, 
                          method = "kendall", exact = FALSE)
      
      list(
        spearman_rho = spearman$estimate,
        spearman_p = spearman$p.value,
        kendall_tau = kendall$estimate,
        kendall_p = kendall$p.value
      )
    }, by = .(Metabolite, Tissue.type, Genotype)]
    
    log_info("Temporal correlation analysis completed")
    return(results)
  }, error = function(e) {
    log_error("Temporal correlation calculation failed: {e$message}")
    return(NULL)
  })
}

#' Create visualization for TSS results
#' @param results TSS results data table
#' @param out_dir Output directory for plots


# Add these exports to your prepare_tss_data function
prepare_tss_data <- function(results) {
  tryCatch({
    log_info("Preparing TSS data for plotting")
    
    # Original TSS distribution data
    tss_distribution <- results$tss[, .(
      Metabolite,
      Genotype,
      TSS
    )]
    fwrite(tss_distribution, file.path(out_dir, "tss_distribution_data.csv"))
    
    # Component scores data
    component_scores <- melt(results$tss,
                             id.vars = c("Metabolite", "Genotype"),
                             measure.vars = c("RTS", "RRS", "MAS"),
                             variable.name = "Component",
                             value.name = "Score")
    
    fwrite(component_scores, file.path(out_dir, "tss_component_scores_data.csv"))
    
    # Basic summary statistics
    summary_stats <- results$tss[, .(
      tss_mean = mean(TSS, na.rm = TRUE),
      tss_median = median(TSS, na.rm = TRUE),
      tss_sd = sd(TSS, na.rm = TRUE),
      tss_iqr = IQR(TSS, na.rm = TRUE),
      rts_mean = mean(RTS, na.rm = TRUE),
      rts_median = median(RTS, na.rm = TRUE),
      rts_sd = sd(RTS, na.rm = TRUE),
      rts_iqr = IQR(RTS, na.rm = TRUE),
      rrs_mean = mean(RRS, na.rm = TRUE),
      rrs_median = median(RRS, na.rm = TRUE),
      rrs_sd = sd(RRS, na.rm = TRUE),
      rrs_iqr = IQR(RRS, na.rm = TRUE),
      mas_mean = mean(MAS, na.rm = TRUE),
      mas_median = median(MAS, na.rm = TRUE),
      mas_sd = sd(MAS, na.rm = TRUE),
      mas_iqr = IQR(MAS, na.rm = TRUE),
      n_metabolites = .N
    ), by = Genotype]
    
    fwrite(summary_stats, file.path(out_dir, "tss_summary_stats.csv"))
    
    # Read original data for additional summaries
    data <- fread(data_path)
    
    # Temporal patterns summary
    temporal_patterns <- data[, .(
      mean_value = mean(Metabolite_Value, na.rm = TRUE),
      sd_value = sd(Metabolite_Value, na.rm = TRUE),
      median_value = median(Metabolite_Value, na.rm = TRUE),
      iqr_value = IQR(Metabolite_Value, na.rm = TRUE),
      n_samples = .N
    ), by = .(Metabolite, Day)]
    fwrite(temporal_patterns, file.path(out_dir, "temporal_patterns.csv"))
    
    # Tissue summaries
    tissue_summaries <- data[, .(
      mean_value = mean(Metabolite_Value, na.rm = TRUE),
      sd_value = sd(Metabolite_Value, na.rm = TRUE),
      median_value = median(Metabolite_Value, na.rm = TRUE),
      iqr_value = IQR(Metabolite_Value, na.rm = TRUE),
      n_samples = .N
    ), by = .(Tissue.type, Metabolite)]
    fwrite(tissue_summaries, file.path(out_dir, "tissue_summaries.csv"))
    
    # Genotype summaries
    genotype_summaries <- data[, .(
      mean_value = mean(Metabolite_Value, na.rm = TRUE),
      sd_value = sd(Metabolite_Value, na.rm = TRUE),
      median_value = median(Metabolite_Value, na.rm = TRUE),
      iqr_value = IQR(Metabolite_Value, na.rm = TRUE),
      n_samples = .N
    ), by = .(Genotype, Metabolite)]
    fwrite(genotype_summaries, file.path(out_dir, "genotype_summaries.csv"))
    
    # Data validation summary
    validation_summary <- data[, .(
      n_samples = .N,
      missing_pct = sum(is.na(Metabolite_Value))/.N * 100,
      unique_timepoints = uniqueN(Day),
      min_value = min(Metabolite_Value, na.rm = TRUE),
      max_value = max(Metabolite_Value, na.rm = TRUE),
      mean_value = mean(Metabolite_Value, na.rm = TRUE),
      median_value = median(Metabolite_Value, na.rm = TRUE)
    ), by = .(Metabolite)]
    fwrite(validation_summary, file.path(out_dir, "data_validation.csv"))
    
    log_info("TSS data preparation completed with additional summaries")
  }, error = function(e) {
    log_error("Data preparation failed: {e$message}")
  })
}


# Enhanced main function with rigorous statistical validation
main <- function() {
  log_info("Starting temporal analysis with TSS")
  
  data <- fread(data_path)
  log_info(sprintf("Read %d rows of data", nrow(data)))
  
  # Validate input data
  required_cols <- c("Metabolite", "Tissue.type", "Genotype", "Day", 
                     "Metabolite_Value", "Replication", "Batch")
  if (!all(required_cols %in% names(data))) {
    stop("Required columns missing from input data")
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- list()
  
  # Original analyses - keep these unchanged
  results$tss <- calculate_tss(data)
  results$fdr <- perform_fdr_correction(data)
  results$temporal_analysis <- perform_temporal_analysis(data)
  results$effect <- calculate_effect_sizes(data)
  results$corr <- calculate_temporal_correlation(data)
  
  # Additional tissue comparison metrics with statistical validation
  results$tissue_comparison <- data[, {
    tryCatch({
      # Bootstrap for response magnitude
      boot_magnitude <- boot(data = .SD$Metabolite_Value,
                             statistic = function(d, i) {
                               sd(d[i]) / mean(abs(d[i]))
                             }, R = 20)
      magnitude_ci <- boot.ci(boot_magnitude, type = "perc")
      
      # Temporal trend with permutation test
      obs_trend <- lm(Metabolite_Value ~ as.numeric(Day))
      obs_slope <- coef(obs_trend)[2]
      
      n_perm <- 10000
      perm_slopes <- replicate(n_perm, {
        perm_data <- .SD$Metabolite_Value[sample(.N)]
        coef(lm(perm_data ~ as.numeric(Day)))[2]
      })
      trend_p <- mean(abs(perm_slopes) >= abs(obs_slope))
      
      list(
        response_magnitude = boot_magnitude$t0,
        magnitude_ci_lower = magnitude_ci$percent[4],
        magnitude_ci_upper = magnitude_ci$percent[5],
        trend_slope = obs_slope,
        trend_p = trend_p,
        trend_r2 = summary(obs_trend)$r.squared
      )
    }, error = function(e) {
      list(
        response_magnitude = NA,
        magnitude_ci_lower = NA,
        magnitude_ci_upper = NA,
        trend_slope = NA,
        trend_p = NA,
        trend_r2 = NA
      )
    })
  }, by = .(Tissue.type, Genotype, Metabolite)]
  
  # Pathway-level analysis with statistical validation
  results$pathway_summary <- data[, {
    tryCatch({
      # Bootstrap for pathway metrics
      boot_response <- boot(data = .SD$Metabolite_Value,
                            statistic = function(d, i) mean(d[i]),
                            R = 20)
      response_ci <- boot.ci(boot_response, type = "perc")
      
      boot_var <- boot(data = .SD$Metabolite_Value,
                       statistic = function(d, i) var(d[i]),
                       R = 20)
      var_ci <- boot.ci(boot_var, type = "perc")
      
      # Temporal range with permutation test
      obs_range <- diff(range(.SD$Metabolite_Value))
      perm_ranges <- replicate(20, {
        diff(range(sample(.SD$Metabolite_Value)))
      })
      range_p <- mean(perm_ranges >= obs_range)
      
      list(
        mean_response = boot_response$t0,
        response_ci_lower = response_ci$percent[4],
        response_ci_upper = response_ci$percent[5],
        response_var = boot_var$t0,
        var_ci_lower = var_ci$percent[4],
        var_ci_upper = var_ci$percent[5],
        temporal_range = obs_range,
        range_p = range_p
      )
    }, error = function(e) {
      list(
        mean_response = NA,
        response_ci_lower = NA,
        response_ci_upper = NA,
        response_var = NA,
        var_ci_lower = NA,
        var_ci_upper = NA,
        temporal_range = NA,
        range_p = NA
      )
    })
  }, by = .(Tissue.type, Genotype, Day)]
  
  # Apply FDR correction to new p-values
  results$tissue_comparison[, trend_p_adj := p.adjust(trend_p, method = "BH")]
  results$pathway_summary[, range_p_adj := p.adjust(range_p, method = "BH")]
  
  return(results)
}

export_results <- function(results) {
  tryCatch({
    log_info("Starting data export")
    
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Original core results
    if (!is.null(results$tss)) {
      fwrite(results$tss, file.path(out_dir, "tss_results.csv"))
      
      # TSS distribution data
      tss_distribution <- results$tss[, .(
        Metabolite,
        Genotype,
        TSS
      )]
      fwrite(tss_distribution, file.path(out_dir, "tss_distribution_data.csv"))
      
      # Component scores data
      component_scores <- melt(results$tss,
                               id.vars = c("Metabolite", "Genotype"),
                               measure.vars = c("RTS", "RRS", "MAS"),
                               variable.name = "Component",
                               value.name = "Score")
      fwrite(component_scores, file.path(out_dir, "tss_component_scores_data.csv"))
      
      # Summary statistics
      summary_stats <- results$tss[, .(
        tss_mean = mean(TSS, na.rm = TRUE),
        tss_median = median(TSS, na.rm = TRUE),
        tss_sd = sd(TSS, na.rm = TRUE),
        tss_iqr = IQR(TSS, na.rm = TRUE),
        n_metabolites = .N
      ), by = Genotype]
      fwrite(summary_stats, file.path(out_dir, "tss_summary_stats.csv"))
    }
    
    # Other original results
    if (!is.null(results$temporal_analysis)) {
      fwrite(results$temporal_analysis, file.path(out_dir, "temporal_analysis.csv"))
    }
    if (!is.null(results$fdr)) {
      fwrite(results$fdr, file.path(out_dir, "fdr_results.csv"))
    }
    if (!is.null(results$effect)) {
      fwrite(results$effect, file.path(out_dir, "effect_sizes.csv"))
    }
    if (!is.null(results$corr)) {
      fwrite(results$corr, file.path(out_dir, "temporal_corr.csv"))
    }
    
    # Read original data for additional summaries
    data <- fread(data_path)
    
    # Temporal patterns summary
    temporal_patterns <- data[, .(
      mean_value = mean(Metabolite_Value, na.rm = TRUE),
      sd_value = sd(Metabolite_Value, na.rm = TRUE),
      median_value = median(Metabolite_Value, na.rm = TRUE),
      iqr_value = IQR(Metabolite_Value, na.rm = TRUE),
      n_samples = .N
    ), by = .(Metabolite, Day)]
    fwrite(temporal_patterns, file.path(out_dir, "temporal_patterns.csv"))
    
    # Tissue summaries
    tissue_summaries <- data[, .(
      mean_value = mean(Metabolite_Value, na.rm = TRUE),
      sd_value = sd(Metabolite_Value, na.rm = TRUE),
      median_value = median(Metabolite_Value, na.rm = TRUE),
      iqr_value = IQR(Metabolite_Value, na.rm = TRUE),
      n_samples = .N
    ), by = .(Tissue.type, Metabolite)]
    fwrite(tissue_summaries, file.path(out_dir, "tissue_summaries.csv"))
    
    # Genotype summaries
    genotype_summaries <- data[, .(
      mean_value = mean(Metabolite_Value, na.rm = TRUE),
      sd_value = sd(Metabolite_Value, na.rm = TRUE),
      median_value = median(Metabolite_Value, na.rm = TRUE),
      iqr_value = IQR(Metabolite_Value, na.rm = TRUE),
      n_samples = .N
    ), by = .(Genotype, Metabolite)]
    fwrite(genotype_summaries, file.path(out_dir, "genotype_summaries.csv"))
    
    # Data validation summary
    validation_summary <- data[, .(
      n_samples = .N,
      missing_pct = sum(is.na(Metabolite_Value))/.N * 100,
      unique_timepoints = uniqueN(Day),
      min_value = min(Metabolite_Value, na.rm = TRUE),
      max_value = max(Metabolite_Value, na.rm = TRUE),
      mean_value = mean(Metabolite_Value, na.rm = TRUE),
      median_value = median(Metabolite_Value, na.rm = TRUE)
    ), by = .(Metabolite)]
    fwrite(validation_summary, file.path(out_dir, "data_validation.csv"))
    
    # Analysis parameters
    params <- data.table(
      parameter = c("fdr_method", "bootstrap_iterations", "correlation_methods",
                    "effect_size_method", "confidence_level", "tss_components"),
      value = c("BH", "20", "spearman,kendall", "cliff.delta", "0.95", 
                "RTS,RRS,MAS")
    )
    fwrite(params, file.path(out_dir, "analysis_params.csv"))
    
    # New statistically validated tissue comparison results
    if (!is.null(results$tissue_comparison)) {
      # Full results with CIs and p-values
      fwrite(results$tissue_comparison, 
             file.path(out_dir, "tissue_comparison_full.csv"))
      
      # Summary with statistical significance
      tissue_summary <- results$tissue_comparison[, .(
        mean_magnitude = mean(response_magnitude, na.rm = TRUE),
        magnitude_ci_width = mean(magnitude_ci_upper - magnitude_ci_lower, na.rm = TRUE),
        significant_trends = sum(trend_p_adj < 0.05, na.rm = TRUE),
        mean_r2 = mean(trend_r2, na.rm = TRUE)
      ), by = .(Tissue.type, Genotype)]
      
      fwrite(tissue_summary, 
             file.path(out_dir, "tissue_comparison_summary.csv"))
    }
    
    # New pathway summary with statistical validation
    if (!is.null(results$pathway_summary)) {
      # Full results with CIs and p-values
      fwrite(results$pathway_summary, 
             file.path(out_dir, "pathway_temporal_full.csv"))
      
      # Statistical summary
      pathway_stats <- results$pathway_summary[, .(
        significant_ranges = sum(range_p_adj < 0.05, na.rm = TRUE),
        mean_response = mean(mean_response, na.rm = TRUE),
        response_ci_width = mean(response_ci_upper - response_ci_lower, na.rm = TRUE),
        mean_variance = mean(response_var, na.rm = TRUE),
        var_ci_width = mean(var_ci_upper - var_ci_lower, na.rm = TRUE)
      ), by = .(Tissue.type, Genotype)]
      
      fwrite(pathway_stats, 
             file.path(out_dir, "pathway_statistics.csv"))
    }
    
    log_info("Results exported successfully")
  }, error = function(e) {
    log_error(paste("Export failed:", e$message))
    stop(e)
  })
}

# Run the analysis
if (!dir.exists(dirname(out_dir))) {
  dir.create(dirname(out_dir), recursive = TRUE)
}

results <- main()
export_results(results)

