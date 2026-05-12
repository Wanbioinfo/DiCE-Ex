#' Create a DiCE centrality rule
#'
#' Helper function (not for users)
#'
#' @param metric Character string specifying the centrality metric.
#' @param threshold_type Character string specifying the cutoff type.
#'   Choose from \code{"percent"}, \code{"rank"}, or \code{"mean"}.
#' @param threshold Numeric cutoff value used when \code{threshold_type} is
#'   \code{"percent"} or \code{"rank"}. Ignored when
#'   \code{threshold_type = "mean"}.
#'
#' @return A list defining one centrality-based DiCE rule.
#' @noRd
dice_centrality_rule <- function(metric,
                                 threshold_type = c("percent", "rank", "mean"),
                                 threshold = NULL) {
  threshold_type <- match.arg(threshold_type)
  
  if (threshold_type %in% c("percent", "rank") && is.null(threshold)) {
    stop("threshold must be provided when threshold_type is 'percent' or 'rank'")
  }
  
  list(
    type = "centrality",
    metric = metric,
    threshold_type = threshold_type,
    threshold = threshold
  )
}

#' Create a DiCE ensemble rule
#'
#' Helper function (not for users)
#'
#' @param threshold_type Character string specifying the cutoff type.
#'   Choose from \code{"percent"} or \code{"rank"}.
#' @param threshold Numeric cutoff value for the ensemble rank.
#'
#' @return A list defining one ensemble-based DiCE rule.
#' @noRd
dice_ensemble_rule <- function(threshold_type = c("percent", "rank"),
                               threshold) {
  threshold_type <- match.arg(threshold_type)
  
  if (missing(threshold) || is.null(threshold)) {
    stop("threshold must be provided for ensemble rules")
  }
  
  list(
    type = "ensemble",
    threshold_type = threshold_type,
    threshold = threshold
  )
}

#' Resolve a numeric rank cutoff
#'
#' Helper function (not for users)
#'
#' @param n Total number of genes.
#' @param threshold_type Character string specifying the cutoff type.
#'   Choose from \code{"percent"} or \code{"rank"}.
#' @param threshold Numeric cutoff value.
#'
#' @return An integer rank cutoff.
#' @noRd
resolve_rank_cutoff <- function(n, threshold_type, threshold) {
  if (threshold_type == "percent") {
    k <- ceiling(n * as.numeric(threshold) / 100)
  } else if (threshold_type == "rank") {
    k <- as.integer(threshold)
  } else {
    stop("threshold_type must be 'percent' or 'rank'")
  }
  
  k <- max(1, min(k, n))
  return(k)
}

#' Evaluate a DiCE rule
#'
#' Helper function (not for users)
#'
#' @param df A data frame containing centrality values and/or
#'   \code{Ensemble_Rank}.
#' @param rule A DiCE rule created using \code{dice_centrality_rule()} or
#'   \code{dice_ensemble_rule()}.
#'
#' @return A logical vector indicating whether each gene passes the rule.
#' @noRd
evaluate_dice_rule <- function(df, rule) {
  
  n <- nrow(df)
  
  if (rule$type == "centrality") {
    pfx <- centrality_prefix(rule$metric)
    
    treat_col <- paste0(pfx, "_treatment")
    ctrl_col  <- paste0(pfx, "_control")
    
    if (!all(c(treat_col, ctrl_col) %in% colnames(df))) {
      stop("Could not find treatment/control columns for centrality: ", rule$metric)
    }
    
    if (rule$threshold_type == "mean") {
      treat_thr <- mean(df[[treat_col]], na.rm = TRUE)
      ctrl_thr  <- mean(df[[ctrl_col]], na.rm = TRUE)
      
      pass_vec <- (df[[treat_col]] >= treat_thr) | (df[[ctrl_col]] >= ctrl_thr)
      pass_vec[is.na(pass_vec)] <- FALSE
      return(pass_vec)
    }
    
    if (rule$threshold_type == "percent") {
      q <- 1 - (as.numeric(rule$threshold) / 100)
      
      treat_thr <- quantile(df[[treat_col]], probs = q, na.rm = TRUE,
                            names = FALSE, type = 7)
      ctrl_thr  <- quantile(df[[ctrl_col]], probs = q, na.rm = TRUE,
                            names = FALSE, type = 7)
      
      pass_vec <- (df[[treat_col]] >= treat_thr) | (df[[ctrl_col]] >= ctrl_thr)
      pass_vec[is.na(pass_vec)] <- FALSE
      return(pass_vec)
    }
    
    if (rule$threshold_type == "rank") {
      cutoff_rank <- resolve_rank_cutoff(
        n = n,
        threshold_type = "rank",
        threshold = rule$threshold
      )
      
      treat_rank <- rank(-df[[treat_col]], ties.method = "min", na.last = "keep")
      ctrl_rank  <- rank(-df[[ctrl_col]], ties.method = "min", na.last = "keep")
      
      pass_vec <- (treat_rank <= cutoff_rank) | (ctrl_rank <= cutoff_rank)
      pass_vec[is.na(pass_vec)] <- FALSE
      return(pass_vec)
    }
    
    stop("Unsupported threshold_type for centrality rule: ", rule$threshold_type)
  }
  
  if (rule$type == "ensemble") {
    if (!"Ensemble_Rank" %in% colnames(df)) {
      stop("Missing Ensemble_Rank column")
    }
    
    if (rule$threshold_type == "percent") {
      cutoff_rank <- ceiling(n * as.numeric(rule$threshold) / 100)
    } else if (rule$threshold_type == "rank") {
      cutoff_rank <- as.integer(rule$threshold)
    } else {
      stop("Ensemble rules support only 'percent' or 'rank'")
    }
    
    cutoff_rank <- max(1, min(cutoff_rank, n))
    
    pass_vec <- df$Ensemble_Rank <= cutoff_rank
    pass_vec[is.na(pass_vec)] <- FALSE
    return(pass_vec)
  }
  
  stop("Unknown rule type: ", rule$type)
}

#' Apply DiCE rules to the results table
#'
#' Helper function (not for users)
#'
#' @param dice_results_df A data frame containing DiCE results, including
#'   centrality columns, \code{Ensemble_Rank}, and \code{Phase}.
#' @param dice_rules A list of DiCE rules created using
#'   \code{dice_centrality_rule()} and/or \code{dice_ensemble_rule()}.
#' @param dice_logic Character string specifying how multiple rules are
#'   combined. Choose from \code{"AND"} or \code{"OR"}.
#'
#' @return The input data frame with added \code{Pass_Count},
#'   \code{Pass_Rules}, \code{DiCE_pass}, and updated \code{Phase} columns.
#' @noRd
apply_dice_rules <- function(dice_results_df, dice_rules, dice_logic = c("AND", "OR")) {
  
  dice_logic <- match.arg(dice_logic)
  
  if (length(dice_rules) == 0) {
    stop("dice_rules cannot be empty")
  }
  
  old_phase <- dice_results_df$Phase
  
  pass_mat <- lapply(dice_rules, function(rule) {
    evaluate_dice_rule(dice_results_df, rule)
  })
  
  pass_mat <- as.data.frame(pass_mat)
  
  rule_names <- vapply(dice_rules, function(rule) {
    if (rule$type == "centrality") {
      centrality_prefix(rule$metric)
    } else if (rule$type == "ensemble") {
      "Ensemble"
    } else {
      "Unknown"
    }
  }, character(1))
  
  colnames(pass_mat) <- rule_names
  
  pass_mat[] <- lapply(pass_mat, function(x) {
    x[is.na(x)] <- FALSE
    as.logical(x)
  })
  
  dice_results_df$Pass_Count <- rowSums(pass_mat, na.rm = TRUE)
  
  dice_results_df$Pass_Rules <- apply(pass_mat, 1, function(x) {
    passed <- colnames(pass_mat)[which(x)]
    if (length(passed) == 0) return(NA_character_)
    paste(passed, collapse = ",")
  })
  
  dice_results_df$DiCE_pass <- if (dice_logic == "AND") {
    apply(pass_mat, 1, all)
  } else {
    apply(pass_mat, 1, any)
  }
  
  dice_results_df$DiCE_pass[is.na(dice_results_df$DiCE_pass)] <- FALSE
  
  dice_results_df$Phase <- old_phase
  dice_results_df$Phase[dice_results_df$DiCE_pass] <- "DiCE"
  
  dice_results_df <- dice_results_df %>%
    dplyr::select(-DiCE_pass) %>% 
    relocate(Pass_Rules, Pass_Count, .after = EigenVector_rank)
  
  return(dice_results_df)
}

#' Format a DiCE rule as text
#'
#' Helper function (not for users)
#'
#' @param rule A DiCE rule created using \code{dice_centrality_rule()} or
#'   \code{dice_ensemble_rule()}.
#'
#' @return A character string describing the rule.
#' @noRd
format_dice_rule <- function(rule) {
  
  if (rule$type == "centrality") {
    metric_name <- centrality_prefix(rule$metric)
    
    if (rule$threshold_type == "mean") {
      return(paste0(metric_name, " >= mean in treatment OR control"))
    }
    
    if (rule$threshold_type == "percent") {
      return(paste0(metric_name, " top ", rule$threshold, "% in treatment OR control"))
    }
    
    if (rule$threshold_type == "rank") {
      return(paste0(metric_name, " top ", rule$threshold, " in treatment OR control"))
    }
  }
  
  if (rule$type == "ensemble") {
    if (rule$threshold_type == "percent") {
      return(paste0("Ensemble top ", rule$threshold, "%"))
    }
    
    if (rule$threshold_type == "rank") {
      return(paste0("Ensemble top ", rule$threshold))
    }
  }
  
  return("Unknown rule")
}

#' Format multiple DiCE rules as text
#'
#' Helper function (not for users)
#'
#' @param dice_rules A list of DiCE rules.
#' @param dice_logic Character string specifying how multiple rules are
#'   combined. Choose from \code{"AND"} or \code{"OR"}.
#'
#' @return A character string describing the combined DiCE rules.
#' @noRd
format_dice_rules_text <- function(dice_rules, dice_logic = c("AND", "OR")) {
  dice_logic <- match.arg(dice_logic)
  
  rule_txt <- vapply(dice_rules, format_dice_rule, character(1))
  
  paste(rule_txt, collapse = paste0(" ", dice_logic, " "))
}