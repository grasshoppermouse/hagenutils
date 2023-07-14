
#' Create a forest plot and table of coefficients and stats from one or more regression models.
#'
#' @param ... One or more model objects.
#' @param modelnames A character vector of names to label each plot and table column.
#' @param varnames A named character vector. Each element is a human friendly version of the variable name, and each name is the original variable name.
#' @param odds.ratio If TRUE, apply exponential transform to coefficients and omit intercept.
#' @param breaks An optional vector of breaks for the plot x-axis.
#' @param intercept If TRUE, include the intercept in the table and plot; otherwise omit.
#' @param facet If TRUE, display models in different facets. Otherwise, distinguish by color.
#' @param dodgewidth For position_dodge(). Default = 0.3.
#' @param sig The number of significant digits for the table.
#' @param draw If TRUE draw the plot; otherwise return the data. Default TRUE
#' @return A ggplot object.
#' @examples
#' data(mtcars)
#' m1 <- lm(mpg ~ hp, mtcars)
#' m2 <- lm(mpg ~ wt, mtcars)
#' mnms <- c('MPG vs HP', 'MPG vs Weight')
#' vnms <- c(mpg='Miles per gallon', hp='Horsepower', wt='Weight')
#' forestplot(m1, m2, modelnames=mnms, varnames=vnms)
#' @export 
forestplot <- function(
    ...,
    modelnames=NULL,
    varnames=NULL,
    stats=NULL,
    odds.ratio=F,
    breaks=NULL,
    intercept=T,
    facet=T,
    dodgewidth = 0.3,
    linewidth = 1,
    size = 1,
    draw = T,
    sig=3
    ){

    models <- list(...)

    if (is.null(modelnames)){
      modelnames = paste('Model', 1:length(models))
    }
    names(models) <- modelnames
    
    tidymodels <- purrr::list_rbind(purrr::map(models, function(x) broom::tidy(x, conf.int=T)), names_to = 'Model')

    if (odds.ratio){
        tidymodels$estimate <- exp(tidymodels$estimate)
        tidymodels$conf.low <- exp(tidymodels$conf.low)
        tidymodels$conf.high <- exp(tidymodels$conf.high)

        # ggplot params
        intr = 1
        xlabel = '\nOdds ratio (95% CI)'
        scle = ggplot2::scale_x_log10(breaks=breaks)
    }
    else {

        # ggplot params
        intr = 0
        xlabel = '\nEstimate (95% CI)'
        scle = ggplot2::scale_x_continuous()

    }

    if (odds.ratio | !intercept){
        # Filter out intercept terms
        tidymodels <- dplyr::filter(tidymodels, term != '(Intercept)')
    }

    if (is.null(varnames)) {
        varnames <- unique(tidymodels$term)
        names(varnames) <- varnames
    }

    tidymodels$term = factor(tidymodels$term, levels=rev(names(varnames)), labels=rev(varnames))
    
    if (!draw) return(tidymodels)
    
    if (facet){
      p <- 
        ggplot2::ggplot(tidymodels, ggplot2::aes(estimate, term, xmin=conf.low, xmax=conf.high)) +
        ggplot2::facet_wrap(~Model, ncol = length(models))
        
    } else
    {
      p <- 
        ggplot2::ggplot(tidymodels, ggplot2::aes(estimate, term, xmin=conf.low, xmax=conf.high, colour=Model, shape=Model))
    }
    
    p <- 
      p +
      ggplot2::geom_pointrange(size = size, linewidth = linewidth, position = position_dodge(width = dodgewidth)) +
      ggplot2::geom_vline(xintercept=intr, linetype='longdash') +
      scle +
      ggplot2::labs(y='', x=xlabel) +
      ggplot2::theme_bw()

  return(p)
}

#' logistic_forestplot
#'
#' @param ... One or more logistic regression model objects.
#' @param modelnames A character vector of names to label each plot and table column.
#' @param varnames A named character vector. Each element is a human friendly version of the variable name, and each name is the original variable name.
#' @param stats A named list of numeric vectors. The names are the names of the stats; the vectors are the stats to include in the table.
#' @param odds.ratio If TRUE, apply exponential transform to coefficients and omit intercept.
#' @param breaks An optional vector of breaks for the plot x-axis.
#' @param intercept If TRUE, include the intercept in the table and plot; otherwise omit.
#' @param sig The number of significant digits for the table.
#' @return A named list. $plot has the plot. $table has the data frame of coefficients and stats.
#' @examples
#' data(mtcars)
#' m1 <- glm(am ~ hp, family=binomial, mtcars)
#' m2 <- glm(mpg ~ wt, family=binomial, mtcars)
#' mnms <- c('AM vs HP', 'AM vs Weight')
#' vnms <- c(mpg='Miles per gallon', hp='Horsepower', wt='Weight')
#' out <- logistic_forestplot(m1, m2, modelsnames=mnms, varnames=vnms)
#' @export 
logistic_forestplot <- function(...,
                                modelnames = NULL,
                                varnames = NULL,
                                stats = NULL,
                                odds.ratio = F,
                                breaks = NULL,
                                intercept = T,
                                sig = 3) {
  # library(binomTools)
  models <- list(...)

  chisqr <- function(m) {
    chi.sqr <- m$null.deviance - m$deviance
    df <- m$df.null - m$df.residual
    p <- pchisq(chi.sqr, df, lower.tail = FALSE)
    stars = ''
    if (p < .05)
      stars = '*'
    if (p < .01)
      stars = '**'
    if (p < .001)
      stars = '***'
    return(paste0('Chisq (', signif(df, 3), ') = ', signif(chi.sqr, 3), stars))
  }

  # GOF stats
  chistats = character()
  tjurD = numeric()
  nullDev = character()
  residDev = character()

  for (m in models){
    chistats = c(chistats, chisqr(m))
    tjurD = c(tjurD, Rsq(m)$D)
    g <- broom::glance(m)
    nullDev = c(nullDev, paste0(signif(g$null.deviance, 3), ' (', g$df.null, ')'))
    residDev = c(residDev, paste0(signif(g$deviance, 3), ' (', g$df.residual, ')'))
  }

  tjurD = signif(tjurD, 2)

  stats <- list(
    "Tjur's D" = tjurD,
    "Null deviance (df)" = nullDev,
    "Residual deviance (df)" = residDev,
    "Chisqr" = chistats
  )

  return (
    forestplot(
      models,
      modelnames = modelnames,
      varnames = varnames,
      stats = stats,
      odds.ratio = odds.ratio,
      breaks = breaks,
      intercept = intercept,
      sig = sig
    )
  )

}

#' Output results of a t-test as a single string to be formatted inline by knitr
#'
#' @param formula Specify t-test using standard formula interface.
#' @param data A data frame containing the variables in formula.
#' @param effsize Include a user-supplied effect size (default = NULL)
#' @param sig The number of signficant digits for the output.
#' @return A character vector with 1 element.
#' @examples
#' data(mtcars)
#' m1 <- lm(mpg ~ hp, mtcars)
#' m2 <- lm(mpg ~ wt, mtcars)
#' mnms <- c('MPG vs HP', 'MPG vs Weight')
#' vnms <- c('Miles per gallon'=mpg, 'Horsepower'=hp, 'Weight'=wt)
#' out <- forestplot(m1, m2, modelsnames=mnms, varnames=vnms)
#' @export
inline_ttest <- function(ttest, effsize=NULL, sig=3){
  
  # ttest <- t.test(formula, data=data)
  tstat <- signif(ttest$statistic, sig)
  dfstat <- round(ttest$parameter, 1)
  p <- signif(ttest$p.value, sig)
  
  if (ttest$method == "Paired t-test"){
    m1 <- signif(ttest$estimate, sig)
    mean_txt <- glue::glue('Mean diff. = {m1}')
  } else {
    m1 <- signif(ttest$estimate[1], sig)
    m2 <- signif(ttest$estimate[2], sig)
    mean_txt <- glue::glue('M = {m2} vs. M = {m1}')
  }
  
  if (p < 2e-16){
    p <- glue::glue('$p<{fmt_sci(2.2e-16)}$')
  } else if (p < 0.001){
    p <- glue::glue('$p={fmt_sci(p)}$')
  } else {
    p <- glue::glue('$p={p}$')
  }
  
  if (is.null(effsize)){
    return(
      glue::glue('{mean_txt}, t({dfstat}) = {tstat}, {p}')
    )
  } else {
    return(
      glue::glue('{mean_txt}, t({dfstat}) = {tstat}, {p}, d = {effsize}')
    )
  }
}

#' Output summary table split by a categorical variable
#'
#' @param df data frame
#' @param vars Named vector. Names are variables to summarize, values are human readable versions.
#' @param facvar Factor variable defining exactly two groups.
#' @param statscol Include statistical test of location difference between two groups.
#' @param sig.digits The number of signficant digits for the output.
#' @return A data frame with the summary stats that can be printed with kable
#' @examples
#' data(mtcars)
#' vars = c(mpg = 'miles per gallon', hp = 'horse power')
#' custom.summarize(mtcars, vars, facvar='am')
#' @export 
custom.summarize = function(df, vars, facvar=NULL, statscol=T, test_type='wilcox', sig.digits=2){
  
  if(is.null(names(vars))) stop('vars variable must be a named variable')
  sum_stats <- function(df, var.names){
    
    d <- data.frame(
      N = sapply(var.names, FUN=function(x) sum(!is.na(df[[x]]))),
      Min = sapply(var.names, FUN=function(x) min(df[[x]], na.rm=T)),
      Max = sapply(var.names, FUN=function(x) max(df[[x]], na.rm=T)),
      Mean = sapply(var.names, FUN=function(x) mean(df[[x]], na.rm=T)),
      SD = sapply(var.names, FUN=function(x) sd(df[[x]], na.rm=T))
    )
    d[-1] <- sapply(d[-1], FUN=function(x) signif(x, sig.digits))
    # d[-1] <- map_df(d[-1], ~ signif(., sig.digits))
    
    # Now format -> char var
    d <- data.frame(
      N = d$N,
      Range = glue::glue_data(d, '{Min}-{Max}'),
      `Mean (SD)` = glue::glue_data(d, '{Mean} ({SD})'),
      stringsAsFactors = F, 
      check.names = F
    )
    row.names(d) <- NULL
    return(d)
  }

  location_test <- function(x, y, test = 't-test'){
    x <- as.numeric(x)
    y <- as.numeric(y)

    if (test == 't-test'){
      out <- t.test(x, y)
      return(list(statistic = out$statistic, p.value = out$p.value))
    }else{
      # wilcox test
      out <- wilcox.test(x, y)
      return(list(statistic = out$statistic, p.value = out$p.value))
    }

  }

  var.names = names(vars)
  d_vars = data.frame(Variable=as.character(vars), stringsAsFactors = F)

  if (!is.null(facvar)){
    facs = unique(df[[facvar]])
    if (length(facs)!=2) stop('facvar must have exactly two levels')
    df_list = list()
    for (i in 1:length(facs)){
      fac = facs[i]
      d_sub = df[df[[facvar]]==fac,]
      df_list[[i]] <- sum_stats(d_sub, var.names)
    }
    d <- do.call(cbind, df_list)
    if (statscol){
      d$d = NA
      # d$stat = NA
      d$p = NA
      for (i in 1:length(var.names)){
        df1 <- df[df[[facvar]]==facs[1],]
        df2 <- df[df[[facvar]]==facs[2],]
        testresults <- location_test(df1[[var.names[i]]], df2[[var.names[i]]], test = test_type)
        # d$stat[i] <- testresults$statistic
        d$p[i] <- signif(testresults$p.value, sig.digits)
        d$d[i] <- signif(cohen_d(df1[[var.names[i]]], df2[[var.names[i]]]), sig.digits) # d$t[i]*2/sqrt(ttest$parameter)
      }

    }
  }else{
    d <- sum_stats(df, var.names)
  }

  # d[-1] = sapply(d[-1], FUN=function(x) signif(x, sig.digits))

  d <- cbind(d_vars, d)
  row.names(d) <- NULL
  return(d)
}

#' @title summaryTable
#' @description Provides summary stats for a data frame, possibly by a two-level factor
#' @param d Data frame
#' @param vars Named vector. Names are variables to summarize, values are human readable versions., Default: NULL
#' @param fac Character value of grouping variable, Default: NULL
#' @param statscols Include Cohen's d and a Mann-Whitney test if a fac is provided, Default: T
#' @param sigfigs Number of significant figures to report, Default: 2
#' @return Either a single kable object, or a list of two kable objects
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  summaryTable(mtcars, fac="am", vars = c(mpg="Miles per gallon", hp="Horsepower"))
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}
#'  \code{\link[tibble]{tibble}}
#'  \code{\link[glue]{glue}}
#'  \code{\link[purrr]{map}}
#'  \code{\link[knitr]{kable}}
#'  \code{\link[kableExtra]{add_header_above}}
#' @rdname summaryTable
#' @export 
#' @importFrom dplyr select
#' @importFrom tibble tibble
#' @importFrom glue glue
#' @importFrom purrr map map_dfr
#' @importFrom knitr kable
#' @importFrom kableExtra add_header_above
summaryTable <- function(d, vars=NULL, fac=NULL, statscols=T, sigfigs=2, ...){
  require(dplyr)
  
  numeric_only <- function(d){
    d <- 
      d %>% 
      dplyr::select(where(~is.numeric(.x) | is.logical(.x)))
    if (!is.null(vars)){
      vdiff <- setdiff(names(vars), names(d))
      if (length(vdiff) > 0){
        vdiff <- paste(vdiff, collapse = ' ')
        stop(glue::glue('These vars not numeric or logical: {vdiff}'))
      }
      d <- d[names(vars)]
      names(d) <- vars
    }
    return(d)
  }
  
  statsummary <- function(v, sigfigs){
    v <- na.omit(v)
    tibble::tibble(
      N = length(v),
      Range = glue::glue('{signif(min(v), sigfigs)}-{signif(max(v), sigfigs)}'),
      `Mean (SD)` = glue::glue('{signif(mean(v), sigfigs)} ({signif(sd(v), sigfigs)})'),
    )
  }
  
  if (!is.null(fac)){
    facs <- sort(as.character(unique(d[[fac]])))
    if (length(facs) != 2) stop('fac must have exactly two unique values')
    dfs <- split(d, d[fac])[facs] # Make sure order in dfs is same as facs
  } else {
    dfs <- list(d)
  }
  
  dfs <- map(dfs, numeric_only)
  
  thesummary <- purrr::map(
    dfs,
    ~purrr::map_dfr(.x, ~statsummary(.x, sigfigs), .id = 'Variable')
  )
  
  if (!is.null(fac) & statscols){
    thestats <- map2_dfr(dfs[[1]], dfs[[2]], ~cohen_d.default(.x, .y, test=T, sig=sigfigs), .id = 'Variable')
    thesummary[[2]] <- left_join(thesummary[[2]], thestats)[-1]
  }
  
  if (is.null(fac)){
    return(knitr::kable(thesummary[[1]], ...))
  } else {
    header1 <- setNames(c(1, 3), c(" ", facs[1]))
    header2 <- setNames(c(3, 2), c(facs[2], " "))
    tables <- list(
      knitr::kable(thesummary[[1]], ...) %>% kableExtra::add_header_above(header1),
      knitr::kable(thesummary[[2]], ...) %>% kableExtra::add_header_above(header2)
    )
    return(tables)
  }
}

#' @title fmt_sci
#' @description Generate latex for scientific notation
#' @param x Numeric vector to be formatted
#' @return latex character vector with x in scientific notation
#' @details For use in inline markdown.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  fmt_sci(2.34e-9)
#'  }
#' }
#' @seealso 
#'  \code{\link[glue]{glue}}
#' @rdname fmt_sci
#' @export 
#' @importFrom glue glue
fmt_sci <- function(x){
  a <- floor(log10(abs(x)))
  b <- x/(10^a)
  return(glue::glue('{b} \\times 10^{{{a}}}'))
}

#' @title fmt_terms
#' @description Format estimate and 95\% CIs for inline markdown
#' @param models A named list of models
#' @param str A glue format string, Default: '$\\\\beta={signif(estimate, 2)}$ ({signif(conf.low, 2)}, {signif(conf.high, 2)})'
#' @return A name list of model statistics, including str, a string version of the term estimate.
#' @details Provide a named list of regression models, and the resulting named list will have a formatted string suitable for inline rmarkdown
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  stats <- fmt_terms(list(m1 = lm(mpg ~ hp, mtcars), m2 = lm(mpg ~ wt, mtcars)))
#'  stats$m1$hp$str
#'  }
#' }
#' @rdname fmt_terms
#' @export
fmt_terms <- function(models, str = "$\\beta={signif(estimate, 2)}$ ({signif(conf.low, 2)}, {signif(conf.high, 2)})"){
  tdy <- function(m){
    tidy(m, conf.int = T) %>% 
      mutate(
        str = glue(str)
      ) %>% 
      split(.$term)
  }
  map(models, tdy)
}

#' @title renumber_tweets
#' @description Replace #. with sequential numbering (e.g., 1., 2.,)
#' @param fn file name of markdown file
#' @return Writes to the input file, fn!
#' @details Use when composing twitter thread in markdown. Number tweets with #. placeholder.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  renumber_tweets("thread.md")
#'  }
#' }
#' @rdname renumber_tweets
#' @export 

renumber_tweets <- function(fn){
  f <- readr::read_lines(fn)
  n = 1
  for (i in 1:length(f)){
    if (stringr::str_starts(f[i], "(\\d+\\. )|(#\\. )")){
      f[i] <- stringr::str_replace(f[i], "(\\d+\\. )|(#\\. )", paste0(n, ". "))
      n <- n+1
    }
  }
  readr::write_lines(f, fn)
}
