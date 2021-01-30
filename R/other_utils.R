
library(tidyverse)

#' Predicted probability from svyglm model and new data, along with 95\% CI.
#'
#' @param mod A svyglm model object.
#' @param df A data frame of new data created with, e.g., expand.grid.
#' @return df with three additional columns:
#' \describe{
#' \item{Probability}{The predicted probability at each value in df}
#' \item{LL}{The lower 95\% confidence limit}
#' \item{UL}{The upper 95\% confidence limit}
#' }
#' @examples
#' m = glm(am ~ mpg + hp + cyl, family=binomial, data=mtcars)
#' newdata = expand.grid(mpg=10:35, hp=50:335, cyl=c(4,6,8))
#' df = predicted_prob(m, newdata)
#' @export 
predicted_prob = function(mod, df){

  df2 = cbind(df, predict(mod, newdata = df, type = "link", se=T))

  return(
    within(df2, {
      Probability = plogis(link)
      LL = plogis(link - (1.96 * SE))
      UL <- plogis(link + (1.96 * SE))
    })
  )
}

#' Cohen's d
#'
#' @param f A two-sided formula: y ~ x.
#' @param d A data frame.
#' @return Cohen's d
#' @examples
#' cohen_d(mpg ~ am, mtcars)
#' @export 
cohen_d = function(f, data, sig = 3){
  mf <- model.frame(f, data = data)
  if (ncol(mf) != 2) stop('Formula must contain only 2 variables: y ~ x')
  mf[[2]] <- c(mf[[2]]) # convert factors to regular vectors
  if (length(table(mf[[2]])) != 2) stop('factor must have exactly two values')
  if (!mode(mf[[1]]) %in% c('integer', 'numeric')) stop('Left hand side must be a numeric variable')
  names(mf) <- c('y', 'x')
  x <-
    mf %>% 
    group_by(x) %>% 
    summarise(n = n(), mean = mean(y, na.rm = T), sd2 = sd(y, na.rm = T)^2)
  d <- (x$mean[1] - x$mean[2])/sqrt(((x$n[1]-1)*x$sd2[1] + (x$n[2]-1)*x$sd2[2])/(sum(x$n)-2))
  return(signif(d, sig))
}

#' pca_loadings_plot
#'
#' A lollipop plot of PCA variable loadings,
#' as provided by prcomp, e.g., PC1, PC2, ...
#'
#' @param obj Output of prcomp
#' @param components A numeric vector of components to plot
#' @param sortby Sort variables by loadings on this component (index to component vector)
#' @return ggplot object
#' @export
#' @examples
#' m <- prcomp(mtcars, scale = T)
#' pca_loadings_plot(m, components = c(1,2))
#' @export 
pca_loadings_plot <- function(obj, components = 1:3, sortby = 1){
    require(dplyr)
    varprop <- summary(obj)$importance[2,components] # Get % variance
    nms <- names(varprop)
    varprop <- paste0(nms, ' (', round(varprop*100, 1), '%)')
    names(varprop) <- nms
    loadings <- obj$rotation[,components, drop = F]
    loadings %>%
    data.frame %>%
    dplyr::mutate(Variable = forcats::fct_reorder(rownames(.), loadings[, sortby])) %>%
    tidyr::gather(key = PC, value = Loading, -Variable) %>%
    dplyr::mutate(PC = varprop[PC]) %>%
    ggplot2::ggplot() +
    ggalt::geom_lollipop(ggplot2::aes(Loading, Variable, colour = Loading), size = 1, horizontal = T) +
    ggplot2::scale_color_gradient2(low = 'red', mid = 'white', 'high' = 'blue') +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::facet_wrap(~PC) +
    ggplot2::theme_bw(15)
}

#' @title tjurD
#' @description Computes Tjur's D
#' @param m A glm object with family = binomial
#' @return A list with two elements, rsq = Tjur's D, and plot = ggplot object.
#' @details Computes Tjur's D and produces a probability plot
#' @examples
#' \dontrun{
#' if(interactive()){
#'  m <- glm(am ~ mpg, family = binomial, mtcars)
#'  out <- tjurD(m)
#'  out$rsq
#'  out$plot
#'  }
#' }
#' @rdname tjurD
#' @export
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap
tjurD <- function(m){
  y <- m$model[[1]] # Outcome var
  probs <- m$fitted.values
  rsq <- mean(probs[y == 1]) - mean(probs[y == 0])
  df <- data.frame(
    Probability = probs,
    y = y
  )
  list(
    rsq = rsq,
    plot = ggplot2::ggplot(df, ggplot2::aes(Probability)) +
      ggplot2::geom_histogram() +
      ggplot2::facet_wrap(~y, ncol = 1)
  )
}


#' svysmooth2df
#'
#' Extracts the x and y values from one or more svysmooth objects,
#' from the survey package, and returns them in a data frame.
#'
#' @param ... One or more svysmooth objects with names.
#'
#' @return A data frame (tibble) with 3 columns: x, y, and a
#' character vector indicating which values belong to which smooth
#' @details If multiple smooths are supplied, the x and y variables
#' must be the same, but computed over different design objects (e.g.,
#' different subsets of a survey).
#' @export
#' @examples
#' library(survey)
#' data(api)
#' dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
#' dclus2<-subset(dclus1, sch.wide=="Yes" & comp.imp=="Yes")
#' m1 <- svysmooth(api00~ell,dclus1)
#' m2 <- svysmooth(api00~ell,dclus2)
#' df <- svysmooth2df(full = m1, subset = m2)
svysmooth2df <-
  function(...){

    smooths <- list(...)
    if(is.null(names(smooths))) stop("Each smooth must have a name")

    x <- purrr::map(smooths, list(1, 'x'))
    y <- purrr::map(smooths, list(1, 'y'))
    Smooth <- rep(names(smooths), times = purrr::map_int(x, length))

    df <-
      tibble::data_frame(
        flatten_dbl(x),
        flatten_dbl(y),
        Smooth
      )

    names(df)[1] <- names(smooths[[1]])
    names(df)[2] <- attr(smooths[[1]], 'ylab')

    return(df)
  }


#' @title ggdotchart
#' @description ggplot version of base::dotchart
#' @param v A named vector, e.g., a table
#' @param threshold A numeric value. Values in v that are greater than this will have different colors and shapes
#' @return a ggplot object
#' @details Uses theme_minimal, and no axis labels.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  ggdotchart(table(mtcars$cyl))
#'  }
#' }
#' @rdname ggdotchart
#' @export
ggdotchart <- function(v, threshold = NULL){
  nms <- names(v)
  d <- tibble(
    x = as.numeric(v),
    y = factor(nms, levels = nms[order(v)])
  ) 
  
  if (!is.null(threshold) & is.numeric(threshold)){
    d$threshold <- d$x > threshold
    p <- ggplot(d, aes(x, y, colour = threshold, shape = threshold))
  } else {
    p <- ggplot(d, aes(x, y))
  }
  
  p +
    geom_point(size = 3) +
    labs(x = "", y = "") +
    theme_minimal(15)
}

#' @title scale_colour_binary
#' @description Provides exactly two colors from the viridis magma palette (red and purple)
#' @param direction Order of colors, 1: purple-red. -1: red-purple, Default: 1
#' @param ... Additional arguments to be passed to ggplot2::discrete_scale
#' @return a colour scale to be used with ggplot2
#' @details The default top level from the discrete viridis scale is bright yellow, which is too bright. This scale uses red for top level.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  ggplot(mtcars, aes(hp, mpg, colour=factor(am))) + geom_point() + scale_colour_binary()
#'  }
#' }
#' @seealso 
#'  \code{\link[ggplot2]{discrete_scale}}
#'  \code{\link[viridisLite]{magma}}
#' @rdname scale_colour_binary
#' @export 
#' @importFrom ggplot2 discrete_scale
#' @importFrom viridisLite magma
scale_colour_binary<- function(direction=1, ...){
  ggplot2::discrete_scale(
    "colour", "binary", 
    function(n){
      if (n != 2) stop('custom binary palette provides exactly 2 colors')
      cols <- viridisLite::magma(11)[c(4,8)]
      if (direction >= 0) cols else rev(cols)
    }, 
    ...
  )
}

#' @rdname scale_colour_binary
#' @export
scale_color_binary <- scale_colour_binary

#' @title scale_fill_binary
#' @description Provides exactly two colors from the viridis magma palette (red and purple)
#' @param direction Order of colors, 1: purple-red. -1: red-purple, Default: 1
#' @param ... Additional arguments to be passed to ggplot2::discrete_scale
#' @return a fill scale to be used with ggplot2
#' @details The default top level from the discrete viridis scale is bright yellow, which is too bright. This scale uses red for top level.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  ggplot(mtcars, aes(mpg, fill=factor(am))) + geom_histogram() + scale_fill_binary()
#'  }
#' }
#' @seealso 
#'  \code{\link[ggplot2]{discrete_scale}}
#'  \code{\link[viridisLite]{magma}}
#' @rdname scale_fill_binary
#' @export 
#' @importFrom ggplot2 discrete_scale
#' @importFrom viridisLite magma
scale_fill_binary<- function(direction=1, ...){
  ggplot2::discrete_scale(
    "fill", "binary", 
    function(n){
      if (n > 2) stop('custom binary palette provides only 2 colors')
      cols <- viridisLite::magma(11)[c(4,8)]
      if (direction >= 0) cols else rev(cols)
    }, 
    ...
  )
}

#' @title ggemmeans
#' @description ggplot of estimated marginal means object, with categorical spec variable sorted by estimate
#' @param em An emmeans object
#' @param reorder If TRUE (default), reorder levels of a categorical predictor by values of the estimate
#' @return a ggplot object
#' @details This function is similar to the emmeans plot function, except that it sorts the levels of a categorical variable by the estimated mean
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  m <- lm(height ~ eye_color, starwars)
#'  em <- emmeans(m, spec = 'eye_color')
#'  ggemmeans(em)
#'  }
#' }
#' @rdname ggemmeans
#' @export 
ggemmeans <- function(em, reorder = T){
  emsum <- summary(em)
  estName <- attr(emsum, 'estName')
  clNames <- attr(emsum, 'clNames')
  var <- attr(emsum, 'pri.vars')[1]
  if (reorder & is.factor(emsum[[var]])){
    emsum[var] <- forcats::fct_reorder(emsum[[var]], emsum[[estName]])
  }
  ggplot(emsum, aes_string(estName, var, xmin = clNames[1], xmax = clNames[2])) +
    geom_errorbarh(height = 0, lwd = 2.5, alpha = 0.3) +
    geom_point() +
    theme_minimal(15)
}

#' @title hagenheat
#' @description Basic heatmap using ggplot, dist, and seriation. No dendrograms are plotted.
#' @param d A data frame, matrix, or dist object. Column 1 of a data frame can be row labels. Remaining columns must be numeric.
#' @param method "seriate" or "hclust". Default: "seriate"
#' @param seriation_method If not specified, 'Spectral' will be used for dist objects, and 'PCA' for matrices and data frames.
#' @param independent Whether to seriate a matrix using matrix methods (F) or separate distance matrices for rows and columns (T). Default: F
#' @param hc_method Agglomeration method from hclust, Default: 'ward.D'
#' @param dist Distance method from dist, Default: 'euclidean'
#' @param scale. Whether to scale rows ('row'), columns ('column'), or neither ('none'), Default: 'none'
#' @param viridis_option One of the viridis color options, 'A', 'B', 'C', 'D', 'E', Default: 'D'
#' @param wrap_labels Logical. Wrap the x-axis labels. Default: T
#' @param rotate_labels Logical. Rotate the x-axis labels by 90 degrees. Default: F
#' @param ann_col A data frame with two variables. The first column must be the names of the columns of 'd'. The second column contains the annotation values.
#' @return A ggplot object
#' @details Produces a very simple ggplot heatmap using viridis colors. 
#' Rows and columns are ordered using the \code{seriation} package.
#' For data frames, first column must be row labels. Remaining columns must be numeric.
#' Scaling, if performed, occurs after seriation, and is therefore only for display purposes.
#' Seriation methods for objects when \code{independent=T}: \code{list_seriation_methods('dist')}.
#' Seriation methods for matrices and data frames when \code{independent=F}: \code{list_seriation_methods('matrix')}.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  hagenheat(mtcars, scale. = 'col')
#'  }
#' }
#' @seealso 
#'  \code{\link[seriation]{seriate}}
#'  \code{\link[seriation]{seriation_methods}}
#'  \code{\link[ggnewscale]{new_scale_fill}}
#' @rdname hagenheat
#' @export 
#' @importFrom seriation seriate get_order
#' @importFrom dplyr bind_cols
#' @importFrom tidyr gather
#' @importFrom scales label_wrap
#' @importFrom ggnewscale new_scale_fill
hagenheat <- function(
  d, 
  method='seriate', 
  seriation_method=NULL, 
  independent = F,
  hc_method = 'ward.D', 
  dist_method = 'euclidean', 
  scale. = "none", 
  viridis_option = 'D',
  wrap_labels = T,
  rotate_labels = F,
  ann_col = NULL
  ){
  
  if (inherits(d, 'dist')){
    
    rwnms <- labels(d)
    if (is.null(seriation_method)) seriation_method <- 'Spectral'
    o <- seriation::seriate(d, method = seriation_method)
    row_order <- seriation::get_order(o)
    col_order <- row_order
    d <- as.matrix(d)
    if (is.null(rwnms)) {
      warning('Lack of labels on dist object might make plot difficult to interpret')
      rwnms <- as.character(1:nrow(d))
    }
    
  } else {
    
    # Convert data frame to matrix
    rwnms <- character(0)
    if (inherits(d, 'data.frame')){
      if (is.character(d[[1]]) | is.factor(d[[1]])){
        rwnms <- d[[1]]
        d <- d[-1]
      }
      d <- as.matrix(d)
    }
    
    # If d is logical, converts to numeric
    if (mode(d) == 'logical') d <- d*1
    if (mode(d) != 'numeric') stop('d must be convertible to a numeric matrix')
    
    if (length(di <- dim(d)) != 2) stop("'d' must have 2 dimensions")
    
    nr <- di[1L]
    nc <- di[2L]
    
    if (nr <= 1 || nc <= 1) stop("'d' must have at least 2 rows and 2 columns")
    
    # Get rownames if haven't gotten them already
    if (length(rwnms) == 0){
      if(length(rownames(d))>1){
        rwnms <- rownames(d)
      } else {
        rwnms <- as.character(1:nrow(d))
      }
    }
    
    # Order the rows and columns
    if (method == 'seriate'){
      
      if (is.null(seriation_method)) seriation_method <- 'PCA'
      
      if (independent){
        
        o <- seriation::seriate(dist(d, method = dist_method), method=seriation_method)
        row_order <- seriation::get_order(o)
        o <- seriation::seriate(dist(t(d), method = dist_method), method=seriation_method)
        col_order <- seriation::get_order(o)
        
      } else {
        
        o <- seriation::seriate(d, method=seriation_method)
        row_order <- seriation::get_order(o, dim = 1)
        col_order <- seriation::get_order(o, dim = 2)
        
      }
      
    } else if (method == 'hclust'){
      
      hclustrows <- hclust(dist(d, method = dist_method), method = hc_method)
      row_order <- hclustrows$order
      hclustcols <- hclust(dist(t(d), method = dist_method), method = hc_method)
      col_order <- hclustcols$order
      
    } else {
      
      stop("method must be 'seriate' or 'hclust'")
      
    }
  }
  
  if (scale. == 'row'){
    d <- as_tibble(t(scale(t(d))))
  } else if (scale. == 'column'){
    d <- as_tibble(scale(d)) 
  } else {
    d <- as_tibble(d)
  }
  
  rwnms <- factor(rwnms, levels = rwnms[row_order])
  d <- dplyr::bind_cols(rowname=rwnms, d)
  
  p <-
    d %>%
    tidyr::gather(key = key, value = value, -1) %>% 
    mutate(
      key = factor(key, levels = colnames(d[-1])[col_order]),
    ) %>%
    ggplot(aes_string('key', 'rowname', fill = 'value')) +
    geom_tile() +
    scale_fill_viridis_c(option = viridis_option) +
    labs(x = "", y = "") +
    theme_minimal()
  
  if (wrap_labels) p <- p + scale_x_discrete(labels = scales::label_wrap(10))
  if (rotate_labels) p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  if (!is.null(ann_col)){
    
    nms <- names(ann_col)
    ann_col$value = 1
    h <- nrow(d)/20
    
    p <- p + 
      ggnewscale::new_scale_fill() +
      geom_tile(data = ann_col, aes_string(x = nms[1], y = nrow(d) + h, fill = nms[2]), height = h, width = 1) +
      coord_cartesian(clip = 'off') +
      theme(
        plot.margin = unit(c(2,1,1,1), "lines")
      )
  }
  return(p)
}
