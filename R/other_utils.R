
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

cohen_d <- function(x, ...) UseMethod("cohen_d")

cohen_d.default <- function(x, y, test=F, sig = 3, ...){
  if (!mode(x) %in% c('integer', 'numeric')) stop('x must be a numeric variable')
  if (!mode(y) %in% c('integer', 'numeric')) stop('y must be a numeric variable')
  
  x <- na.omit(x)
  nx <- length(x)
  if (nx < 2) stop("x must have 2 or more non-missing values")
  
  y <- na.omit(y)
  ny <- length(y)
  if (ny < 2) stop("y must have 2 or more non-missing values")
  
  d <- signif((mean(x) - mean(y))/sqrt(((nx-1)*sd(x)^2 + (ny-1)*sd(y)^2)/(sum(nx, ny)-2)), sig)
  if (test){
    p.value <- suppressWarnings(wilcox.test(x, y)$p.value)
    p.value <- signif(p.value, sig)
    return(list(d=d, 'p-value'=p.value))
  }
  return(d)
}

#' Cohen's d
#'
#' @param f A two-sided formula: y ~ x.
#' @param d A data frame.
#' @param test Whether to include a Wilcoxon (Mann-Whitney) test of location difference (default: F)
#' @param sig The number of significant digits to report (default: 3)
#' @return Cohen's d
#' @examples
#' cohen_d(mpg ~ am, mtcars)
#' @export 
cohen_d.formula = function(f, data, ...){
  
  mf <- model.frame(f, data = data)
  if (ncol(mf) != 2) stop('Formula must contain only 2 variables: y ~ x')
  facs <- unique(mf[[2]])
  if (length(facs) != 2) stop('factor must have exactly two values')
  if (!mode(mf[[1]]) %in% c('integer', 'numeric')) stop('Left hand side must be a numeric variable')
  v1 <- mf[[1]][mf[[2]]==facs[1]]
  v2 <- mf[[1]][mf[[2]]==facs[2]]
  do.call("cohen_d", list(v1, v2, ...))
}

#' pca_loadings_plot
#'
#' A lollipop plot of PCA variable loadings,
#' as provided by prcomp, e.g., PC1, PC2, ...
#'
#' @param obj Output of prcomp
#' @param components A numeric vector of components to plot
#' @param sortby Sort variables by loadings on this component (index to component vector)
#' @param threshold Remove variables with sum(abs(loadings)) below threshold (default=0)
#' @param reverse Vector of components that will be multiplied by -1
#' @return ggplot object
#' @examples
#' m <- prcomp(mtcars, scale = T)
#' pca_loadings_plot(m, components = c(1,2))
#' @export 
pca_loadings_plot <- function(obj, components = 1:3, sortby = 1, threshold = 0, reverse=NULL){
  require(dplyr)
  if (!all(reverse %in% components)) stop("reverse vector must be a subset of components vector")
  pca_summ <- summary(obj)
  comp_nms <- colnames(pca_summ$importance)
  varprop <- pca_summ$importance[2,components] # Get % variance
  nms <- comp_nms[components]
  varprop <- paste0(nms, ' (', round(varprop*100, 1), '%)')
  names(varprop) <- nms
  
  pcs <- colnames(obj$rotation)[components]
  loadings <- 
    obj$rotation[,components, drop = F] %*% diag(obj$sdev[components]) %>%
    data.frame %>% 
    setNames(pcs)
  
  if(!is.null(reverse)) loadings[reverse] <- -loadings[reverse]
  
  loadings <- loadings[rowSums(abs(loadings)>threshold)>0,] 
  
  loadings %>%
    dplyr::mutate(Variable = forcats::fct_reorder(rownames(.), loadings[, sortby])) %>%
    tidyr::gather(key = PC, value = Loading, -Variable) %>%
    dplyr::mutate(PC = varprop[PC]) %>%
    ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x=0, xend=Loading, y=Variable, yend=Variable, colour = Loading), size=2) +
    ggplot2::geom_point(ggplot2::aes(x=Loading, y=Variable, colour = Loading), size=3) +
    ggplot2::scale_color_gradient2(low = viridisLite::magma(11)[8], mid = 'white', 'high' = viridisLite::magma(11)[4]) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::facet_wrap(~PC) +
    ggplot2::labs(x = '\nStandardized loadings', y = '') +
    ggplot2::theme_bw(15)
}

#' @title pca_biplot
#' @description Create a ggplot biplot of a prcomp obj
#' @param obj prcomp obj
#' @param components The PCs to plot, specified as a numeric vector of length 2, Default: c(1, 2)
#' @param data Optional data frame with the same number of rows as that used to compute pca, Default: NULL
#' @param threshold Only display variables with loadings exceeding this value, Default: 0
#' @param reverse Vector of components that will be multiplied by -1
#' @param label_size text size of the eigenvector labels, Default: 5
#' @param geom_point whether to add geom_point to plot, Default: T
#' @return ggplot obj
#' @details Creates a ggplot obj from a prcomp obj. Use threshold to omit variables with small loadings. Add a data frame to create a more complex ggplot. Possibly set geom_point=F to add a geom_point with custom aesthetics.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  m <- prcomp(mtcars, scale. = T)
#'  pca_biplot(m)
#'  pca_biplot(m, data = mtcars, geom_point=F) + geom_point(aes(size = mpg))
#'  }
#' }
#' @rdname pca_biplot
#' @export 
pca_biplot <- function(obj, components = c(1,2), data = NULL, threshold = 0, reverse=NULL, label_size = 5, geom_point=T) {
  
  if (class(obj) != 'prcomp') stop('obj class must be prcomp')
  if (length(components) != 2 | mode(components) != 'numeric') stop('components is not a numeric vector of length 2')
  if (max(components > ncol(obj$x)) | min(components) < 1) stop(paste('components must be between 1 and ', num_pc))
  if (!all(reverse %in% components)) stop("reverse vector must be a subset of components vector")
  
  pcs <- colnames(obj$x)
  if (!is.null(reverse)){
    obj$x[,reverse] <- -obj$x[,reverse]
    obj$rotation[,reverse] <- -obj$rotation[,reverse]
  }
  
  pcX <- pcs[components[1]]
  pcY <- pcs[components[2]]
  
  if (!is.null(data)){
    if (!'data.frame' %in% class(data)) stop('data must be a data frame')
    if (nrow(data) != nrow(obj$x)) stop('the number of rows of data, ', nrow(data), ', does not equal the number of rows of x, ', nrow(obj$x))
    if (pcX %in% names(data) | pcY %in% names(data)) stop(paste('data cannot contain variables', pcX, 'or', pcY))
    data[pcX] <- obj$x[,pcX]
    data[pcY] <- obj$x[,pcY]
  } else {
    data = data.frame(obj$x[,components])
  }
  
  pct_var <- paste0('(', round(100 * summary(obj)$importance[2,components], 1), '%)')
  pct_var <- paste(colnames(obj$x)[components], pct_var)
  
  plot <- 
    ggplot(data, aes_string(x=pcX, y=pcY)) + 
    geom_hline(aes(yintercept = 0), size=.2) + 
    geom_vline(aes(xintercept = 0), size=.2) +
    labs(x = pct_var[1], y = pct_var[2])
  
  if (geom_point){
    plot <- plot + geom_point()
  }
  
  datapc <- data.frame(varnames=rownames(obj$rotation), obj$rotation[,components])
  datapc <- datapc[rowSums(abs(datapc[-1])>threshold)>0,] 
  mult <- min(
    (max(data[,pcY]) - min(data[,pcY])/(max(datapc[,pcY])-min(datapc[,pcY]))),
    (max(data[,pcX]) - min(data[,pcX])/(max(datapc[,pcX])-min(datapc[,pcX])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(pcX)),
                      v2 = .7 * mult * (get(pcY))
  )
  plot + 
    coord_equal() + 
    geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = label_size, vjust=1, color="red") + 
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red") +
    theme_minimal()
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
  v <- sort(unique(y)) # binary levels of y
  probs <- m$fitted.values
  rsq <- mean(probs[y == v[2]]) - mean(probs[y == v[1]])
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
        purrr::flatten_dbl(x),
        purrr::flatten_dbl(y),
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
  d <- tibble::tibble(
    x = c(v), # remove table class
    y = forcats::fct_reorder(names(v), v)
  ) 
  
  if (!is.null(threshold) & is.numeric(threshold)){
    d$threshold <- d$x > threshold
    p <- ggplot2::ggplot(d, ggplot2::aes(x, y, colour = threshold, shape = threshold))
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(x, y))
  }
  
  if (min(v) >= 0) p <- p + ggplot2::xlim(0, NA)
  
  p +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme_minimal(15)
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
#' @param by String indicating "by" term (conditioning variable)
#' @param reorder If TRUE (default), reorder levels of a categorical predictor by values of the estimate
#' @return a ggplot object
#' @details This function is similar to the emmeans plot function, except that it sorts the levels of a categorical variable by the estimated mean
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  data(starwars, package='dplyr')
#'  m <- lm(height ~ eye_color, starwars)
#'  em <- emmeans(m, specs = 'eye_color')
#'  ggemmeans(em)
#'  }
#' }
#' @rdname ggemmeans
#' @export 
ggemmeans <- function(em, by = NULL, reorder = T){
  emsum <- summary(em)
  estName <- attr(emsum, 'estName')
  clNames <- attr(emsum, 'clNames')
  var <- attr(emsum, 'pri.vars')[1]
  if (reorder & is.factor(emsum[[var]])){
    emsum[var] <- forcats::fct_reorder(emsum[[var]], emsum[[estName]])
  }
  
  if (is.null(by)){
    p <- ggplot(emsum, aes_string(estName, var, xmin = clNames[1], xmax = clNames[2]))
  } else {
    p <- ggplot(emsum, aes_string(estName, var, xmin = clNames[1], xmax = clNames[2], colour = by))
  }
  
  p +
    geom_pointrange(lwd = 2.5, alpha = 0.7, position = position_dodge(width = 0.25), fatten = 1) +
    # geom_point(position = position_dodge(width = 0.25)) +
    guides(colour = guide_legend(override.aes = list(size = 1))) +
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
#' @param rev_col Reverse column order. Default: F
#' @param rev_row Reverse row order. Default: F
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
#'  hagenheat(mtcars, scale. = 'column')
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
  ann_col = NULL,
  rev_col = F,
  rev_row = F
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
  
  if (rev_row) row_order <- rev(row_order)
  if (rev_col) col_order <- rev(col_order)
  
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

#' @title regressiontable
#' @description Create gt table of regression coefficients
#' @param models Named list of lm and glm models.
#' @param caption Optional caption for table, Default: NULL
#' @param sigfig The number of significant digits to display, Default: 3
#' @return A gt table
#' @details Currently only supports lm and glm models
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  models <- list(
#'  "model 1" = lm(mpg ~ hp, mtcars),
#'  "model 2" = glm(am ~ hp, family=binomial, mtcars)
#'  )
#'  regressiontable(models)
#'  }
#' }
#' @seealso 
#'  \code{\link[broom]{tidy}}
#'  \code{\link[gt]{gt}}
#' @rdname regressiontable
#' @export 
#' @importFrom purrr map_df map
#' @importFrom broom tidy glance
#' @importFrom stringr str_glue_data
#' @importFrom gt gt cols_label fmt_number fmt_scientific tab_footnote tab_style cells_row_groups cell_text cells_body
regressiontable <- function(models, caption = NULL, sigfig = 3){
  
  for (m in models){
    if (!class(m)[1] %in% c('lm', 'glm', 'lmerModLmerTest')) stop('only lm, glm, and lmerModLmerTest models are supported')
  }
  
  tidy2 <- function(m){
    x <- broom::tidy(m, conf.int = T)
    if (class(m)[1] == 'lmerModLmerTest'){
      x <- x[x$effect == 'fixed',]
      x$effect <- NULL
      x$group <- NULL
      x$df <- NULL
    }
    return(x)
  }
  
  glance2 <- function(m){
    x <- broom::glance(m)
    if (class(m)[1] == 'lmerModLmerTest'){
      x$nobs <- nobs(m)
    } 
    return(x)
  }
  
  model_stats <- purrr::map_df(models, ~tidy2(.), .id = 'Model')
  names(model_stats) <- c('Model', 'Variable', 'Estimate', 'Std.Err', 'Statistic', 'P-value', 'Lower 95% CI', 'Upper 95% CI')
  
  glue_dict <- c(
    'lm' = 'N={nobs}; Rsq={r.squared}; Adj.Rsq={adj.r.squared}; F({df},{df.residual})={statistic}; p={p.value}',
    'glm' = 'N={nobs}; Null deviance={null.deviance} on {df.null} df; Residual deviance={deviance} on {df.residual} df',
    'lmerModLmerTest' = 'N={nobs}; Sigma={sigma}; df.residual={df.residual}'
  )
  
  model_summaries <-
    purrr::map(models, ~stringr::str_glue_data(signif(glance2(.), sigfig), glue_dict[class(.)[1]]))
  names(model_summaries) <- names(models)

  thetable <- model_stats %>%
    gt::gt(groupname_col = 'Model', caption = caption) %>%
    gt::cols_label(Variable = '') %>% 
    gt::fmt(columns = c(3:5, 7, 8), fns = function(x) signif(x, sigfig)) %>% 
    # gt::fmt_number(columns = c(4:5, 7:8), n_sigfig = sigfig) %>%
    gt::fmt_scientific(columns = 6, decimals = 1) %>%
    gt::tab_style(
      style = gt::cell_text(indent = gt::px(40)),
      locations = gt::cells_body(columns = 'Variable')
    )
  
  for(i in 1:length(model_summaries)){
    thetable <- gt::tab_footnote(thetable, model_summaries[[i]], location=gt::cells_row_groups(groups = names(model_summaries[i])))
  }
  return(thetable)
}

#' @title ggmediation
#' @description Plot ACME and ADE from objects produced by mediate function in mediation package
#' @param obj An object produced by the mediate function from the mediation package.
#' @return A ggplot
#' @details Produces a ggplot version of the mediation plot from the mediation package.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  library(mediation)
#'  m_med <- lm(hp ~ wt, mtcars)
#'  m_out <- lm(mpg ~ hp + wt, mtcars)
#'  m <- mediate(m_med, m_out, treat='wt', mediator='hp')
#'  ggmediation(m)
#'  }
#' }
#' @rdname ggmediation
#' @export 
ggmediation <- function(obj){
  d <-
    tibble(
      ACME.point = obj$d.avg,
      ACME.low = obj$d.avg.ci[1],
      ACME.high = obj$d.avg.ci[2],
      ADE.point = obj$z.avg,
      ADE.low = obj$z.avg.ci[1],
      ADE.high = obj$z.avg.ci[2],
      Total.point = obj$tau.coef,
      Total.low = obj$tau.ci[1],
      Total.high = obj$tau.ci[2]
    ) %>%
    pivot_longer(everything()) %>%
    separate(name, into = c('Stat', 'type'), sep = '\\.') %>%
    pivot_wider(names_from = type, values_from = value)
  
  p <-
    ggplot(d, aes(point, Stat, xmin=low, xmax=high)) +
    geom_pointrange(lwd=1.5, fatten = 2) +
    labs(caption = paste0('Proportion mediated: ', signif(100*obj$n.avg, 2), '%'), x='', y='') +
    theme_minimal(15)
  
  if (all(c(d$low, d$high) < 0)) p <- p + xlim(c(NA, 0))
  if (all(c(d$low, d$high) > 0)) p <- p + xlim(c(0, NA))
  return(p)
}

codebib <- function(){
  
  addkey <- function(x){
    x$key <- paste0(x$author$family[[1]], x$year[[1]], sample(letters,1), sample(letters,1))
    return(x)
  }
  
  d <- renv::dependencies() %>% dplyr::filter(Package != 'R')
  map(d$Package, ~toBibtex(addkey(citation(.x))))
}

#' @title ggxtabs
#' @description Mosaic ggplot of xtabs object
#' @param xtab xtabs object of exactly two categorical variables
#' @param cell_counts Display cell counts, Default: F
#' @param stats Display Chisq stats in caption, Default: F
#' @param viridis_option Viridis color palette, Default: 'G'
#' @param border_color Color of border between cells, Default: 'white'
#' @return ggplot
#' @details Takes a table of exactly two categorical variables produced by xtabs and produces a ggplot mosaic plot
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  ggxtabs(xtabs(~color+clarity, diamonds), viridis_option = 'H')
#'  }
#' }
#' @import dplyr tidyr tibble ggplot2
#' @rdname ggxtabs
#' @export 
ggxtabs <- function(xtab, cell_counts = F, stats = F, viridis_option = 'G', border_color = 'white'){
  if (!"xtabs" %in% class(xtab)) stop('xtab must be an xtabs object')
  d <- as_tibble(xtab)
  if (ncol(d) != 3) stop('Two categorical variables only')
  
  # put levels in order
  d[[1]] <- factor(d[[1]], levels = unique(d[[1]]))
  d[[2]] <- factor(d[[2]], levels = unique(d[[2]]))
  
  nms <- names(d)
  d2 <-
    d %>% 
    group_by(across(2)) %>% 
    summarise(sum = sum(n)) %>% 
    mutate(
      xmax = cumsum(sum),
      xmin = c(0, xmax[1:length(xmax)-1]),
      xmid = (xmin + xmax)/2
    )
  d <- 
    left_join(d, d2) %>% 
    group_by(across(2)) %>% 
    mutate(
      ymax = cumsum(n)/sum(n),
      ymin = c(0, ymax[1:length(ymax)-1]),
      ymid = (ymin + ymax)/2
    )
  
  d3 <- unique(d[c(2, 7)])
  
  # Dealing with ggplot weirdness (or my ignorance)
  d3$xmin = 0
  d3$xmax = 0
  d3$ymin = 0
  d3$ymax = 0
  
  p <- 
    ggplot(d, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) + 
    geom_rect(aes_string(fill=nms[1]), color=border_color) 
  
  if (cell_counts){
    p <- p + geom_label(aes(x=xmid, y = ymid, label = n), label.size=0)
  }  
  
  if (stats){
    st <- summary(xtab)
    stat_str <- glue::glue("N={st$n.cases}, χ²={signif(st$statistic, 3)}, df={st$parameter}, p={signif(st$p.value, 2)}")
    p <- p + labs(caption=stat_str)
  }
  
  p +
    geom_text(data=d3, aes_string(x='xmid', label=nms[2]), y = 1.05) +
    scale_fill_viridis_d(option = viridis_option) +
    scale_y_continuous(limits = c(0, 1.05)) +
    guides(fill = guide_legend(reverse = T)) +
    coord_cartesian(clip = 'off') +
    xlab(nms[2]) +
    theme_void() +
    theme(axis.title.x = element_text())
}

#' @title model_stats
#' @description Return a named list of estimates, standard errors, and p-values for one or more regression models.
#' @param ... Either a singal regression model object, or multiple named models, e.g., m1 = m1, m2 = m2, ...
#' @return A named list of model parameters and statistics
#' @details Uses the tidy functions from broom and broom.mixed to return a named list of model stats that can be easily used inline in rmarkdown documents.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  m1 <- lm(mpg ~ hp, mtcars)
#'  m2 <- lm(mpg ~ wt, mtcars)
#'  x <- model_stats(m1 = m1, m2 = m2)
#'  
#'  x$m1$hp$estimate
#'  x$m1$hp$str # formatted beta coefficient with 95% CI.
#'  }
#' }
#' @seealso 
#'  \code{\link[broom.mixed]{reexports}}
#'  \code{\link[glue]{glue}}
#' @rdname model_stats
#' @export 
#' @importFrom broom.mixed tidy
#' @importFrom glue glue
model_stats <- function(...){
  tdy <- function(m){
    broom.mixed::tidy(m, conf.int = T) %>% 
      mutate(
        str = glue::glue("$\\beta={signif(estimate, 2)}$ ({signif(conf.low, 2)}, {signif(conf.high, 2)})")
      ) %>% 
      split(.$term)
  }
  models <- list(...)
  if (length(models) == 1) return(tdy(models[[1]]))
  nms <- names(models)
  if (is.null(nms) | "" %in% nms) stop('all models must be named, e.g., m1 = m, ...')
  map(models, tdy)
}

#' @title twittersize
#' @description Add borders to an image to force an aspect ratio for twitter post, and write it out in png format. Does not resize original image.
#' @param path Path to image file
#' @param out_prefix Prepend to name of output image, Default = "twitter"
#' @param border_color Border color, Default = "white"
#' @param aspect_ratio Desired aspect ratio, Default = "16:9" (single image in tweet). Use "7:8" if tweeting two images.
#' @return path to output image
#' @export
#' @examples
#' \dontrun{
#' if(interactive()){
#'  twittersize('~/Desktop/image.png')
#'  }
#' }
twittersize <- function(path, aspect_ratio = "16:9", out_prefix = "twitter ", border_color = 'white'){
  ar <- as.numeric(strsplit(aspect_ratio, ':')[[1]])
  width <- ar[1]
  height <- ar[2]
  img <- magick::image_read(path)
  info <- magick::image_info(img)
  if (info$width > width*info$height/height){
    border_height <- round(((height*info$width/width) - info$height)/2)
    geometry <- paste0("0x", border_height)
  } else {
    border_width <- round(((width*info$height/height) - info$width)/2)
    geometry <- paste0(border_width, 'x0')
  }
  img <- magick::image_border(img, geometry = geometry, color = border_color)
  fn <- paste0(out_prefix, basename(tools::file_path_sans_ext(path)), '.png')
  path_out <- file.path(dirname(path), fn)
  magick::image_write(img, path = path_out, format = "png")
}