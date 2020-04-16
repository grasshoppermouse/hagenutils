
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
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @return Cohen's d
#' @examples
#' cohen_d(x, y)
#' @export 
cohen_d = function(x, y){
  x <- na.omit(x)
  y <- na.omit(y)
  nx = length(x)
  ny = length(y)
  pooled_sd = sqrt(((nx-1)*sd(x)^2 + (ny-1)*sd(y)^2)/(nx + ny))
  return((mean(x)-mean(y))/pooled_sd)
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
ggdotchart <- function(v){
  nms <- names(v)
  tibble(
    x = as.numeric(v),
    y = factor(nms, levels = nms[order(v)])
  ) %>%
    ggplot(aes(x, y)) +
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
      if (n != 2) stop('custom binary palette provides exactly 2 colors')
      cols <- viridisLite::magma(11)[c(4,8)]
      if (direction >= 0) cols else rev(cols)
    }, 
    ...
  )
}

#' @title ggemmeans
#' @description ggplot of estimated marginal means object, with categorical spec variable sorted by estimate
#' @param em An emmeans object
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
ggemmeans <- function(em){
  emsum <- emmeans::summary.emmGrid(em)
  estName <- attr(emsum, 'estName')
  clNames <- attr(emsum, 'clNames')
  var <- attr(emsum, 'pri.vars')[1]
  if (is.character(emsum[var])){
    emsum[var] <- forcats::fct_reorder(emsum[var], emsum[estName])
  }
  ggplot(emsum, aes_string(estName, var, xmin = clNames[1], xmax = clNames[2])) +
    geom_errorbarh(height = 0, lwd = 2.5, alpha = 0.3) +
    geom_point() +
    theme_minimal(15)
}
