
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
