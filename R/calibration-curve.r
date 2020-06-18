#' Calibration Curve
#'
#' `calibration_curve()` computes the true and predicted probabilities for a
#' calibration curve.
#'
#' The function takes on inputs coming from a binary classifier.
#'
#' Calibration curve is also known as reliability diagram. This function is
#' named as so to be akin to scikit-learn's calibration_curve method.
#'
#' Quotes from Niculescu-Mizil & Caruana (2005) with minor modifications: First,
#' the predicted values (probabilities) is discretized into ten bins (default,
#' can be changed). Cases with predicted values between 0 and 0.1 fall in the
#' first bin, between 0.1 and 0.2 in the second bin, etc. For each bin, the mean
#' predicted value is plotted against the true fraction of positive cases.
#'
#' There is a [ggplot2::autoplot()] method for quickly visualising the curve.
#' his works for binary and multiclass output, and also works with grouped data
#' (i.e. from resamples).
#'
#' @param truth The column identifier for the true class results (that is a
#'     `factor`). This should be an unquoted column name although this argument
#'     is passed by expression and supports quasiquotation (you can unquote
#'     column names). For `_vec()` functions, a factor vector.
#' @param estimate The column identifier for the predicted results (that is also
#'     `numeric`). As with `truth` this can be specified different ways but the
#'     primary method is to use an unquoted variable name. For `_vec()`
#'     functions, a `numeric` vector.
#' @param n_bins Number of bins to discretize the `[0,1]` interval. Default is
#'     10.
#' @param scale_estimate A `logical` value indicating whether `estimate` should
#'     be normalised into the `[0,1]` interval.
#' @param discretise_strategy Strategy used to define the widths of the bins
#'     which is either 'uniform' (default) or 'quantile'. If 'uniform', the bins
#'     have idential widths. If 'quantile', the bins have the same number of
#'     samples.
#' @param na_rm A `logical` value indicating whether NA values should be
#'     stripped before the computation proceeds.
#' @param object The `clbr_df` data frame returned from `calibration_curve()`
#'
#' @family curve metrics
#' @template multiclass-curve
#' @template event_first
#' 
#'
#' @return
#' A tibble with `clbr_df` or `clbr_grouped_df` having columns `.frac_positive`
#' and `.mean_predicted`
#'
#' @author An Chu
#' @examples
#'
#' \dontrun{
#'
#' library(dplyr)
#' library(ggplot2)
#'
#' data("two_class_example", package = "yardstick")
#'
#' two_class_example %>%
#'     calibration_curve(truth, Class1)
#' 
#' two_class_example %>%
#'    calibration_curve(truth, Class1) %>%
#'    autoplot()
#'
#' }
#'
#' @export
calibration_curve = function(data, ...) {
    UseMethod("calibration_curve")
}



#' @export
#' @rdname clbr_curve
#' @import yardstick
#' @import ggplot2
#' @importFrom rlang enquo eval_tidy enquos `!!!`
#' @importFrom dplyr do is_grouped_df as_tibble
#' @importFrom stats relevel
calibration_curve.data.frame = function(data,
                                        truth,
                                        ...,
                                        n_bins = 10L,
                                        scale_estimate = FALSE,
                                        discretise_strategy = c("uniform", "quantile"),
                                        na_rm = TRUE) {

    estimate = yardstick::dots_to_estimate(data, !!!rlang::enquos(...))
    truth = rlang::enquo(truth)
    
    yardstick:::validate_not_missing(truth)

    # Explicit handling of length 1 character vectors as column names
    truth = yardstick:::handle_chr_names(truth, colnames(data))
    estimate = yardstick:::handle_chr_names(estimate, colnames(data))

    discretise_strategy = match.arg(discretise_strategy)

    res = dplyr::do(
                     data,
                     calibration_curve_vec(
                         truth = rlang::eval_tidy(truth, data = .),
                         estimate = rlang::eval_tidy(estimate, data = .),
                         n_bins = n_bins,
                         scale_estimate = scale_estimate,
                         discretise_strategy = discretise_strategy,
                         na_rm = na_rm
                     )
                 )

    if (dplyr::is_grouped_df(res)) {
        class(res) = c("clbr_grouped_df", "clbr_df", class(res))
    } else {
        class(res) = c("clbr_df", class(res))
    }

    res
}



calibration_curve_vec = function(truth,
                                 estimate,
                                 n_bins = 10L,
                                 scale_estimate = FALSE,
                                 discretise_strategy = c("uniform", "quantile"),
                                 na_rm = TRUE,
                                 ...) {

    estimator = yardstick::finalize_estimator(truth, metric_class = "clbr_curve")

    calibration_curve_impl = function(truth,
                                      estimate,
                                      n_bins,
                                      scale_estimate,
                                      discretise_strategy) {
        
        calibration_curve_estimator_impl(truth, estimate, n_bins, scale_estimate, discretise_strategy, estimator)

    }

    yardstick::metric_vec_template(
        metric_impl = calibration_curve_impl,
        truth = truth,
        estimate = estimate,
        na_rm = na_rm,
        estimator = estimator,
        cls = c("factor", "numeric"),
        ...,
        n_bins = n_bins,
        scale_estimate = scale_estimate,
        discretise_strategy = discretise_strategy
    )

}


calibration_curve_estimator_impl = function(truth,
                                            estimate,
                                            n_bins,
                                            scale_estimate,
                                            discretise_strategy,
                                            estimator) {
    
    if (yardstick:::is_binary(estimator)) {
        calibration_curve_binary(truth, estimate, n_bins, scale_estimate, discretise_strategy)
    } else {
        calibration_curve_multiclass(truth, estimate, n_bins, scale_estimate, discretise_strategy)
    }
    
}

calibration_curve_binary = function(truth,
                                    estimate,
                                    n_bins,
                                    scale_estimate,
                                    discretise_strategy) {


    if (!getOption("yardstick.event_first", default = TRUE)) {
        lvls = levels(truth)
        truth = stats::relevel(truth, lvls[2L])
    }

    truth = as.integer(truth)
    truth = ifelse(truth == 1L, 1L, 0L)

    res = calibration_curve_binary_impl(truth, estimate, n_bins, scale_estimate, discretise_strategy)
    res = dplyr::as_tibble(res)

    res
}

calibration_curve_multiclass = function(truth,
                                        estimate,
                                        n_bins,
                                        scale_estimate,
                                        discretise_strategy) {

    yardstick:::one_vs_all_with_level(calibration_curve_binary, truth, estimate, n_bins, scale_estimate, discretise_strategy)

}

calibration_curve_binary_impl = function(truth,
                                         estimate,
                                         n_bins,
                                         scale_estimate,
                                         discretise_strategy) {

    ## normalise predicted values into interval [0,1]
    if (isTRUE(scale_estimate)) {
        estimate = (estimate - min(estimate)) / (max(estimate) - min(estimate))
    } else if (min(estimate) < 0 || max(estimate) > 1) {
        stop("`estimate` has values outsite interval [0,1] and `scale_estimate` is set to FALSE.", call. = FALSE)
    }

    ## choosing bin strategy
    if (discretise_strategy == "uniform") {
        bins = ggplot2::cut_number(estimate, n_bins)
    } else if (discretise_strategy == "quantile") {
        bins = ggplot2::cut_interval(estimate, n_bins)
    } else {
        stop("Invalid value of `discretise_strategy` which must be either 'uniform' or 'quantile'.", call. = FALSE)
    }

    bins_lvls = levels(bins)
    frac_pos = numeric(n_bins)
    mean_pred = numeric(n_bins)

    for (i in seq_len(n_bins)) {
        match_idx = which(bins == bins_lvls[i])
        frac_pos[i] = sum(truth[match_idx]) / length(truth[match_idx])
        mean_pred[i] = mean(estimate[match_idx])
    }

    out = data.frame(.mean_predicted = mean_pred, .frac_positive = frac_pos)
    out
}


#' @rdname clbr_curve
autoplot.clbr_df = function(object, ...) {

    `%+%` = ggplot2::`%+%`

    ## base chart
    clbr_chart = ggplot2::ggplot(data = object)

    if (inherits(object, "clbr_grouped_df")) {
        grps = dplyr::groups(object)
        grps_chr = paste0(dplyr::group_vars(object), collapse = "_")
        interact_expr = list(color = rlang::expr(interaction(!!!grps, sep = "_")))
        clbr_chart = clbr_chart %+%
            ggplot2::labs(color = grps_chr)
    } else {
        interact_expr = list()
    }

    aes_spliced = ggplot2::aes(x = .mean_predicted,
                               y = .frac_positive,
                               !!!interact_expr)

    ## build the graph
    clbr_chart = clbr_chart %+%
        ggplot2::geom_path(mapping = aes_spliced) %+%
        ggplot2::geom_point(mapping = aes_spliced) %+%
        ggplot2::geom_abline(linetype = 3) %+%
        ggplot2::coord_equal() %+%
        ggplot2::theme_bw() %+%
        ggplot2::labs(x = "Mean Predicted Values", y = "Fraction of Positives")

    ## for multiclass case
    if (".level" %in% colnames(object)) {
        clbr_chart = clbr_chart %+%
            ggplot2::facet_wrap( ~ .level)
    }

    clbr_chart
}
