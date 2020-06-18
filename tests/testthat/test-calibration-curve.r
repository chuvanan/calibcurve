

context("Test against sklearn's calibration_curve()")

read_pydata = function(py_path) {
    ## depending on if testing one at a time or running
    ## CTRL+SHIFT+T
    if(interactive()) {
        ## interactive + local use. devtools::test() interactively
        if (Sys.getenv("NOT_CRAN") == "true") {
            path = "../sklearn-compare/"
        } else {
            path = "tests/sklearn-compare/"
        }
    } else {
        path = "../sklearn-compare/"
    }
    read.csv(paste0(path, py_path))
}


py_unf_clbr = read_pydata("unf_df.csv")
py_unf_clbr = py_unf_clbr[, c(2, 1)]
py_qtl_clbr = read_pydata("qtl_df.csv")
py_qtl_clbr = py_qtl_clbr[, c(2, 1)]

data(two_class_example, package = "yardstick")

r_unf_clbr = calibration_curve(two_class_example, truth, Class1, discretise_strategy = "uniform")
r_qtl_clbr = calibration_curve(two_class_example, truth, Class1, discretise_strategy = "quantile")

expect_equal(r_unf_clbr, py_unf_clbr)
expect_equal(r_qtl_clbr, py_qtl_clbr)
