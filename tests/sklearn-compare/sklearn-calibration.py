#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.calibration import calibration_curve

two_class_example = pd.read_csv("two_class_example.csv")

y_truth = np.where(two_class_example["truth"] == "Class1", 1, 0)
y_prob = two_class_example["Class1"]

unf_frac_positive, unf_mean_predicted = calibration_curve(y_truth, y_prob, n_bins = 10, strategy="uniform")
qtl_frac_positive, qtl_mean_predicted = calibration_curve(y_truth, y_prob, n_bins = 10, strategy="quantile")

unf_df = pd.DataFrame({"unf_frac_pos":unf_frac_positive, "unf_mean_pre":unf_mean_predicted})
unf_df.to_csv("unf_df.csv", index=False)

qtl_df = pd.DataFrame({"qtl_frac_pos":qtl_frac_positive, "qtl_mean_pre":qtl_mean_predicted})
qtl_df.to_csv("qtl_df.csv", index=False)
