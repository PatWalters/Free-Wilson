#!/usr/bin/env python

import sys
import re

import numpy as np
import pandas as pd
from pyfancy import pyfancy
import joblib
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score


def read_and_combine_dataframes(descriptor_file_name, activity_file_name, log_scale=False):
    """
    Read dataframes and merge
    @param descriptor_file_name: dataframe with R-group descriptor_values
    @param activity_file_name: csv file with activity data, format is "Name,Act" where Act is activity
    @param log_scale: log transform the data
    @return: combined dataframe, activity dataframe, descriptor_names
    """
    desc_df = pd.read_csv(descriptor_file_name)
    desc_names = desc_df.columns[1:]
    act_df = pd.read_csv(activity_file_name)
    act_df.dropna(inplace=True)
    if "Act" not in act_df.columns:
        error_msg = """Error:\nThe activity file %s does not have a column called "Act"\nPlease fix %s""" % (
            activity_file_name, activity_file_name)
        pyfancy.pyfancy().yellow().bold().add(error_msg).output()
        sys.exit(0)
    if log_scale:
        act_df['Act'] = -np.log10(act_df['Act'] * 1e-6)
        print("Transformed data to log scale")
    combo_df = pd.merge(act_df, desc_df, on="Name")
    return combo_df, act_df, desc_names


def run_regression(x_val, y_val, prefix):
    """
    Do ridge regression on descriptors vs activity
    @param x_val: descriptors
    @param y_val: activity values
    @param prefix: prefix for the output file names
    @return: predicted values, regression model
    """
    pyfancy.pyfancy().red().bold().add("Running Ridge Regression").output()
    lm = Ridge()
    lm.fit(x_val, y_val)
    # write the model out for future predicitons
    model_file_name = prefix + "_lm.pkl"
    joblib.dump(lm, model_file_name)
    print(f"wrote model to {model_file_name}")
    # get the predicted values
    pred_val = lm.predict(x_val)
    return pred_val, lm


def write_output(prediction_list, lm, combo_df, act_df, desc_names,
                 prefix):
    """
    Write the output files
    @param prediction_list: list of predicted values
    @param lm: model created in run_regression
    @param combo_df: dataframe with descriptors and experimental activities
    @param act_df: dataframe with activities
    @param desc_names: descriptor names
    @param prefix:
    @return: No return value
    """
    # get the y values again - kind of wasteful
    y = combo_df.values[0::, 1]
    # create a dictionary to hold the number of occurrences of each R-group
    r_group_count_dict = dict(combo_df.sum())
    combo_df.to_csv("combo.csv", index=False)
    # add the predictions to the activity dataframe
    act_df["Pred"] = prediction_list
    # print out some stats
    print("Intercept = %.2f" % lm.intercept_)
    print("R**2 = %.2f" % r2_score(y, prediction_list))
    # create a dataframe from regression coefficients
    coefficient_df = pd.DataFrame(list(zip(desc_names, lm.coef_)), columns=["R-Group SMILES", "Coefficient"])
    # create a list of R-group numbers and add to the coefficients dataframe
    feature_list = list(coefficient_df['R-Group SMILES'])
    rg_re = re.compile("\[\*:([0-9]+)\]")
    coefficient_df['R-group'] = [int(rg_re.search(x).group(1)) for x in feature_list]
    # add R-group counts to the coefficients dataframe
    r_group_names = list(coefficient_df['R-Group SMILES'])
    r_group_counts = [r_group_count_dict[x] for x in r_group_names]
    coefficient_df['Count'] = r_group_counts
    # write the comparison and coefficient output files
    comparison_csv = prefix + "_comparison.csv"
    act_df.to_csv(comparison_csv, index=False, float_format="%.2f")
    print(f"wrote comparison to {comparison_csv}")
    coefficient_csv = prefix + "_coefficients.csv"
    coefficient_df.to_csv(coefficient_csv, index=False, float_format="%.3f")
    print(f"wrote coefficients to {coefficient_csv}")


def free_wilson_regression(descriptor_file_name, activity_file_name, prefix, log_scale):
    """
    Driver function
    @param descriptor_file_name: file with R-group descriptor values
    @param activity_file_name: file with activity data, format is "Name,Act" where "Act" is the activity data
    @param prefix: prefix for output file names
    @param log_scale: log transform the data
    @return:
    """
    # read the data from the input files
    combo_df, act_df, desc_names = read_and_combine_dataframes(descriptor_file_name, activity_file_name, log_scale)
    # remove any rows in the activity file that are not in the descriptor file
    act_df = act_df[act_df['Name'].isin(combo_df['Name'])]
    # extract x and y variable from the data frame
    y = combo_df.values[0::, 1]
    x = combo_df.values[0::, 2::]
    # do the regression and write the output files
    pred_vals, lm = run_regression(x, y, prefix)
    write_output(pred_vals, lm, combo_df, act_df, desc_names, prefix)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: %s descriptors.csv activities.csv prefix" % sys.argv[0])
        sys.exit(0)
    free_wilson_regression(sys.argv[1], sys.argv[2], sys.argv[3], False)
