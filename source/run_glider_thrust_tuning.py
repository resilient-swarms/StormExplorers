import pathlib
import sys
module_dir = pathlib.Path(__file__).parent.resolve()
root_dir = module_dir.parent
model_dir = root_dir.joinpath("model")
asvlite_wrapper_dir = root_dir.joinpath("dependency", "ASVLite", "wrapper", "cython")
sys.path.insert(0, str(asvlite_wrapper_dir))

import math
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from xgboost import XGBRegressor
from autosklearn.regression import AutoSklearnRegressor
from joblib import dump

import glider_thrust_factor

# Tune thrust
tuning = glider_thrust_factor._Thrust_tuning()
tuning.tune_wg_thrust()

# Create the models
# Load the training data
training_data_file = "../results/glider_thrust/tuning/tuning_factors_for_training.csv"
cols_to_use = ["Timestamp(UTC)",
               "distance(m)",
               "delta_T(s)",
               "delta_T_simulated(s)",
               "speed(knots)",
               "speed_simulated(knots)",
               "simulated_wave_height(m)",
               "current(knots)",
               "relative_current_direction",
               "tuning_factor",
               "error_msg"]
df = pd.read_csv(training_data_file, usecols = cols_to_use)
df["relative_current_direction(abs)"] = df["relative_current_direction"].abs()
df = df[df["tuning_factor"].notna()]
df.reset_index(drop=True, inplace=True)

# Group segments based on simuated wave heights
for i in range(len(df)):
    df.loc[i, "tuning_group"] = int(math.ceil(df.loc[i, "simulated_wave_height(m)"]))

# Remove outliers and replot
q_1  = df["tuning_factor"].quantile(0.01)
q_99 = df["tuning_factor"].quantile(0.99)
df = df[(df["tuning_factor"] < q_99) & (df["tuning_factor"] > q_1)]

# Test and training data
X = df[["simulated_wave_height(m)", "current(knots)", "relative_current_direction(abs)"]].values
y = df["tuning_factor"].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

# Model 1 - Linear regression based only on wave height
model_1 = LinearRegression().fit(np.array([row[0] for row in X_train]).reshape(-1, 1), np.array(y_train))
# Model 2 - Linear regression using wave height, current speed and direction
model_2 = LinearRegression().fit(X_train, y_train)
# Model 3 - XGBoost
model_3 = XGBRegressor().fit(X_train, y_train)
# Model 4 - AutoML
model_4 = AutoSklearnRegressor().fit(X_train, y_train)
# Model 5 - Linear regression with uniform population
# Create a new training data with a max of 200 instances from each group
max_count = 200
df_new = pd.DataFrame()
for group_name, group in df.groupby("tuning_group"):
    if len(group) > max_count:
        df_new = pd.concat([df_new, group.sample(max_count)])
    else:
        df_new = pd.concat([df_new, group])
X_train = df_new[["simulated_wave_height(m)", "current(knots)", "relative_current_direction(abs)"]].values
y_train = df_new["tuning_factor"].values
model_5 = LinearRegression().fit(np.array([row[0] for row in X_train]).reshape(-1, 1), np.array(y_train))
# Model 6 - Linear regression with uniform population and taking wave and current
model_6 = LinearRegression().fit(X_train, y_train)

# Save models to file
file_name = model_dir.joinpath("thrust_tuning_lin_reg_1.joblib")
print("saving {}".format(file_name))
dump(model_1, str(file_name))
file_name = model_dir.joinpath("thrust_tuning_lin_reg_2.joblib")
print("saving {}".format(file_name))
dump(model_2, str(file_name))
file_name = model_dir.joinpath("thrust_tuning_xgboost.joblib")
print("saving {}".format(file_name))
dump(model_3, str(file_name))
file_name = model_dir.joinpath("thrust_tuning_automl.joblib")
print("saving {}".format(file_name))
dump(model_4, str(file_name))
file_name = model_dir.joinpath("thrust_tuning_lin_reg_3.joblib")
print("saving {}".format(file_name))
dump(model_5, str(file_name))
file_name = model_dir.joinpath("thrust_tuning_lin_reg_4.joblib")
print("saving {}".format(file_name))
dump(model_6, str(file_name))

# Reload the glider_thrust_factor module and generate results for Benjamin
import importlib
importlib.reload(glider_thrust_factor)
tuning = glider_thrust_factor._Thrust_tuning("HPC")
tuning.run_benjamin()