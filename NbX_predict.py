import glob
import numpy as np
import pandas as pd
import joblib

# Load
feature_df = pd.read_csv("NbX_feature.csv")
features = np.array(feature_df.loc[:,"epitope_positive_count":"chainH_MSWHIM3.1"])

model_list = glob.glob("./model/model*")

model_list.sort(reverse=True)

for model in model_list:

  # Predict
  xgb = joblib.load(model)
  
  pred_proba = xgb.predict_proba(features)[:,1]

  # Assign
  model_name = model.split("/")[2]

  feature_df.insert(1, model_name + "_predicted_CAPRI_binary_proba", pred_proba)

# Assign (mean)
proba_col = [col for col in feature_df.columns if '_predicted_CAPRI_binary_proba' in col]
mean_proba = feature_df[proba_col].mean(axis=1)
feature_df.insert(1, "mean_predicted_CAPRI_binary_proba", mean_proba)

# Save
feature_df.to_csv("NbX_prediction.csv")
