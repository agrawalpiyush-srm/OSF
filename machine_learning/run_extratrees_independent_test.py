import joblib
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

# ---------------- PATHS ----------------
MODEL_PATH = "ExtraTrees_FINAL_MODEL.joblib"
FEATURE_PATH = "FEATURE_NAMES.joblib"
TEST_DATA_PATH = "independent_test_features.csv"
OUTPUT_PATH = "ExtraTrees_independent_predictions.csv"

# ---------------- LOAD MODEL ----------------
print("Loading model and features...")
model = joblib.load(MODEL_PATH)
feature_names = joblib.load(FEATURE_PATH)

# ---------------- LOAD TEST DATA ----------------
print("Loading independent dataset...")
test_df = pd.read_csv(TEST_DATA_PATH)

# ---------------- FEATURE CHECK ----------------
missing_features = set(feature_names) - set(test_df.columns)
if missing_features:
    raise ValueError(f"Missing features in test data: {missing_features}")

X_test = test_df[feature_names]

# ---------------- PREDICTION ----------------
print("Running predictions...")
test_proba = model.predict_proba(X_test)[:, 1]
test_pred = model.predict(X_test)

# ---------------- SAVE OUTPUT ----------------
output_df = test_df.copy()
output_df["Predicted_Label"] = test_pred
output_df["Predicted_Probability"] = test_proba

output_df.to_csv(OUTPUT_PATH, index=False)
print("Predictions saved to:", OUTPUT_PATH)

# ---------------- OPTIONAL EVALUATION ----------------
if "Label" in test_df.columns:
    print("Evaluating performance...")
    y_true = test_df["Label"]
    print("AUROC:", roc_auc_score(y_true, test_proba))
    print("AUPRC:", average_precision_score(y_true, test_proba))
else:
    print("No ground truth labels found. Skipping evaluation.")

