

###############################################################################################
###### Identify discriminatory leaf chemicals which separates stressed and control condition ##
###############################################################################################
import pandas as pd
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler, LabelEncoder
from sklearn.model_selection import GridSearchCV, cross_val_score, StratifiedKFold
from sklearn.metrics import mean_squared_error, confusion_matrix, ConfusionMatrixDisplay
from sklearn.base import BaseEstimator, TransformerMixin
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
import plotly.express as px
import plotly.graph_objs as go

# 1. Initial Setup
# ----------------
# Font sizes
font_sizes = {
    'title': 16,
    'label': 14,
    'tick': 12,
    'legend': 12,
}

# 2. Data Loading and Preprocessing
# ---------------------------------
# Load the leaf tissue dataset
npl_data = pd.read_csv(r'C:\Users\ms\Desktop\data_chem_3_10\data\data\n_p_l.csv')

# Filter for leaf tissue
npl_data = npl_data[npl_data['Tissue.type'] == 'L']

# Extract metabolite data and treatment labels
metabolite_columns = [col for col in npl_data.columns if 'Cluster' in col]
X_leaf = npl_data[metabolite_columns]
Y_leaf = npl_data['Treatment']

# Ensure all data in X is numeric and handle missing values
X_leaf = X_leaf.apply(pd.to_numeric, errors='coerce')
X_leaf.fillna(X_leaf.mean(), inplace=True)

# Encode the Treatment labels to numeric values
label_encoder = LabelEncoder()
Y_leaf_encoded = label_encoder.fit_transform(Y_leaf)

# Display some information about the data
print(f"Number of samples: {X_leaf.shape[0]}")
print(f"Number of metabolites: {X_leaf.shape[1]}")
print(f"Treatment classes: {np.unique(Y_leaf)}")

# 3. Feature Selection
# --------------------
# Fit an initial PLS model to calculate VIP scores
pls_initial = PLSRegression(n_components=5)  # Initial number of components
pls_initial.fit(X_leaf, Y_leaf_encoded)

# VIP Calculation
def calculate_vip(model, X):
    t = model.x_scores_
    w = model.x_weights_
    q = model.y_loadings_

    p, h = w.shape
    s = np.sum(t**2, axis=0) * np.sum(q**2, axis=0)
    total_s = np.sum(s)

    vips = np.sqrt(p * np.sum((w**2 * s), axis=1) / total_s)
    return vips

# Calculate VIP scores
vip_scores_leaf = calculate_vip(pls_initial, X_leaf)

# Filter for VIP > 1
vip_scores_leaf_df = pd.DataFrame({'Metabolite': X_leaf.columns, 'VIP_Score': vip_scores_leaf})
vip_scores_leaf_filtered = vip_scores_leaf_df[vip_scores_leaf_df['VIP_Score'] > 1]
X_leaf_filtered = X_leaf[vip_scores_leaf_filtered['Metabolite']]

# Ensure X_leaf_filtered retains feature names
X_leaf_filtered.columns = X_leaf_filtered.columns.astype(str)

# 4. Model Selection
# ------------------
# Define scalers and parameter grid
scalers = {
    'Standard': StandardScaler(),
    'MinMax': MinMaxScaler(),
    'Robust': RobustScaler()
}

param_grid = {
    'scaler': ['Standard', 'MinMax', 'Robust'],
    'n_components': range(1, min(11, X_leaf_filtered.shape[1], X_leaf_filtered.shape[0]))
}

# Custom Pipeline
class ScalerPLSPipeline(BaseEstimator, TransformerMixin):
    def __init__(self, scaler=None, n_components=2):
        self.scaler = scaler
        self.n_components = n_components
        self.pls = PLSRegression(n_components=self.n_components)

    def fit(self, X, y):
        self.scaler.fit(X)
        X_scaled = self.scaler.transform(X)
        self.pls.fit(X_scaled, y)
        return self

    def predict(self, X):
        X_scaled = self.scaler.transform(X)
        return self.pls.predict(X_scaled)
    
    def score(self, X, y):
        X_scaled = self.scaler.transform(X)
        return self.pls.score(X_scaled, y)

# Nested Cross-Validation
outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Grid Search with a custom model
best_mse = float('inf')
best_model = None
cv_scores = []

for scaler_name, scaler in scalers.items():
    for n_comp in param_grid['n_components']:
        pipeline = ScalerPLSPipeline(scaler=scaler, n_components=n_comp)
        scores = cross_val_score(pipeline, X_leaf_filtered, Y_leaf_encoded, cv=outer_cv, scoring='neg_mean_squared_error')
        mean_mse = -scores.mean()
        cv_scores.extend(-scores)  # Store all CV scores
        if mean_mse < best_mse:
            best_mse = mean_mse
            best_model = pipeline

# Calculate and print nested CV results
cv_mean = np.mean(cv_scores)
cv_std = np.std(cv_scores)
print(f"Nested Cross-Validation MSE: {cv_mean:.4f} (+/- {cv_std:.4f})")

# 5. Model Fitting and Evaluation
# -------------------------------
# Fitting the final best model
best_model.fit(X_leaf_filtered, Y_leaf_encoded)
y_pred = best_model.predict(X_leaf_filtered)
mse = mean_squared_error(Y_leaf_encoded, y_pred)
r2 = best_model.score(X_leaf_filtered, Y_leaf_encoded)

print(f"Mean Squared Error (Best Model): {mse:.4f}")
print(f"R-squared (Best Model): {r2:.4f}")
print(f"Best Mean Squared Error: {best_mse:.4f}")

# Get scores from the best model
X_leaf_scores = best_model.pls.transform(X_leaf_filtered)

# 6. Explained Variance Calculation
# ---------------------------------
# Calculate explained variance for all components
total_variance = np.var(X_leaf_filtered, axis=0).sum()
explained_variance_ratio = np.var(X_leaf_scores, axis=0) / total_variance
cumulative_explained_variance = np.cumsum(explained_variance_ratio)

# 7. Plotting
# -----------
# 7.1 PLS-DA Scores Plot
plt.figure(figsize=(8, 6))
colors = ['#78c679', '#238443']  # Shades of green to bluish-green

for treatment in np.unique(Y_leaf_encoded):
    idx = Y_leaf_encoded == treatment
    plt.scatter(X_leaf_scores[idx, 0], X_leaf_scores[idx, 1], 
                label=f'Treatment {label_encoder.inverse_transform([treatment])[0]}', 
                color=colors[treatment])

plt.xlabel(f'PLS Component 1 ({explained_variance_ratio[0]:.2f} explained variance)', fontsize=font_sizes['label'])
plt.ylabel(f'PLS Component 2 ({explained_variance_ratio[1]:.2f} explained variance)', fontsize=font_sizes['label'])
plt.title('PLS-DA Scores Plot (Control vs Stressed) - Leaf Tissue', fontsize=font_sizes['title'])
plt.legend(fontsize=font_sizes['legend'])
plt.xticks(fontsize=font_sizes['tick'])
plt.yticks(fontsize=font_sizes['tick'])
plt.tight_layout()
plt.show()

# 7.2 PLS-DA Loading Plot
loadings_leaf = pd.DataFrame(best_model.pls.x_loadings_, 
                             index=X_leaf_filtered.columns, 
                             columns=[f'Component_{i+1}' for i in range(best_model.n_components)])
loadings_leaf['Metabolite'] = loadings_leaf.index  # Add metabolite names for plotly

plt.figure(figsize=(10, 8))
plt.scatter(loadings_leaf['Component_1'], loadings_leaf['Component_2'], s=100, color='#238443')
plt.xlabel('PLS Component 1', fontsize=font_sizes['label'])
plt.ylabel('PLS Component 2', fontsize=font_sizes['label'])
plt.title('PLS-DA Loading Plot - Leaf Tissue', fontsize=font_sizes['title'])

for i, metabolite in enumerate(loadings_leaf.index):
    plt.text(loadings_leaf['Component_1'][i], loadings_leaf['Component_2'][i], metabolite, fontsize=8)

plt.tight_layout()
plt.show()

# 7.3 Confusion Matrix and Classification Metrics
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, ConfusionMatrixDisplay

# Convert continuous PLS-DA predictions to binary class predictions
y_pred_class = (y_pred > 0.5).astype(int)

# Calculate additional metrics
accuracy = accuracy_score(Y_leaf_encoded, y_pred_class)
f1 = f1_score(Y_leaf_encoded, y_pred_class, average='weighted')

# Generate confusion matrix
conf_matrix = confusion_matrix(Y_leaf_encoded, y_pred_class)

# Plot confusion matrix
plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Greens',
            xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
plt.title('Confusion Matrix (PLS-DA Model)', fontsize=font_sizes['title'])
plt.xlabel('Predicted', fontsize=font_sizes['label'])
plt.ylabel('Actual', fontsize=font_sizes['label'])

# Add text annotations for accuracy and F1-score
plt.text(0.5, -0.15, f"Accuracy: {accuracy:.4f}", ha='center', va='center', transform=plt.gca().transAxes, fontsize=font_sizes['label'])
plt.text(0.5, -0.25, f"F1-Score: {f1:.4f}", ha='center', va='center', transform=plt.gca().transAxes, fontsize=font_sizes['label'])

plt.tight_layout()
plt.show()

# Print classification metrics
print(f"Accuracy: {accuracy:.4f}")
print(f"F1-Score: {f1:.4f}")

# Calculate and print class-wise metrics
for class_label in np.unique(Y_leaf_encoded):
    class_accuracy = accuracy_score(Y_leaf_encoded == class_label, y_pred_class == class_label)
    class_f1 = f1_score(Y_leaf_encoded == class_label, y_pred_class == class_label, average='binary')
    print(f"Class {label_encoder.inverse_transform([class_label])[0]}:")
    print(f"  Accuracy: {class_accuracy:.4f}")
    print(f"  F1-Score: {class_f1:.4f}")

# 7.4 VIP Scores Plot
plt.figure(figsize=(10, 8))
vip_scores_leaf_filtered.sort_values('VIP_Score', ascending=True).plot(kind='barh', x='Metabolite', y='VIP_Score', color='#238443')
plt.title('VIP Scores of Selected Features', fontsize=font_sizes['title'])
plt.xlabel('VIP Score', fontsize=font_sizes['label'])
plt.ylabel('Metabolites', fontsize=font_sizes['label'])
plt.tight_layout()
plt.show()

# 7.5 Explained Variance Plot
plt.figure(figsize=(10, 6))
individual_color = '#66c2a5'  # Light green
cumulative_color = '#238b45'  # Dark green
plt.bar(range(1, len(explained_variance_ratio) + 1), explained_variance_ratio, 
        alpha=0.7, align='center', label='Individual explained variance', color=individual_color)
plt.step(range(1, len(cumulative_explained_variance) + 1), cumulative_explained_variance, 
         where='mid', label='Cumulative explained variance', color=cumulative_color, linewidth=2)
plt.ylabel('Explained variance ratio', fontsize=font_sizes['label'])
plt.xlabel('Principal components', fontsize=font_sizes['label'])
plt.legend(loc='best', fontsize=font_sizes['legend'])
plt.title('Explained Variance by Principal Components', fontsize=font_sizes['title'])
plt.xticks(fontsize=font_sizes['tick'])
plt.yticks(fontsize=font_sizes['tick'])
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()


# 7.6 Interactive PLS-DA Loading Plot
vip_threshold = 1.5
loadings_leaf['VIP_Score'] = vip_scores_leaf_filtered['VIP_Score'].values

# Custom color scale from light yellowish green to dark green, with more gradual transitions
custom_colorscale = [
    [0, "#077487"],    # Very light green (low VIP score)
    [0.2, "#15857d"],  # Light green
    [0.4, "#1c853f"],  # Bluish green
    [0.6, "#2a7a28"],  # Green-blue (mid score)
    [0.8, "#70b840"],  # Dark bluish green
    [1, "#d1f263"]     # Dark green (high VIP score)
]

# Create a more extreme size scaling
min_size = 5
max_size = 50
loadings_leaf['Marker_Size'] = min_size + (max_size - min_size) * (loadings_leaf['VIP_Score'] - loadings_leaf['VIP_Score'].min()) / (loadings_leaf['VIP_Score'].max() - loadings_leaf['VIP_Score'].min())

# Apply a power function to make the size difference more extreme
loadings_leaf['Marker_Size'] = loadings_leaf['Marker_Size'] ** 1.5

# Create the interactive scatter plot
fig = px.scatter(loadings_leaf, 
                 x='Component_1', y='Component_2', 
                 text=loadings_leaf.apply(lambda row: row['Metabolite'] if row['VIP_Score'] > vip_threshold else "", axis=1),
                 size='Marker_Size', color='VIP_Score',
                 size_max=max_size,  # Set the maximum size of the markers
                 title='Interactive PLS-DA Loading Plot - Leaf Tissue (VIP Colored)',
                 labels={
                     'Component_1': 'PLS Component 1',
                     'Component_2': 'PLS Component 2',
                     'VIP_Score': 'VIP Score'
                 },
                 hover_data={'Metabolite': True, 'VIP_Score': ':.2f'},
                 color_continuous_scale=custom_colorscale)

# Update trace appearance with grey text for metabolite names and opacity based on VIP score
fig.update_traces(
    textposition='top center', 
    textfont=dict(color='#274357'),
    marker=dict(opacity=0.6 + 0.4 * (loadings_leaf['VIP_Score'] - loadings_leaf['VIP_Score'].min()) / 
                (loadings_leaf['VIP_Score'].max() - loadings_leaf['VIP_Score'].min()))
)

# Update layout for a cleaner, more modern look
fig.update_layout(
    title_font_size=18,
    xaxis_title='PLS Component 1',
    yaxis_title='PLS Component 2',
    hoverlabel=dict(bgcolor="white", font_size=12, font_family="Rockwell"),
    font=dict(size=14),
    coloraxis_colorbar=dict(title="VIP Score"),
    plot_bgcolor='white',        # Set the plot background to white
    paper_bgcolor='white',       # Set the paper background to white
    showlegend=False             # Hide the legend as it's not needed anymore
)

# Display the plot
fig.show()

# Save the interactive plot as an HTML file
fig.write_html("interactive_loading_plot_leaf.html")

# 8. Save Results
# ---------------
# Save loadings
loadings_leaf.to_csv(r'C:\Users\ms\Desktop\data_chem_3_10\output\results\updated_PLS\pls_loadings_leaf.csv')

# Save VIP scores (all scores, not just filtered)
vip_scores_leaf_df.to_csv(r'C:\Users\ms\Desktop\data_chem_3_10\output\results\updated_PLS\vip_scores_all_leaf.csv', index=False)

# Save filtered VIP scores (VIP > 1)
vip_scores_leaf_filtered.to_csv(r'C:\Users\ms\Desktop\data_chem_3_10\output\results\updated_PLS\vip_scores_filtered_leaf.csv', index=False)

# Save explained variance data
pd.DataFrame({
    'Component': range(1, len(explained_variance_ratio) + 1),
    'Explained_Variance_Ratio': explained_variance_ratio,
    'Cumulative_Explained_Variance': cumulative_explained_variance
}).to_csv(r'C:\Users\ms\Desktop\data_chem_3_10\output\results\updated_PLS\explained_variance.csv', index=False)

print("Results saved successfully:")
print("1. PLS Loadings: pls_loadings_leaf.csv")
print("2. All VIP Scores: vip_scores_all_leaf.csv")
print("3. Filtered VIP Scores (VIP > 1): vip_scores_filtered_leaf.csv")
print("4. Explained Variance: explained_variance.csv")
print("5. Interactive Loading Plot: interactive_loading_plot.html")






###############################################################################################
###### Identify discriminatory root chemicals which separates stressed and control condition ##
###############################################################################################

import pandas as pd
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler, LabelEncoder
from sklearn.model_selection import GridSearchCV, cross_val_score, StratifiedKFold
from sklearn.metrics import mean_squared_error, confusion_matrix, ConfusionMatrixDisplay
from sklearn.base import BaseEstimator, TransformerMixin
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
import plotly.graph_objs as go

# 1. Initial Setup
# ----------------
# Font sizes
font_sizes = {
    'title': 16,
    'label': 14,
    'tick': 12,
    'legend': 12,
}

# 2. Data Loading and Preprocessing
# ---------------------------------
# Load the root tissue dataset
npr_data = pd.read_csv(r'C:\Users\ms\Desktop\data_chem_3_10\data\data\n_p_r.csv')

# Filter for root tissue
npr_data = npr_data[npr_data['Tissue.type'] == 'R']

# Extract metabolite data and treatment labels
metabolite_columns = [col for col in npr_data.columns if 'Cluster' in col]
X_root = npr_data[metabolite_columns]
Y_root = npr_data['Treatment']

# Ensure all data in X is numeric and handle missing values
X_root = X_root.apply(pd.to_numeric, errors='coerce')
X_root.fillna(X_root.mean(), inplace=True)

# Encode the Treatment labels to numeric values
label_encoder = LabelEncoder()
Y_root_encoded = label_encoder.fit_transform(Y_root)

# Display some information about the data
print(f"Number of samples: {X_root.shape[0]}")
print(f"Number of metabolites: {X_root.shape[1]}")
print(f"Treatment classes: {np.unique(Y_root)}")

# 3. Feature Selection
# --------------------
# Fit an initial PLS model to calculate VIP scores
pls_initial = PLSRegression(n_components=5)  # Initial number of components
pls_initial.fit(X_root, Y_root_encoded)

# VIP Calculation
def calculate_vip(model, X):
    t = model.x_scores_
    w = model.x_weights_
    q = model.y_loadings_

    p, h = w.shape
    s = np.sum(t**2, axis=0) * np.sum(q**2, axis=0)
    total_s = np.sum(s)

    vips = np.sqrt(p * np.sum((w**2 * s), axis=1) / total_s)
    return vips

# Calculate VIP scores
vip_scores_root = calculate_vip(pls_initial, X_root)

# Filter for VIP > 1
vip_scores_root_df = pd.DataFrame({'Metabolite': X_root.columns, 'VIP_Score': vip_scores_root})
vip_scores_root_filtered = vip_scores_root_df[vip_scores_root_df['VIP_Score'] > 1]
X_root_filtered = X_root[vip_scores_root_filtered['Metabolite']]

# Ensure X_root_filtered retains feature names
X_root_filtered.columns = X_root_filtered.columns.astype(str)

# 4. Model Selection
# ------------------
# Define scalers and parameter grid
scalers = {
    'Standard': StandardScaler(),
    'MinMax': MinMaxScaler(),
    'Robust': RobustScaler()
}

param_grid = {
    'scaler': ['Standard', 'MinMax', 'Robust'],
    'n_components': range(1, min(11, X_root_filtered.shape[1], X_root_filtered.shape[0]))
}

# Custom Pipeline
class ScalerPLSPipeline(BaseEstimator, TransformerMixin):
    def __init__(self, scaler=None, n_components=2):
        self.scaler = scaler
        self.n_components = n_components
        self.pls = PLSRegression(n_components=self.n_components)

    def fit(self, X, y):
        self.scaler.fit(X)
        X_scaled = self.scaler.transform(X)
        self.pls.fit(X_scaled, y)
        return self

    def predict(self, X):
        X_scaled = self.scaler.transform(X)
        return self.pls.predict(X_scaled)
    
    def score(self, X, y):
        X_scaled = self.scaler.transform(X)
        return self.pls.score(X_scaled, y)

# Nested Cross-Validation
outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Grid Search with a custom model
best_mse = float('inf')
best_model = None
cv_scores = []

for scaler_name, scaler in scalers.items():
    for n_comp in param_grid['n_components']:
        pipeline = ScalerPLSPipeline(scaler=scaler, n_components=n_comp)
        scores = cross_val_score(pipeline, X_root_filtered, Y_root_encoded, cv=outer_cv, scoring='neg_mean_squared_error')
        mean_mse = -scores.mean()
        cv_scores.extend(-scores)  # Store all CV scores
        if mean_mse < best_mse:
            best_mse = mean_mse
            best_model = pipeline

# Calculate and print nested CV results
cv_mean = np.mean(cv_scores)
cv_std = np.std(cv_scores)
print(f"Nested Cross-Validation MSE: {cv_mean:.4f} (+/- {cv_std:.4f})")

# 5. Model Fitting and Evaluation
# -------------------------------
# Fitting the final best model
best_model.fit(X_root_filtered, Y_root_encoded)
y_pred = best_model.predict(X_root_filtered)
mse = mean_squared_error(Y_root_encoded, y_pred)
r2 = best_model.score(X_root_filtered, Y_root_encoded)

print(f"Mean Squared Error (Best Model): {mse:.4f}")
print(f"R-squared (Best Model): {r2:.4f}")
print(f"Best Mean Squared Error: {best_mse:.4f}")

# Get scores from the best model
X_root_scores = best_model.pls.transform(X_root_filtered)

# 6. Explained Variance Calculation
# ---------------------------------
# Calculate explained variance for all components
total_variance = np.var(X_root_filtered, axis=0).sum()
explained_variance_ratio = np.var(X_root_scores, axis=0) / total_variance
cumulative_explained_variance = np.cumsum(explained_variance_ratio)

# 7. Plotting
# -----------
# 7.1 PLS-DA Scores Plot
plt.figure(figsize=(8, 6))
colors = ['#78c679', '#238443']  # Shades of green to bluish-green

for treatment in np.unique(Y_root_encoded):
    idx = Y_root_encoded == treatment
    plt.scatter(X_root_scores[idx, 0], X_root_scores[idx, 1], 
                label=f'Treatment {label_encoder.inverse_transform([treatment])[0]}', 
                color=colors[treatment])

plt.xlabel(f'PLS Component 1 ({explained_variance_ratio[0]:.2f} explained variance)', fontsize=font_sizes['label'])
plt.ylabel(f'PLS Component 2 ({explained_variance_ratio[1]:.2f} explained variance)', fontsize=font_sizes['label'])
plt.title('PLS-DA Scores Plot (Control vs Stressed) - Root Tissue', fontsize=font_sizes['title'])
plt.legend(fontsize=font_sizes['legend'])
plt.xticks(fontsize=font_sizes['tick'])
plt.yticks(fontsize=font_sizes['tick'])
plt.tight_layout()
plt.show()

# 7.2 PLS-DA Loading Plot
loadings_root = pd.DataFrame(best_model.pls.x_loadings_, 
                             index=X_root_filtered.columns, 
                             columns=[f'Component_{i+1}' for i in range(best_model.n_components)])
loadings_root['Metabolite'] = loadings_root.index  # Add metabolite names for plotly

plt.figure(figsize=(10, 8))
plt.scatter(loadings_root['Component_1'], loadings_root['Component_2'], s=100, color='#238443')
plt.xlabel('PLS Component 1', fontsize=font_sizes['label'])
plt.ylabel('PLS Component 2', fontsize=font_sizes['label'])
plt.title('PLS-DA Loading Plot - Root Tissue', fontsize=font_sizes['title'])

for i, metabolite in enumerate(loadings_root.index):
    plt.text(loadings_root['Component_1'][i], loadings_root['Component_2'][i], metabolite, fontsize=8)

plt.tight_layout()
plt.show()

# 7.3 Confusion Matrix and Classification Metrics
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, ConfusionMatrixDisplay

# Convert continuous PLS-DA predictions to binary class predictions
y_pred_class = (y_pred > 0.5).astype(int)

# Calculate additional metrics
accuracy = accuracy_score(Y_root_encoded, y_pred_class)
f1 = f1_score(Y_root_encoded, y_pred_class, average='weighted')

# Generate confusion matrix
conf_matrix = confusion_matrix(Y_root_encoded, y_pred_class)

# Plot confusion matrix
plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Greens',
            xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
plt.title('Confusion Matrix (PLS-DA Model) - Root Tissue', fontsize=font_sizes['title'])
plt.xlabel('Predicted', fontsize=font_sizes['label'])
plt.ylabel('Actual', fontsize=font_sizes['label'])

# Add text annotations for accuracy and F1-score
plt.text(0.5, -0.15, f"Accuracy: {accuracy:.4f}", ha='center', va='center', transform=plt.gca().transAxes, fontsize=font_sizes['label'])
plt.text(0.5, -0.25, f"F1-Score: {f1:.4f}", ha='center', va='center', transform=plt.gca().transAxes, fontsize=font_sizes['label'])

plt.tight_layout()
plt.show()

# Print classification metrics
print(f"Accuracy: {accuracy:.4f}")
print(f"F1-Score: {f1:.4f}")

# Calculate and print class-wise metrics
for class_label in np.unique(Y_root_encoded):
    class_accuracy = accuracy_score(Y_root_encoded == class_label, y_pred_class == class_label)
    class_f1 = f1_score(Y_root_encoded == class_label, y_pred_class == class_label, average='binary')
    print(f"Class {label_encoder.inverse_transform([class_label])[0]}:")
    print(f"  Accuracy: {class_accuracy:.4f}")
    print(f"  F1-Score: {class_f1:.4f}")

# 7.4 VIP Scores Plot
plt.figure(figsize=(10, 8))
vip_scores_root_filtered.sort_values('VIP_Score', ascending=True).plot(kind='barh', x='Metabolite', y='VIP_Score', color='#238443')
plt.title('VIP Scores of Selected Features - Root Tissue', fontsize=font_sizes['title'])
plt.xlabel('VIP Score', fontsize=font_sizes['label'])
plt.ylabel('Metabolites', fontsize=font_sizes['label'])
plt.tight_layout()
plt.show()

# 7.5 Explained Variance Plot
plt.figure(figsize=(10, 6))
individual_color = '#66c2a5'  # Light green
cumulative_color = '#238b45'  # Dark green
plt.bar(range(1, len(explained_variance_ratio) + 1), explained_variance_ratio, 
        alpha=0.7, align='center', label='Individual explained variance', color=individual_color)
plt.step(range(1, len(cumulative_explained_variance) + 1), cumulative_explained_variance, 
         where='mid', label='Cumulative explained variance', color=cumulative_color, linewidth=2)
plt.ylabel('Explained variance ratio', fontsize=font_sizes['label'])
plt.xlabel('Principal components', fontsize=font_sizes['label'])
plt.legend(loc='best', fontsize=font_sizes['legend'])
plt.title('Explained Variance by Principal Components - Root Tissue', fontsize=font_sizes['title'])
plt.xticks(fontsize=font_sizes['tick'])
plt.yticks(fontsize=font_sizes['tick'])
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# 7.6 Interactive PLS-DA Loading Plot
vip_threshold = 1.5
loadings_root['VIP_Score'] = vip_scores_root_filtered['VIP_Score'].values

# Custom color scale from light yellowish green to dark green, with more gradual transitions
custom_colorscale = [
    [0, "#077487"],    # Very light green (low VIP score)
    [0.2, "#15857d"],  # Light green
    [0.4, "#1c853f"],  # Bluish green
    [0.6, "#2a7a28"],  # Green-blue (mid score)
    [0.8, "#70b840"],  # Dark bluish green
    [1, "#d1f263"]     # Dark green (high VIP score)
]

# Create a more extreme size scaling
min_size = 5
max_size = 50
loadings_root['Marker_Size'] = min_size + (max_size - min_size) * (loadings_root['VIP_Score'] - loadings_root['VIP_Score'].min()) / (loadings_root['VIP_Score'].max() - loadings_root['VIP_Score'].min())

# Apply a power function to make the size difference more extreme
loadings_root['Marker_Size'] = loadings_root['Marker_Size'] ** 1.5

# Create the interactive scatter plot
fig = px.scatter(loadings_root, 
                 x='Component_1', y='Component_2', 
                 text=loadings_root.apply(lambda row: row['Metabolite'] if row['VIP_Score'] > vip_threshold else "", axis=1),
                 size='Marker_Size', color='VIP_Score',
                 size_max=max_size,  # Set the maximum size of the markers
                 title='Interactive PLS-DA Loading Plot - Root Tissue (VIP Colored)',
                 labels={
                     'Component_1': 'PLS Component 1',
                     'Component_2': 'PLS Component 2',
                     'VIP_Score': 'VIP Score'
                 },
                 hover_data={'Metabolite': True, 'VIP_Score': ':.2f'},
                 color_continuous_scale=custom_colorscale)

# Update trace appearance with grey text for metabolite names and opacity based on VIP score
fig.update_traces(
    textposition='top center', 
    textfont=dict(color='#274357'),
    marker=dict(opacity=0.6 + 0.4 * (loadings_root['VIP_Score'] - loadings_root['VIP_Score'].min()) / 
                (loadings_root['VIP_Score'].max() - loadings_root['VIP_Score'].min()))
)

# Update layout for a cleaner, more modern look
fig.update_layout(
    title_font_size=18,
    xaxis_title='PLS Component 1',
    yaxis_title='PLS Component 2',
    hoverlabel=dict(bgcolor="white", font_size=12, font_family="Rockwell"),
    font=dict(size=14),
    coloraxis_colorbar=dict(title="VIP Score"),
    plot_bgcolor='white',        # Set the plot background to white
    paper_bgcolor='white',       # Set the paper background to white
    showlegend=False             # Hide the legend as it's not needed anymore
)

# Display the plot
fig.show()

# Save the interactive plot as an HTML file
fig.write_html("interactive_loading_plot_root.html")

# 8. Save Results
# ---------------
# Save loadings
loadings_root.to_csv(r'C:\Users\ms\Desktop\data_chem_3_10\output\results\updated_PLS\pls_loadings_root.csv')

# Save VIP scores (all scores, not just filtered)
vip_scores_root_df.to_csv(r'C:\Users\ms\Desktop\data_chem_3_10\output\results\updated_PLS\vip_scores_all_root.csv', index=False)

# Save filtered VIP scores (VIP > 1)
vip_scores_root_filtered.to_csv(r'C:\Users\ms\Desktop\data_chem_3_10\output\results\updated_PLS\vip_scores_filtered_root.csv', index=False)

# Save explained variance data
pd.DataFrame({
    'Component': range(1, len(explained_variance_ratio) + 1),
    'Explained_Variance_Ratio': explained_variance_ratio,
    'Cumulative_Explained_Variance': cumulative_explained_variance
}).to_csv(r'C:\Users\ms\Desktop\data_chem_3_10\output\results\updated_PLS\explained_variance_root.csv', index=False)

print("Results saved successfully:")
print("1. PLS Loadings: pls_loadings_root.csv")
print("2. All VIP Scores: vip_scores_all_root.csv")
print("3. Filtered VIP Scores (VIP > 1): vip_scores_filtered_root.csv")
print("4. Explained Variance: explained_variance_root.csv")
print("5. Interactive Loading Plot: interactive_loading_plot_root.html")



##############################################################
############ VIP SCORE PLOT with portion #####################
############ uses all VIP score ##############################
##############################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load VIP scores from local paths
vip_scores_leaf_path = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\vip_bonferroni\vip_scores_all_leaf.csv'
vip_scores_root_path = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\vip_bonferroni\vip_scores_all_root.csv'

vip_scores_leaf = pd.read_csv(vip_scores_leaf_path)
vip_scores_root = pd.read_csv(vip_scores_root_path)

# Merge VIP scores on metabolite names to compare
vip_scores_combined = vip_scores_leaf[['Metabolite', 'VIP_Score']].merge(
    vip_scores_root[['Metabolite', 'VIP_Score']], on='Metabolite', 
    how='outer', suffixes=('_leaf', '_root')
)

# Ensure VIP_Score_leaf and VIP_Score_root are numeric (if needed, force conversion)
vip_scores_combined['VIP_Score_leaf'] = pd.to_numeric(vip_scores_combined['VIP_Score_leaf'], errors='coerce')
vip_scores_combined['VIP_Score_root'] = pd.to_numeric(vip_scores_combined['VIP_Score_root'], errors='coerce')

# Define the plot function
def create_combined_vip_plot(vip_scores_combined, thresholds, colors, font_sizes, dot_size=15):
    fig = plt.figure(figsize=(20, 5))
    gs = fig.add_gridspec(2, 2, width_ratios=[3, 1], height_ratios=[1, 1])

    # Scatter plot
    ax_scatter = fig.add_subplot(gs[:, 0])
    x = np.arange(len(vip_scores_combined))
    ax_scatter.scatter(x, vip_scores_combined['VIP_Score_leaf'], color=colors[0], label='Leaf', alpha=0.7, s=dot_size)
    ax_scatter.scatter(x, vip_scores_combined['VIP_Score_root'], color=colors[1], label='Root', alpha=0.7, s=dot_size)

    # Add horizontal lines for thresholds
    for threshold, color in zip(thresholds, ['r', 'g']):
        ax_scatter.axhline(y=threshold, color=color, linestyle='--', linewidth=2)
        ax_scatter.text(len(vip_scores_combined), threshold, f'Threshold: {threshold}', 
                        va='bottom', ha='right', color=color, fontsize=font_sizes['legend'])

    ax_scatter.set_ylabel('VIP Score', fontsize=font_sizes['label'])
    ax_scatter.set_xlabel('Molecular Features', fontsize=font_sizes['label'])
    ax_scatter.set_xticks([])

    # Calculate max y-axis value based only on numeric columns
    y_max = max(vip_scores_combined[['VIP_Score_leaf', 'VIP_Score_root']].max().max() * 1.1, max(thresholds) * 1.1)
    ax_scatter.set_ylim(0, y_max)
    ax_scatter.legend(loc='upper right', fontsize=font_sizes['legend'], markerscale=2)

    # Pie charts for thresholds
    for i, threshold in enumerate(thresholds):
        ax_pie = fig.add_subplot(gs[i, 1])
        
        leaf_count = (vip_scores_combined['VIP_Score_leaf'] > threshold).sum()
        root_count = (vip_scores_combined['VIP_Score_root'] > threshold).sum()
        
        sizes = [leaf_count, root_count]
        labels = ['Leaf', 'Root']
        
        wedges, texts, autotexts = ax_pie.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax_pie.set_title(f'VIP Scores > {threshold}', fontsize=font_sizes['title'])
        
        plt.setp(texts, size=font_sizes['pie_label'])
        plt.setp(autotexts, size=font_sizes['pie_pct'])

    # Annotations for VIP counts
    annotation_text = []
    for threshold in thresholds:
        leaf_count = (vip_scores_combined['VIP_Score_leaf'] > threshold).sum()
        root_count = (vip_scores_combined['VIP_Score_root'] > threshold).sum()
        annotation_text.append(f"VIP > {threshold}:")
        annotation_text.append(f"Leaf: {leaf_count}")
        annotation_text.append(f"Root: {root_count}")
    
    ax_scatter.text(0.02, 0.98, '\n'.join(annotation_text), transform=ax_scatter.transAxes, 
                    va='top', ha='left', fontsize=font_sizes['annotation'],
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    plt.tight_layout()
    plt.savefig('vip_scores_plot_with_pie.png', dpi=300, bbox_inches='tight')
    plt.show()

# Define color scheme and font sizes
colors = ['#72e883', '#54baaf']
font_sizes = {
    'title': 20,
    'label': 18,
    'tick': 16,
    'legend': 16,
    'annotation': 16,
    'pie_label': 15,
    'pie_pct': 15,
}

# Create combined plot with different thresholds
thresholds = [1, 2]
create_combined_vip_plot(vip_scores_combined, thresholds, colors, font_sizes, dot_size=30)







##########################################################
###### VIP SCORE PLOT with portion #######################
############ uses filtered VIP score #####################
############ VIP + FDR + bonferroni ######################
##########################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load VIP scores from local paths
vip_scores_leaf_path = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\vip_bonferroni\VIP_mann_whitney_bonferroni_fdr_leaf.csv'
vip_scores_root_path = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\vip_bonferroni\VIP_mann_whitney_bonferroni_fdr_root.csv'

vip_scores_leaf = pd.read_csv(vip_scores_leaf_path, index_col=0)
vip_scores_root = pd.read_csv(vip_scores_root_path, index_col=0)

# Merge VIP scores on metabolite names to compare
vip_scores_combined = vip_scores_leaf[['VIP_Score']].join(vip_scores_root[['VIP_Score']], lsuffix='_leaf', rsuffix='_root', how='outer')

# Define the plot function
def create_combined_vip_plot(vip_scores_combined, thresholds, colors, font_sizes, dot_size=15):
    fig = plt.figure(figsize=(20, 5))
    gs = fig.add_gridspec(2, 2, width_ratios=[3, 1], height_ratios=[1, 1])

    # Scatter plot
    ax_scatter = fig.add_subplot(gs[:, 0])
    x = np.arange(len(vip_scores_combined))
    ax_scatter.scatter(x, vip_scores_combined['VIP_Score_leaf'], color=colors[0], label='Leaf', alpha=0.7, s=dot_size)
    ax_scatter.scatter(x, vip_scores_combined['VIP_Score_root'], color=colors[1], label='Root', alpha=0.7, s=dot_size)

    for threshold, color in zip(thresholds, ['r', 'g']):
        ax_scatter.axhline(y=threshold, color=color, linestyle='--', linewidth=2)
        ax_scatter.text(len(vip_scores_combined), threshold, f'Threshold: {threshold}', 
                        va='bottom', ha='right', color=color, fontsize=font_sizes['legend'])

    ax_scatter.set_ylabel('VIP Score', fontsize=font_sizes['label'])
    ax_scatter.set_xlabel('Molecular Features', fontsize=font_sizes['label'])
    ax_scatter.set_xticks([])
    y_max = max(vip_scores_combined.max().max() * 1.1, max(thresholds) * 1.1)
    ax_scatter.set_ylim(0, y_max)
    ax_scatter.legend(loc='upper right', fontsize=font_sizes['legend'], markerscale=2)

    # Pie charts for thresholds
    for i, threshold in enumerate(thresholds):
        ax_pie = fig.add_subplot(gs[i, 1])
        
        leaf_count = (vip_scores_combined['VIP_Score_leaf'] > threshold).sum()
        root_count = (vip_scores_combined['VIP_Score_root'] > threshold).sum()
        
        sizes = [leaf_count, root_count]
        labels = ['Leaf', 'Root']
        
        wedges, texts, autotexts = ax_pie.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax_pie.set_title(f'VIP Scores > {threshold}', fontsize=font_sizes['title'])
        
        plt.setp(texts, size=font_sizes['pie_label'])
        plt.setp(autotexts, size=font_sizes['pie_pct'])

    # Annotations
    annotation_text = []
    for threshold in thresholds:
        leaf_count = (vip_scores_combined['VIP_Score_leaf'] > threshold).sum()
        root_count = (vip_scores_combined['VIP_Score_root'] > threshold).sum()
        annotation_text.append(f"VIP > {threshold}:")
        annotation_text.append(f"Leaf: {leaf_count}")
        annotation_text.append(f"Root: {root_count}")
    
    ax_scatter.text(0.02, 0.98, '\n'.join(annotation_text), transform=ax_scatter.transAxes, 
                    va='top', ha='left', fontsize=font_sizes['annotation'],
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    plt.tight_layout()
    plt.savefig('vip_scores_plot_with_pie.png', dpi=300, bbox_inches='tight')
    plt.show()

# Define color scheme and font sizes
colors = ['#72e883', '#54baaf']
font_sizes = {
    'title': 20,
    'label': 18,
    'tick': 16,
    'legend': 16,
    'annotation': 16,
    'pie_label': 15,
    'pie_pct': 15,
}

# Create combined plot with different thresholds
thresholds = [1, 2]
create_combined_vip_plot(vip_scores_combined, thresholds, colors, font_sizes, dot_size=30)
