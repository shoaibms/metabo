"""
Missing At Random (MAR) Analysis Script
-------------------------------------

This script tests whether missing data in a metabolomics dataset follows a 
Missing At Random (MAR) pattern using logistic regression. MAR is a condition 
where the probability of a value being missing depends only on other observed 
variables and not on the missing value itself.

Methodology:
1.  For each metabolite column (prefixed with 'N_Cluster_'), a binary indicator 
    variable is created (1 for missing, 0 for present).
2.  A logistic regression model is fitted for each metabolite to predict the 
    missingness indicator based on a set of experimental factors (e.g., 
    'Genotype', 'TMT', 'Day', 'Rep').
3.  The script uses one-hot encoding for categorical predictors and standardizes
    all predictor variables before fitting the model.
4.  If the model shows a significant association between the experimental factors 
    and the probability of missingness, it suggests that the data is MAR.

Input:
- A CSV file containing the metabolomics data. The file must include columns 
  for the experimental factors and metabolite measurements.

Output:
- A single CSV file containing the logistic regression results. Each row 
  corresponds to one metabolite, and the columns include the model's intercept 
  and the coefficients for each predictor variable. These coefficients indicate 
  the strength and direction of the association between each factor and the 
  likelihood of missingness for that metabolite.
"""

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

def load_data(file_path: str) -> pd.DataFrame:
    """
    Load data from a CSV file.

    Args:
        file_path (str): The path to the input CSV file.

    Returns:
        pd.DataFrame: A DataFrame containing the loaded data.
    """
    return pd.read_csv(file_path)

def create_missingness_indicator(data: pd.DataFrame, metabolite_column: str) -> pd.Series:
    """
    Create a binary indicator for missing values in a specified column.

    Args:
        data (pd.DataFrame): The input DataFrame.
        metabolite_column (str): The name of the metabolite column to check for missingness.

    Returns:
        pd.Series: A Series with binary values (1 for missing, 0 for present).
    """
    return data[metabolite_column].isnull().astype(int)

def prepare_predictors(data: pd.DataFrame, predictor_columns: list[str]) -> pd.DataFrame:
    """
    Prepare predictor variables for modeling using one-hot encoding.

    Categorical variables are converted into dummy/indicator variables.
    `drop_first=True` is used to avoid multicollinearity by removing one category per feature.

    Args:
        data (pd.DataFrame): The input DataFrame.
        predictor_columns (list[str]): A list of column names to be used as predictors.

    Returns:
        pd.DataFrame: A DataFrame of one-hot encoded predictor variables.
    """
    return pd.get_dummies(data[predictor_columns], drop_first=True, dtype=int)

def fit_logistic_model(X_train: pd.DataFrame, X_test: pd.DataFrame, y_train: pd.Series, y_test: pd.Series) -> LogisticRegression:
    """
    Fit a logistic regression model with standardized features.

    Features are scaled using StandardScaler to improve model convergence and performance.

    Args:
        X_train: Training predictor variables.
        X_test: Test predictor variables.
        y_train: Training response variable.
        y_test: Test response variable.

    Returns:
        LogisticRegression: The fitted logistic regression model.
    """
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    
    model = LogisticRegression(solver='liblinear', random_state=42) # Added for reproducibility and performance on smaller datasets
    model.fit(X_train_scaled, y_train)
    
    return model

def extract_model_results(model: LogisticRegression, metabolite_column: str, feature_names: list[str]) -> pd.DataFrame:
    """
    Extract regression coefficients and format them into a DataFrame.

    Args:
        model: The fitted logistic regression model.
        metabolite_column (str): The name of the metabolite being analyzed.
        feature_names (list[str]): The list of predictor feature names.

    Returns:
        pd.DataFrame: A single-row DataFrame with the model results.
    """
    coefficients = dict(zip(feature_names, model.coef_[0]))
    results = {
        'Metabolite': metabolite_column,
        'Intercept': model.intercept_[0],
        **coefficients
    }
    return pd.DataFrame([results])

def create_null_results(metabolite_column: str, feature_names: list[str]) -> pd.DataFrame:
    """
    Create a DataFrame with null results for cases with no variation in missingness.
    This occurs when a metabolite is either fully observed or fully missing.

    Args:
        metabolite_column (str): The name of the metabolite.
        feature_names (list[str]): The list of predictor feature names.

    Returns:
        pd.DataFrame: A DataFrame with null values for model results.
    """
    null_data = {
        'Metabolite': [metabolite_column],
        'Intercept': [None]
    }
    for feature in feature_names:
        null_data[feature] = [None]
    return pd.DataFrame(null_data)

def main():
    """
    Main function to orchestrate the MAR analysis.
    """
    # File paths
    input_path = r'C:\Users\USER\Desktop\data_chem\data\old_1\n_column_data_r_clean.csv'
    output_path = r'C:\Users\USER\Desktop\data_chem\result\n_column_data_r_clean_results.csv'
    
    data = load_data(input_path)
    
    metabolite_columns = [col for col in data.columns if col.startswith('N_Cluster_')]
    if not metabolite_columns:
        print("No metabolite columns (prefixed with 'N_Cluster_') found to analyze.")
        return

    predictor_cols = ['Genotype', 'TMT', 'Day', 'Rep']
    X = prepare_predictors(data, predictor_cols)
    
    results_list = []
    
    for column in metabolite_columns:
        y = create_missingness_indicator(data, column)
        
        # A model can only be fit if there's variation in the outcome variable.
        if len(y.unique()) <= 1:
            result = create_null_results(column, X.columns)
        else:
            try:
                X_train, X_test, y_train, y_test = train_test_split(
                    X, y, test_size=0.3, random_state=42, stratify=y
                )
                model = fit_logistic_model(X_train, X_test, y_train, y_test)
                result = extract_model_results(model, column, X.columns)
            except Exception as e:
                print(f"Error fitting model for {column}: {e}")
                # Create null results if model fitting fails for any reason
                result = create_null_results(column, X.columns)
        
        results_list.append(result)
    
    # Combine all results into a single DataFrame.
    # `concat` handles different columns across results by filling missing values with NaN.
    if results_list:
        final_results = pd.concat(results_list, ignore_index=True)
        final_results.to_csv(output_path, index=False)
        print(f"Analysis complete. Results saved to: {output_path}")

if __name__ == "__main__":
    main()
