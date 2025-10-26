"""
Missing At Random (MAR) Analysis for Metabolomics Data
---------------------------------------------------

This script tests whether missing values in metabolomics data follow a Missing At Random (MAR) 
pattern by analyzing the relationship between missingness and observed variables using 
logistic regression.

Statistical Methodology:
- For each metabolite, creates a binary indicator (0/1) for missing values.
- Uses logistic regression to model the relationship between missingness and predictors.
- Significant relationships with predictors suggest a MAR pattern.
- Non-significant relationships suggest a MCAR (Missing Completely At Random) pattern.

Predictors used:
- Genotype
- TMT (Treatment)
- Day
- Replication

Output:
- A CSV file containing the logistic regression coefficients for each metabolite.
- The coefficients indicate the strength and direction of the relationships between the 
  predictors and the likelihood of missingness.
"""
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from typing import List, Tuple, Optional

# --- Configuration ---
CONFIG = {
    "input_path": r'C:\Users\ms\Desktop\data_chem\combine\p_column_data_r_clean.csv',
    "output_path": r'C:\Users\ms\Desktop\data_chem\combine\missing_data_analysis_results.csv',
    "predictor_columns": ['Genotype', 'TMT', 'day', 'Rep'],
    "metabolite_prefix": 'P_Cluster_',
    "test_size": 0.3,
    "random_state": 42
}

def load_data(file_path: str) -> pd.DataFrame:
    """
    Load metabolomics data from a CSV file.
    
    Args:
        file_path (str): The path to the input CSV file.
        
    Returns:
        pd.DataFrame: A DataFrame containing the loaded data.
    """
    return pd.read_csv(file_path)

def create_missingness_indicator(data: pd.DataFrame, column: str) -> pd.Series:
    """
    Create a binary indicator for missing values in a specified column.
    
    Args:
        data (pd.DataFrame): The input DataFrame.
        column (str): The name of the column to analyze for missingness.
        
    Returns:
        pd.Series: A binary Series where 1 indicates a missing value and 0 indicates a present value.
    """
    return data[column].isnull().astype(int)

def prepare_predictors(data: pd.DataFrame, predictor_columns: List[str]) -> pd.DataFrame:
    """
    Prepare predictor variables for regression by applying one-hot encoding.
    
    Args:
        data (pd.DataFrame): The input DataFrame.
        predictor_columns (list): A list of column names to be used as predictors.
        
    Returns:
        pd.DataFrame: A DataFrame of one-hot encoded predictor variables.
    """
    return pd.get_dummies(data[predictor_columns], drop_first=True)

def fit_logistic_model(X_train: pd.DataFrame, X_test: pd.DataFrame, 
                       y_train: pd.Series, y_test: pd.Series) -> Tuple[LogisticRegression, pd.DataFrame, pd.DataFrame]:
    """
    Fit a logistic regression model with standardized features.
    
    Args:
        X_train (pd.DataFrame): Training feature matrix.
        X_test (pd.DataFrame): Testing feature matrix.
        y_train (pd.Series): Training target vector.
        y_test (pd.Series): Testing target vector.
        
    Returns:
        tuple: A tuple containing the fitted model, the scaled training data, and the scaled test data.
    """
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    model = LogisticRegression()
    model.fit(X_train_scaled, y_train)
    
    return model, X_train_scaled, X_test_scaled

def create_result_entry(column: str, feature_names: List[str], model: Optional[LogisticRegression] = None) -> pd.DataFrame:
    """
    Create a DataFrame entry for the regression results for a single metabolite.

    Args:
        column (str): The name of the metabolite column.
        feature_names (List[str]): The names of the predictor variables.
        model (Optional[LogisticRegression]): The fitted logistic regression model. If None, results will be null.

    Returns:
        pd.DataFrame: A single-row DataFrame containing the results.
    """
    result = {'Metabolite': column}
    if model:
        result['Intercept'] = model.intercept_[0]
        result.update(dict(zip(feature_names, model.coef_[0])))
    else:
        result['Intercept'] = None
        for name in feature_names:
            result[name] = None
            
    return pd.DataFrame([result])

def analyze_metabolite(data: pd.DataFrame, column: str, predictors: List[str]) -> Optional[pd.DataFrame]:
    """
    Analyze the missing value patterns for a single metabolite.
    
    Args:
        data (pd.DataFrame): The input DataFrame.
        column (str): The name of the metabolite column to analyze.
        predictors (list): A list of predictor variable names.
        
    Returns:
        Optional[pd.DataFrame]: A DataFrame with the analysis results for the metabolite, or None if an error occurs.
    """
    y = create_missingness_indicator(data, column)
    X = prepare_predictors(data, predictors)
    
    if len(y.unique()) <= 1:
        return create_result_entry(column, X.columns.tolist())
    
    try:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=CONFIG["test_size"], random_state=CONFIG["random_state"]
        )
        model, _, _ = fit_logistic_model(
            X_train, X_test, y_train, y_test
        )
        return create_result_entry(column, X.columns.tolist(), model)
        
    except Exception as e:
        print(f"Error analyzing {column}: {e}")
        return None

def main():
    """
    Main function to run the MAR analysis pipeline.
    """
    data = load_data(CONFIG["input_path"])
    
    results = []
    metabolite_columns = [col for col in data.columns if col.startswith(CONFIG["metabolite_prefix"])]
    
    for column in metabolite_columns:
        result = analyze_metabolite(data, column, CONFIG["predictor_columns"])
        if result is not None:
            results.append(result)
    
    if results:
        final_results = pd.concat(results, ignore_index=True)
        final_results.to_csv(CONFIG["output_path"], index=False)
        print(f"Results saved to: {CONFIG['output_path']}")
    else:
        print("No valid results were generated to save.")

if __name__ == "__main__":
    main()
