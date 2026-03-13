
def bootstrap_resilience(*args, **kwargs):
    """Bootstrap resilience estimation (stub)."""
    return None

def classify_resilience(*args, **kwargs):
    """Classify resilience categories (stub)."""
    return None


import numpy as np

def bootstrap_correlation(data, col_x, col_y, n_bootstraps=3000, confidence_level=0.95):
    """
    Estimate the correlation between two variables with bootstrapping to provide confidence intervals.
    
    Parameters:
    - data: DataFrame containing the data.
    - col_x: Column name for the first variable.
    - col_y: Column name for the second variable.
    - n_bootstraps: Number of bootstrap samples to generate.
    - confidence_level: Confidence level for the interval (e.g., 0.95 for 95% CI).
    
    Returns:
    - A dictionary with the mean correlation and confidence interval bounds.
    """
    correlations = []
    n = len(data)
    
    # Bootstrap sampling
    for _ in range(n_bootstraps):
        sample = data.sample(n=n, replace=True)
        try:
            correlation = sample[col_x].corr(sample[col_y])
            correlations.append(correlation)
        except KeyError:
            print(f"One of the columns {col_x} or {col_y} is missing in the data.")
            continue  # Skip if one of the columns is missing

    # Calculate the confidence interval
    lower_bound = np.percentile(correlations, ((1 - confidence_level) / 2) * 100)
    upper_bound = np.percentile(correlations, (1 - (1 - confidence_level) / 2) * 100)
    
    # Return the result as a dictionary
    result = {
        'mean_correlation': np.mean(correlations),
        'ci_lower': lower_bound,
        'ci_upper': upper_bound
    }
    return result


import pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

def time_effect(data, time_column, value_column, group_column=None, alpha=0.05):
    """
    Estimate the effect of time (non-parametric) across different days.
    
    Parameters:
    - data: DataFrame containing the data.
    - time_column: Column in the DataFrame representing time points (e.g., Day 1, Day 2, Day 3).
    - value_column: Column for the variable whose time effect is being assessed.
    - group_column: Optional column representing groupings (e.g., treatment) if time effect within groups is needed.
    - alpha: Significance level for multiple testing correction.

    Returns:
    - results: Dictionary containing p-values and significant day comparisons.
    """
    results = {}
    
    # Group analysis if specified
    if group_column:
        grouped_data = data.groupby(group_column)
    else:
        grouped_data = [(None, data)]  # no grouping, single group analysis

    for group, group_df in grouped_data:
        time_points = group_df[time_column].unique()
        
        # Prepare data for the Kruskal-Wallis test
        samples = [group_df[group_df[time_column] == time][value_column].dropna() for time in time_points]
        
        # Perform Kruskal-Wallis test across time points
        kruskal_result = kruskal(*samples)
        p_value = kruskal_result.pvalue
        results[f"{group}_overall_pvalue"] = p_value

        # Pairwise testing if overall time effect is significant
        if p_value < alpha:
            pairwise_p_values = []
            time_combinations = [(t1, t2) for i, t1 in enumerate(time_points) for t2 in time_points[i + 1:]]
            for t1, t2 in time_combinations:
                sample1 = group_df[group_df[time_column] == t1][value_column].dropna()
                sample2 = group_df[group_df[time_column] == t2][value_column].dropna()
                test_stat, pair_p = kruskal(sample1, sample2)
                pairwise_p_values.append(pair_p)
            
            # Apply multiple testing correction to pairwise results
            corrected_results = multipletests(pairwise_p_values, alpha=alpha, method='bonferroni')
            significant_pairs = [(time_combinations[i], p) for i, p in enumerate(corrected_results[1]) if corrected_results[0][i]]
            results[f"{group}_significant_pairs"] = significant_pairs

    return results

# Section 0: Import necessary libraries
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os
from itertools import combinations
import multiprocessing as mp
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
from functools import partial

# Setup
warnings.filterwarnings('ignore')
tqdm.pandas()  # Enable progress bar for pandas operations

# File paths
leaf_data_path = r'C:\Users\USER\Desktop\data_chem_3_10\data\data\n_p_l.csv'
root_data_path = r'C:\Users\USER\Desktop\data_chem_3_10\data\data\n_p_r.csv'
output_dir = r'C:\Users\USER\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6e'
summary_output_dir = r'C:\Users\USER\Desktop\data_chem_3_10\output\results\initial_stat\initial_stat6e\Summary'

# Create output directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(summary_output_dir, exist_ok=True)

# Load and validate data
try:
    leaf_data = pd.read_csv(leaf_data_path)
    root_data = pd.read_csv(root_data_path)
    print(f"Data loaded successfully - Leaf: {leaf_data.shape}, Root: {root_data.shape}")
except Exception as e:
    print(f"Error loading data: {str(e)}")
    raise

# Utility Functions
def get_metabolite_cols(df):
    """Get metabolite column names."""
    return [col for col in df.columns if col.startswith('N_Cluster_') or col.startswith('P_Cluster_')]

def bootstrap_ci(data, statistic, n_iterations=5000, ci_level=0.95):
    """Calculate bootstrap confidence intervals."""
    if isinstance(data, np.ndarray):
        data = pd.Series(data)
    
    bootstrap_stats = []
    for _ in range(n_iterations):
        resampled = data.sample(n=len(data), replace=True)
        try:
            stat = statistic(resampled)
            bootstrap_stats.append(stat)
        except Exception:
            continue
            
    if not bootstrap_stats:
        return np.nan, np.nan
        
    lower_percentile = ((1 - ci_level) / 2) * 100
    upper_percentile = (1 - ((1 - ci_level) / 2)) * 100
    
    return (np.percentile(bootstrap_stats, lower_percentile),
            np.percentile(bootstrap_stats, upper_percentile))

def calculate_nonparametric_power(effect_size, n1, n2=None, alpha=0.05, n_sim=5000):
    """
    Calculate statistical power for non-parametric tests using simulation.
    
    Parameters:
    -----------
    effect_size : float
        Observed effect size (Cliff's Delta)
    n1, n2 : int
        Sample sizes
    alpha : float
        Significance level
    n_sim : int
        Number of simulations
    
    Returns:
    --------
    float : Estimated power
    """
    if n2 is None:
        n2 = n1
    
    significant_tests = 0
    
    for _ in range(n_sim):
        # Generate samples based on effect size
        x = np.random.randn(n1)
        y = np.random.randn(n2) + effect_size
        
        # Perform Mann-Whitney U test
        _, p_value = stats.mannwhitneyu(x, y, alternative='two-sided')
        
        if p_value < alpha:
            significant_tests += 1
    
    return significant_tests / n_sim

def cliff_delta(group1, group2):
    """Calculate Cliff's Delta effect size."""
    x = np.asarray(group1)
    y = np.asarray(group2)
    
    nx = len(x)
    ny = len(y)
    
    if nx < 2 or ny < 2:
        return np.nan
        
    dominance = 0
    for i in x:
        for j in y:
            if i > j:
                dominance += 1
            elif i < j:
                dominance -= 1
                
    return dominance / (nx * ny)

def wilcoxon_r(statistic, n):
    """Calculate effect size r for Wilcoxon test."""
    return abs(statistic) / np.sqrt(n)

def adjust_pvalues(df, method='fdr_bh'):
    """Adjust p-values for multiple comparisons."""
    if 'P_value' not in df.columns:
        return df
        
    df['Adjusted_P_value'] = multipletests(df['P_value'].fillna(1), 
                                         method=method)[1]
    return df

# Save and Summary Functions
def save_results(df, analysis_name):
    """Save analysis results to CSV."""
    output_path = os.path.join(output_dir, f'{analysis_name}.csv')
    df.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")

def save_summary(df, analysis_name):
    """Save summary results to CSV."""
    output_path = os.path.join(summary_output_dir, f'{analysis_name}_summary.csv')
    df.to_csv(output_path, index=False)
    print(f"Summary saved to {output_path}")


# Section 1: Tissue Comparison
def tissue_comparison():
    """Perform tissue comparison analysis."""
    results = []
    metabolite_pairs = list(zip(get_metabolite_cols(leaf_data), 
                               get_metabolite_cols(root_data)))
    
    total_comparisons = len(leaf_data['Genotype'].unique()) * \
                       len(leaf_data['Day'].unique()) * \
                       len(metabolite_pairs)
                       
    with tqdm(total=total_comparisons, desc="Tissue Comparison") as pbar:
        for genotype in leaf_data['Genotype'].unique():
            for day in leaf_data['Day'].unique():
                for l_met, r_met in metabolite_pairs:
                    try:
                        # Get data for comparison
                        leaf_values = leaf_data[(leaf_data['Genotype'] == genotype) & 
                                              (leaf_data['Day'] == day)][l_met]
                        root_values = root_data[(root_data['Genotype'] == genotype) & 
                                              (root_data['Day'] == day)][r_met]
                        
                        if len(leaf_values) >= 3 and len(root_values) >= 3:
                            # Calculate statistics
                            statistic, p_value = stats.mannwhitneyu(
                                leaf_values, root_values, 
                                alternative='two-sided'
                            )
                            effect_size = cliff_delta(leaf_values, root_values)
                            
                            # Calculate confidence intervals
                            ci_lower, ci_upper = bootstrap_ci(
                                pd.concat([leaf_values, root_values]),
                                lambda x: np.median(x[:len(leaf_values)]) - \
                                         np.median(x[len(leaf_values):])
                            )
                            
                            # Calculate power using non-parametric method
                            power = calculate_nonparametric_power(
                                effect_size,
                                len(leaf_values),
                                len(root_values)
                            )
                            
                            results.append({
                                'Genotype': genotype,
                                'Day': day,
                                'Leaf_Metabolite': l_met,
                                'Root_Metabolite': r_met,
                                'Statistic': statistic,
                                'P_value': p_value,
                                'Effect_size': effect_size,
                                'CI_lower': ci_lower,
                                'CI_upper': ci_upper,
                                'Power': power,
                                'N_leaf': len(leaf_values),
                                'N_root': len(root_values)
                            })
                    except Exception as e:
                        print(f"Error in comparison {l_met}-{r_met}: {str(e)}")
                    finally:
                        pbar.update(1)
    
    # Create DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    results_df = adjust_pvalues(results_df)
    
    return results_df



# Section 2: Genotype Comparison
def genotype_comparison():
    """
    Perform comprehensive genotype comparison analysis.
    
    Returns:
    --------
    pandas.DataFrame : Results of genotype comparison analysis
    """
    results = []
    total_comparisons = 0
    
    # Calculate total comparisons for progress bar
    for tissue, data in [('L', leaf_data), ('R', root_data)]:
        total_comparisons += len(data['Day'].unique()) * \
                           len(data['Treatment'].unique()) * \
                           len(get_metabolite_cols(data))
    
    with tqdm(total=total_comparisons, desc="Genotype Comparison") as pbar:
        for tissue, data in [('L', leaf_data), ('R', root_data)]:
            metabolites = get_metabolite_cols(data)
            
            for day in data['Day'].unique():
                for treatment in data['Treatment'].unique():
                    for met in metabolites:
                        try:
                            # Get data for comparison
                            g1_values = data[(data['Genotype'] == 'G1') & 
                                           (data['Day'] == day) & 
                                           (data['Treatment'] == treatment)][met]
                            g2_values = data[(data['Genotype'] == 'G2') & 
                                           (data['Day'] == day) & 
                                           (data['Treatment'] == treatment)][met]
                            
                            if len(g1_values) >= 3 and len(g2_values) >= 3:
                                # Calculate statistics
                                statistic, p_value = stats.mannwhitneyu(
                                    g1_values, g2_values, 
                                    alternative='two-sided'
                                )
                                effect_size = cliff_delta(g1_values, g2_values)
                                
                                # Calculate confidence intervals
                                ci_lower, ci_upper = bootstrap_ci(
                                    pd.concat([g1_values, g2_values]),
                                    lambda x: np.median(x[:len(g1_values)]) - \
                                             np.median(x[len(g1_values):])
                                )
                                
                                # Calculate power using non-parametric method
                                power = calculate_nonparametric_power(
                                    effect_size,
                                    len(g1_values),
                                    len(g2_values)
                                )
                                
                                results.append({
                                    'Tissue': tissue,
                                    'Day': day,
                                    'Treatment': treatment,
                                    'Metabolite': met,
                                    'Statistic': statistic,
                                    'P_value': p_value,
                                    'Effect_size': effect_size,
                                    'CI_lower': ci_lower,
                                    'CI_upper': ci_upper,
                                    'Power': power,
                                    'N_G1': len(g1_values),
                                    'N_G2': len(g2_values),
                                    'G1_median': np.median(g1_values),
                                    'G2_median': np.median(g2_values)
                                })
                                
                        except Exception as e:
                            print(f"Error in comparison for {met}: {str(e)}")
                        finally:
                            pbar.update(1)
    
    # Create DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    results_df = adjust_pvalues(results_df)
    
    # Add fold change
    results_df['Log2_Fold_Change'] = np.log2(results_df['G2_median'] / results_df['G1_median'])
    
    return results_df



# Section 3: Treatment Effect Analysis
def treatment_effect():
    """
    Analyze treatment effects using paired comparisons.
    
    Returns:
    --------
    pandas.DataFrame : Results of treatment effect analysis
    """
    results = []
    total_comparisons = 0
    
    # Calculate total comparisons
    for tissue, data in [('L', leaf_data), ('R', root_data)]:
        total_comparisons += len(data['Genotype'].unique()) * \
                           len(data['Day'].unique()) * \
                           len(get_metabolite_cols(data))
    
    with tqdm(total=total_comparisons, desc="Treatment Effect Analysis") as pbar:
        for tissue, data in [('L', leaf_data), ('R', root_data)]:
            metabolites = get_metabolite_cols(data)
            
            for genotype in data['Genotype'].unique():
                for day in data['Day'].unique():
                    for met in metabolites:
                        try:
                            # Get control and treatment data
                            control = data[(data['Genotype'] == genotype) & 
                                         (data['Day'] == day) & 
                                         (data['Treatment'] == 0)][met]
                            treated = data[(data['Genotype'] == genotype) & 
                                         (data['Day'] == day) & 
                                         (data['Treatment'] == 1)][met]
                            
                            if len(control) >= 3 and len(treated) >= 3:
                                # Calculate statistics
                                try:
                                    statistic, p_value = stats.wilcoxon(control, treated)
                                    effect_size = wilcoxon_r(statistic, len(control))
                                except:
                                    # Fallback to Mann-Whitney U test if Wilcoxon fails
                                    statistic, p_value = stats.mannwhitneyu(
                                        control, treated,
                                        alternative='two-sided'
                                    )
                                    effect_size = cliff_delta(control, treated)
                                
                                # Calculate confidence intervals
                                ci_lower, ci_upper = bootstrap_ci(
                                    pd.concat([control, treated]),
                                    lambda x: np.median(x[len(control):]) - \
                                             np.median(x[:len(control)])
                                )
                                
                                results.append({
                                    'Tissue': tissue,
                                    'Genotype': genotype,
                                    'Day': day,
                                    'Metabolite': met,
                                    'Statistic': statistic,
                                    'P_value': p_value,
                                    'Effect_size': effect_size,
                                    'CI_lower': ci_lower,
                                    'CI_upper': ci_upper,
                                    'N_control': len(control),
                                    'N_treated': len(treated),
                                    'Control_median': np.median(control),
                                    'Treated_median': np.median(treated)
                                })
                        
                        except Exception as e:
                            print(f"Error in treatment effect analysis for {met}: {str(e)}")
                        finally:
                            pbar.update(1)
    
    # Create DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    results_df = adjust_pvalues(results_df)
    
    # Add fold change
    results_df['Log2_Fold_Change'] = np.log2(
        results_df['Treated_median'] / results_df['Control_median']
    )
    
    return results_df



# Section 4: Stress Level Comparison
def stress_level_comparison():
    """
    Compare stress levels between conditions.
    
    Returns:
    --------
    pandas.DataFrame : Results of stress level comparison
    """
    results = []
    total_comparisons = 0
    
    # Calculate total comparisons
    for tissue, data in [('L', leaf_data), ('R', root_data)]:
        total_comparisons += len(data['Genotype'].unique()) * \
                           len(data['Day'].unique()) * \
                           len(get_metabolite_cols(data))
    
    with tqdm(total=total_comparisons, desc="Stress Level Analysis") as pbar:
        for tissue, data in [('L', leaf_data), ('R', root_data)]:
            metabolites = get_metabolite_cols(data)
            
            for genotype in data['Genotype'].unique():
                for day in data['Day'].unique():
                    for met in metabolites:
                        try:
                            # Get harsh and mild stress data
                            harsh = data[(data['Genotype'] == genotype) & 
                                       (data['Day'] == day) & 
                                       (data['Batch'] == 'B1')][met]
                            mild = data[(data['Genotype'] == genotype) & 
                                      (data['Day'] == day) & 
                                      (data['Batch'] == 'B2')][met]
                            
                            if len(harsh) >= 3 and len(mild) >= 3:
                                # Calculate statistics
                                statistic, p_value = stats.mannwhitneyu(
                                    harsh, mild,
                                    alternative='two-sided'
                                )
                                effect_size = cliff_delta(harsh, mild)
                                
                                # Calculate confidence intervals
                                ci_lower, ci_upper = bootstrap_ci(
                                    pd.concat([harsh, mild]),
                                    lambda x: np.median(x[:len(harsh)]) - \
                                             np.median(x[len(harsh):])
                                )
                                
                                # Calculate power using non-parametric method
                                power = calculate_nonparametric_power(
                                    effect_size,
                                    len(harsh),
                                    len(mild)
                                )
                                
                                results.append({
                                    'Tissue': tissue,
                                    'Genotype': genotype,
                                    'Day': day,
                                    'Metabolite': met,
                                    'Statistic': statistic,
                                    'P_value': p_value,
                                    'Effect_size': effect_size,
                                    'CI_lower': ci_lower,
                                    'CI_upper': ci_upper,
                                    'Power': power,
                                    'N_harsh': len(harsh),
                                    'N_mild': len(mild),
                                    'Harsh_median': np.median(harsh),
                                    'Mild_median': np.median(mild)
                                })
                        
                        except Exception as e:
                            print(f"Error in stress level comparison for {met}: {str(e)}")
                        finally:
                            pbar.update(1)
    
    # Create DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    results_df = adjust_pvalues(results_df)
    
    # Add stress response ratio
    results_df['Stress_Response_Ratio'] = np.log2(
        results_df['Harsh_median'] / results_df['Mild_median']
    )
    
    return results_df



# Section 5: Temporal Analysis
def classify_temporal_pattern(row):
    """
    Classify temporal patterns based on median values and statistical significance.
    
    Parameters:
    -----------
    row : pandas.Series
        Row containing temporal metabolite data
    
    Returns:
    --------
    str : Pattern classification
    """
    if row['P_value'] >= 0.05:
        return 'No Significant Change'
    
    d1, d2, d3 = row['Day1_median'], row['Day2_median'], row['Day3_median']
    
    # Calculate fold changes
    changes = [
        ((d2 - d1) / abs(d1)) * 100 if d1 != 0 else float('inf'),
        ((d3 - d2) / abs(d2)) * 100 if d2 != 0 else float('inf')
    ]
    
    threshold = 20  # 20% change threshold for biological significance
    
    # Pattern classification
    if all(abs(c) < threshold for c in changes):
        return 'Stable'
    
    if d1 < d2 < d3:
        return 'Continuous Increase'
    if d1 > d2 > d3:
        return 'Continuous Decrease'
    
    if d2 > d1 and d2 > d3:
        return 'Peak Response'
    if d2 < d1 and d2 < d3:
        return 'Valley Response'
    
    if abs(d2 - d1) > abs(d3 - d2):
        return 'Early Response'
    
    return 'Late Response'

def temporal_analysis():
    """
    Analyze temporal patterns in metabolite responses including:
    - Overall temporal changes (Friedman test)
    - Effect sizes (Kendall's W)
    - Temporal trends
    - Pattern classification
    - Statistical power analysis
    
    Returns:
    --------
    pandas.DataFrame : Comprehensive temporal analysis results
    """
    results = []
    
    # Calculate total combinations for progress tracking
    total_combinations = len(leaf_data['Genotype'].unique()) * \
                        len(leaf_data['Treatment'].unique()) * \
                        (len(get_metabolite_cols(leaf_data)) + len(get_metabolite_cols(root_data)))
    
    with tqdm(total=total_combinations, desc="Temporal Analysis") as pbar:
        for tissue, data in [('L', leaf_data), ('R', root_data)]:
            metabolites = get_metabolite_cols(data)
            
            # Process in batches for memory efficiency
            batch_size = 50
            for batch_start in range(0, len(metabolites), batch_size):
                batch_mets = metabolites[batch_start:batch_start + batch_size]
                
                for genotype in data['Genotype'].unique():
                    for treatment in data['Treatment'].unique():
                        for met in batch_mets:
                            try:
                                # Extract time series data
                                day_groups = [
                                    data[(data['Genotype'] == genotype) & 
                                         (data['Treatment'] == treatment) & 
                                         (data['Day'] == day)][met] 
                                    for day in [1, 2, 3]
                                ]
                                
                                if all(len(group) >= 3 for group in day_groups):
                                    # Friedman test for temporal changes
                                    statistic, p_value = stats.friedmanchisquare(*day_groups)
                                    
                                    # Calculate Kendall's W for effect size
                                    k = len(day_groups)
                                    n = len(day_groups[0])
                                    kendall_w = statistic / (n * k * (k - 1))
                                    
                                    # Calculate temporal dynamics
                                    def time_effect(data):
                                        groups = np.array_split(data, 3)
                                        return np.max([np.median(g) for g in groups]) - \
                                               np.min([np.median(g) for g in groups])
                                    
                                    # Bootstrap confidence intervals
                                    ci_lower, ci_upper = bootstrap_ci(
                                        pd.concat(day_groups),
                                        time_effect
                                    )
                                    
                                    # Calculate temporal trend using robust regression
                                    median_values = [np.median(group) for group in day_groups]
                                    trend = np.polyfit(range(1, 4), median_values, 1)[0]
                                    
                                    # Calculate non-parametric power
                                    power = calculate_nonparametric_power(
                                        kendall_w,
                                        len(day_groups[0]),
                                        len(day_groups[0]),
                                        n_sim=5000
                                    )
                                    
                                    # Store results
                                    results.append({
                                        'Tissue': tissue,
                                        'Genotype': genotype,
                                        'Treatment': treatment,
                                        'Metabolite': met,
                                        'Statistic': statistic,
                                        'P_value': p_value,
                                        'Effect_size': kendall_w,
                                        'CI_lower': ci_lower,
                                        'CI_upper': ci_upper,
                                        'Temporal_trend': trend,
                                        'Power': power,
                                        'Day1_median': median_values[0],
                                        'Day2_median': median_values[1],
                                        'Day3_median': median_values[2],
                                        'N_per_timepoint': n
                                    })
                            
                            except Exception as e:
                                print(f"Error in temporal analysis for {met}: {str(e)}")
                            finally:
                                pbar.update(1)
    
    # Create DataFrame
    results_df = pd.DataFrame(results)
    
    # Adjust p-values for multiple testing
    results_df = adjust_pvalues(results_df)
    
    # Add temporal pattern classification
    results_df['Pattern'] = results_df.apply(classify_temporal_pattern, axis=1)
    
    # Add additional temporal metrics
    results_df['Response_Magnitude'] = results_df.apply(
        lambda x: max(abs(x['Day2_median'] - x['Day1_median']),
                     abs(x['Day3_median'] - x['Day2_median'])), axis=1)
    
    results_df['Response_Speed'] = results_df.apply(
        lambda x: abs(x['Day2_median'] - x['Day1_median']) /
                 abs(x['Day3_median'] - x['Day1_median'])
        if abs(x['Day3_median'] - x['Day1_median']) > 0 else 0, axis=1)
    
    return results_df



# Section 6: Integrated Response Analysis
def integrated_response():
    """
    Analyze integrated metabolic responses.
    
    Returns:
    --------
    pandas.DataFrame : Results of integrated response analysis
    """
    results = []
    for tissue, data in [('L', leaf_data), ('R', root_data)]:
        metabolites = get_metabolite_cols(data)
        for day in data['Day'].unique():
            for treatment in data['Treatment'].unique():
                subset = data[(data['Day'] == day) & (data['Treatment'] == treatment)]
                
                if len(subset) > 0:
                    # Calculate integrated response using ranked data
                    ranked_data = subset[metabolites].rank()
                    integrated_ranks = ranked_data.mean(axis=1)
                    
                    # Split by genotype
                    g1_ranks = integrated_ranks[subset['Genotype'] == 'G1']
                    g2_ranks = integrated_ranks[subset['Genotype'] == 'G2']
                    
                    if len(g1_ranks) >= 3 and len(g2_ranks) >= 3:
                        # Calculate statistics
                        statistic, p_value = stats.mannwhitneyu(g1_ranks, g2_ranks, alternative='two-sided')
                        effect_size = cliff_delta(g1_ranks, g2_ranks)
                        
                        # Bootstrap confidence intervals
                        ci_lower, ci_upper = bootstrap_ci(
                            pd.concat([g1_ranks, g2_ranks]),
                            lambda x: np.median(x[:len(g1_ranks)]) - np.median(x[len(g1_ranks):])
                        )
                        
                        # Calculate non-parametric power
                        power = calculate_nonparametric_power(
                            effect_size,
                            len(g1_ranks),
                            len(g2_ranks),
                            n_sim=5000
                        )
                        
                        results.append({
                            'Tissue': tissue,
                            'Day': day,
                            'Treatment': treatment,
                            'Statistic': statistic,
                            'P_value': p_value,
                            'Effect_size': effect_size,
                            'CI_lower': ci_lower,
                            'CI_upper': ci_upper,
                            'Power': power,
                            'N_G1': len(g1_ranks),
                            'N_G2': len(g2_ranks),
                            'G1_median': np.median(g1_ranks),
                            'G2_median': np.median(g2_ranks)
                        })
    
    # Create DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    results_df = adjust_pvalues(results_df)
    
    return results_df



# Section 8: Resilience Index Analysis
def resilience_index():
    """
    Calculate resilience indices for metabolites.
    
    Returns:
    --------
    pandas.DataFrame : Results of resilience index analysis
    """
    results = []
    
    total_comparisons = 0
    for tissue, data in [('L', leaf_data), ('R', root_data)]:
        total_comparisons += len(data['Genotype'].unique()) * \
                           len(data['Day'].unique()) * \
                           len(get_metabolite_cols(data))
    
    with tqdm(total=total_comparisons, desc="Resilience Analysis") as pbar:
        for tissue, data in [('L', leaf_data), ('R', root_data)]:
            metabolites = get_metabolite_cols(data)
            for genotype in data['Genotype'].unique():
                for day in data['Day'].unique():
                    for met in metabolites:
                        try:
                            control = data[(data['Genotype'] == genotype) & 
                                         (data['Day'] == day) & 
                                         (data['Treatment'] == 0)][met]
                            treated = data[(data['Genotype'] == genotype) & 
                                         (data['Day'] == day) & 
                                         (data['Treatment'] == 1)][met]
                            
                            if len(control) >= 3 and len(treated) >= 3:
                                # Calculate resilience using robust statistics
                                control_median = np.median(control)
                                treated_median = np.median(treated)
                                
                                if control_median != 0:  # Avoid division by zero
                                    resilience = treated_median / control_median
                                    
                                    # Calculate bootstrap confidence intervals
                                    def bootstrap_resilience(data):
                                        ctrl = data[:len(control)]
                                        trt = data[len(control):]
                                        return np.median(trt) / np.median(ctrl)
                                    
                                    ci_lower, ci_upper = bootstrap_ci(
                                        pd.concat([control, treated]),
                                        bootstrap_resilience
                                    )
                                    
                                    # Calculate robust effect size
                                    effect_size = cliff_delta(treated/np.median(control), 
                                                            control/np.median(control))
                                    
                                    # Calculate recovery score
                                    recovery_score = 1 - abs(1 - resilience)
                                    
                                    # Calculate non-parametric power for resilience
                                    power = calculate_nonparametric_power(
                                        effect_size,
                                        len(control),
                                        len(treated),
                                        n_sim=5000
                                    )
                                    
                                    results.append({
                                        'Tissue': tissue,
                                        'Genotype': genotype,
                                        'Day': day,
                                        'Metabolite': met,
                                        'Resilience_Index': resilience,
                                        'Effect_size': effect_size,
                                        'Power': power,
                                        'CI_lower': ci_lower,
                                        'CI_upper': ci_upper,
                                        'Recovery_Score': recovery_score,
                                        'N_control': len(control),
                                        'N_treated': len(treated),
                                        'Control_median': control_median,
                                        'Treated_median': treated_median
                                    })
                        
                        except Exception as e:
                            print(f"Error in resilience analysis for {met}: {str(e)}")
                        finally:
                            pbar.update(1)
    
    results_df = pd.DataFrame(results)
    
    # Add resilience classification
    results_df['Resilience_Class'] = results_df.apply(classify_resilience, axis=1)
    
    return results_df



def classify_resilience(row):
    """Classify resilience based on index and recovery score."""
    if row['Recovery_Score'] >= 0.9:
        return 'Highly Resilient'
    elif row['Recovery_Score'] >= 0.7:
        return 'Moderately Resilient'
    elif row['Recovery_Score'] >= 0.5:
        return 'Slightly Resilient'
    else:
        return 'Not Resilient'


# Right before this section in initial_analysis_V6c.py

def save_checkpoint(partial_results, analysis_name):
    """Save intermediate results during long computations"""
    checkpoint_dir = os.path.join(output_dir, 'checkpoints')
    os.makedirs(checkpoint_dir, exist_ok=True)
    checkpoint_path = os.path.join(checkpoint_dir, f'{analysis_name}_checkpoint.pkl')
    
    with open(checkpoint_path, 'wb') as f:
        pickle.dump(partial_results, f)



# Main execution function
def run_all_analyses():
    """Execute analyses with checkpoint capability and progress tracking"""
    analyses = {
        'tissue_comparison': tissue_comparison,
        'genotype_comparison': genotype_comparison,
        'treatment_effect': treatment_effect,
        'stress_level_comparison': stress_level_comparison,
        'temporal_analysis': temporal_analysis,
        'integrated_response': integrated_response,
        'resilience_index': resilience_index
    }

    # Check for completed analyses
    completed_analyses = []
    for name in analyses.keys():
        result_path = os.path.join(output_dir, f'{name}.csv')
        if os.path.exists(result_path):
            completed_analyses.append(name)
            print(f"Found completed analysis: {name}")
    
    # Create checkpoint directory
    checkpoint_dir = os.path.join(output_dir, 'checkpoints')
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    # Run remaining analyses
    for name, func in analyses.items():
        if name in completed_analyses:
            print(f"Skipping completed analysis: {name}")
            continue
            
        print(f"Running {name}...")
        try:
            # Check for checkpoint
            checkpoint_path = os.path.join(checkpoint_dir, f'{name}_checkpoint.pkl')
            if os.path.exists(checkpoint_path):
                print(f"Resuming {name} from checkpoint...")
                with open(checkpoint_path, 'rb') as f:
                    results = pickle.load(f)
            else:
                results = func()
            
            # Save results and cleanup checkpoint
            save_results(results, name)
            if os.path.exists(checkpoint_path):
                os.remove(checkpoint_path)
            print(f"{name} completed successfully.")
            
        except Exception as e:
            print(f"Error in {name}: {str(e)}")
            # Save checkpoint on error
            if 'results' in locals() and results is not None:
                with open(checkpoint_path, 'wb') as f:
                    pickle.dump(results, f)
                print(f"Saved checkpoint for {name}")
            return False
    
    print("\nAll analyses completed successfully.")
    return True

def save_analysis_progress(completed_name):
    """Track completed analyses"""
    progress_file = os.path.join(output_dir, 'analysis_progress.txt')
    with open(progress_file, 'a') as f:
        f.write(f"{completed_name}\t{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Execution guard
if __name__ == "__main__":
    # Enable multiprocessing support for Windows
    mp.freeze_support()
    
    # Required imports
    import pickle
    from datetime import datetime
    
    # Start analysis
    print("Starting metabolomics data analysis...")
    print(f"Input data shapes - Leaf: {leaf_data.shape}, Root: {root_data.shape}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        success = run_all_analyses()
        if success:
            print("Analysis pipeline completed successfully.")
        else:
            print("Analysis pipeline completed with errors. Check logs and restart to resume.")
    except Exception as e:
        print(f"Critical error in analysis: {str(e)}")
    finally:
        print("\nAnalysis pipeline execution finished.")