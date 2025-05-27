# 🌾 Tissue-specific Metabolomic Networks Orchestrate Osmotic Stress Adaptation in Wheat

*Uncovering the architectural principles of drought tolerance through integrated metabolomic-network analysis*

## 🎯 Project Overview
This repository contains the analytical pipeline used to investigate how wheat plants adapt to drought stress through tissue-specific metabolic networks. Our study reveals that drought-tolerant wheat varieties maintain distinct molecular organisations in leaves versus roots - leaves show highly integrated networks optimised for rapid photosynthetic responses, while roots display modular networks suited for localised environmental adaptation. These architectural differences help explain how some wheat varieties better withstand drought conditions.

Key findings:
- Identified fundamental differences in how leaves and roots organise their molecular responses to drought
- Discovered that drought-tolerant wheat has ~40% denser leaf networks compared to roots
- Found that leaf-root coordination changes over time as drought stress continues
- Validated findings using rigorous statistical approaches

### Workflow Overview
```mermaid
graph TD
    A["Clean Dataset<br/>2,471 features"] --> B("PLS-DA Feature Selection<br/>VIP scores > 1")
    B --> C("High-Confidence Features<br/>964 selected")
    
    C --> D{"Network Analysis Strategy"}
    
    D --> E
    D --> F
    D --> G

    subgraph AnalysisLayers ["Network Analysis Layers"]
        E["LAYER 1: Correlation Networks<br/>Spearman |r| > 0.7, FDR < 0.05<br/>Purpose: Establish metabolite co-regulation patterns"]
        F["LAYER 2: Topology Analysis<br/>Density, Transitivity, Modularity, Hubs<br/>Purpose: Quantify tissue-specific architecture"]
        G["LAYER 3: Bayesian Networks<br/>Hill-climbing DAG, Bootstrap n=5,000<br/>Purpose: Infer directional dependencies"]
    end
    
    E --> H("Integration and Interpretation")
    F --> H
    G --> H
    
    H --> I{"VALIDATION FRAMEWORK"}
    
    I --> I1
    I --> I2
    I --> I3

    subgraph ValidationAspects ["Statistical Robustness"]
        I1["Permutation Testing<br/>n=1,000-10,000 iterations<br/>Purpose: Assess non-random organization"]
        I2["Bootstrap Validation<br/>n=5,000 resamples<br/>Purpose: Confidence interval estimation"]
        I3["Cross-Validation<br/>Temporal stability + Module preservation<br/>Purpose: Assess dynamic network integrity"]
    end
    
    I1 --> J(("Complementary Insights<br/>Each layer provides distinct non-redundant information<br/>P < 0.001 statistical significance"))
    I2 --> J
    I3 --> J

    style A fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
    style B fill:#a5d6a7,stroke:#2e7d32,stroke-width:2px
    style C fill:#81c784,stroke:#2e7d32,stroke-width:2px
    style D fill:#2e7d32,stroke:#1b5e20,stroke-width:3px,color:#fff
    
    style E fill:#4caf50,stroke:#2e7d32,stroke-width:2px
    style F fill:#66bb6a,stroke:#2e7d32,stroke-width:2px
    style G fill:#81c784,stroke:#2e7d32,stroke-width:2px
    
    style H fill:#a5d6a7,stroke:#2e7d32,stroke-width:2px
    style I fill:#2e7d32,stroke:#1b5e20,stroke-width:3px,color:#fff
    
    style I1 fill:#4caf50,stroke:#2e7d32,stroke-width:2px
    style I2 fill:#66bb6a,stroke:#2e7d32,stroke-width:2px
    style I3 fill:#81c784,stroke:#2e7d32,stroke-width:2px
    
    style J fill:#1b5e20,stroke:#1b5e20,stroke-width:2px,color:#fff
```

### Core Capabilities
- Analysis of 2,471 molecular features across tissues
- Advanced network topology analysis
- Bayesian network validation framework
- Comprehensive temporal dynamics assessment

## 🔑 Keywords
`metabolomics` `network-analysis` `drought-tolerance` `wheat` `systems-biology` `temporal-dynamics` `tissue-specific-metabolism` `LC-MS` `bioinformatics` `plant-science` `osmotic-stress` `metabolic-networks`

## 🏗️ System Architecture

### Analysis Pipeline Structure
```
📦 Metabolomics-Analysis-Pipeline
 ├── 📂 data                          
 │
 ├── 📂 src                           
 │   ├── 📂 1_data_preprocessing    # Data pre-processing scripts
 │   │   ├── feature_filter.py       # Initial feature filtering
 │   │   ├── missing_vis.py          # Missing value visualisation
 │   │   ├── mar_test.py             # Missing at Random test
 │   │   ├── logistic_test.py        # Logistic regression test
 │   │   ├── mcar_test.py            # Missing Completely at Random test
 │   │   ├── median_impute.py        # Median imputation
 │   │   ├── rf_impute.R             # Random Forest imputation
 │   │   ├── ml_impute.py            # Machine learning imputation
 │   │   ├── impute_validate.py      # Imputation validation
 │   │   ├── impute_dist_check.py    # Distribution check after imputation
 │   │   ├── isolation_forest.py      # Isolation Forest for outliers
 │   │   ├── dim_reduce_outliers.py   # Dimensionality reduction for outliers
 │   │   ├── outlier_vis.py           # Outlier visualisation
 │   │   ├── transform_data.py        # Data transformation
 │   │   ├── normality_test.py        # Normality testing
 │   │   ├── normality_vis.py         # Normality visualisation
 │   │   ├── transform_metrics.py      # Transformation metrics
 │   │   ├── transform_eval.py         # Transformation evaluation
 │   │   ├── variance_calc.py          # Variance calculation
 │   │   └── diversity_metrics.py      # Diversity metrics calculation
 │   │
 │   ├── 📂 2_analysis                # Main analysis scripts
 │   │   ├── ini_analysis.py            # Final preprocessing
 │   │   ├── ini_analysis_summary.py    # Preprocessing summary
 │   │   ├── stat_tests.py            # Statistical tests
 │   │   ├── pls_tissue.py            # PLS analysis by tissue
 │   │   ├── spearman_network.py      # Spearman correlation network
 │   │   ├── network_decay.py         # Network decay analysis
 │   │   ├── network_summary.py       # Network summary
 │   │   ├── baysian_network.R        # baysian network
 │   │   ├── tissue_analysis.R        # Tissue-specific analysis
 │   │   └── tissue_summary.R         # Tissue analysis summary
 │   │
 │   ├── 📂 3_visualisation           # Plotting scripts
 │   │   ├── 📂 figure1              
 │   │   │   ├── network_vis.py      # Network visualisation
 │   │   │   ├── radar_plot.R        # Radar plot
 │   │   │   ├── bayesian_crosstalk.R # Bayesian network analysis
 │   │   │   ├── hub_dist.R          # Hub distribution
 │   │   │   ├── hub_decay.R         # Hub decay
 │   │   │   ├── module_org.R        # Module organisation
 │   │   │   ├── module_stability.R  # Module stability
 │   │   │   └── temporal_stability.R # Temporal stability
 │   │   │
 │   │   ├── 📂 figure2
 │   │   │   ├── tissue_temporal.R   # Tissue temporal analysis
 │   │   │   ├── temporal_corr.R     # Temporal correlation
 │   │   │   └── tissue_plot.R       # Tissue plotting
 │   │   │
 │   │   └── 📂 figure3
 │   │       └── validation_vis.R     # Validation visualisation
 │   │
 │   └── 📂 4_chemical_identification # Chemical ID scripts
 │       ├── hmdb_annotate.py         # HMDB annotation
 │       ├── gnps_annotate.py         # GNPS annotation
 │       ├── struct_classify.py       # Structural classification
 │       └── func_group.py            # Functional group analysis
 │
 ├── 📂 3D_figures                    # Interactive plot
 ├── 📂 images                         # Images
 ├── requirements.txt                  # Dependencies
 ├── environment.yaml                  # Configuration 
 └── README.md                         # Project overview


```

## 📊 Key Findings

| Network Property | Leaves | Roots |
|-----------------|--------|--------|
| Network Density | 0.354 | 0.192 |
| Transitivity | 0.740-0.804 | 0.686-0.714 |
| Modularity | 0.097-0.162 | 0.213-0.288 |
| Components | 6 | 18-21 |




## 🚀 Technical Stack

### Core Analysis Pipeline
- **Data Processing**: 
  - Pandas/NumPy (data manipulation)
  - scikit-learn (machine learning)
  - RDKit (chemical analysis)
  
- **Network Analysis**:
  - NetworkX (network construction)
  - igraph (community detection)
  - bnlearn (Bayesian networks)

- **Visualization**:
  - Matplotlib/Seaborn
  - ggplot2
  - Plotly (interactive plots)

## 🛠️ Installation & Setup


### Quick Start

1. **Clone Repository**
   ```bash
   git clone https://github.com/shoaibms/metabo-net.git
   cd metabo-net
   ```

2. **Setup Python Environment**
   ```bash
   # Create and activate environment
   conda env create -f environment.yaml
   conda activate my_environment
   
   # Install additional requirements
   pip install -r requirements.txt
   ```

3. **Setup R Environment**
   ```bash
   # Create and activate R environment
   conda env create -f environment_r.yaml
   conda activate r_env
   
   # Install R packages
   Rscript -e "source('requirements_r.txt')"
   ```

### Verify Installation
```bash
# Test Python setup
python src/1_data_preprocessing/feature_filter.py --test

# Test R setup
Rscript src/2_analysis/tissue_analysis.R --test
```

### Troubleshooting
Common issues and solutions:
- **Python package conflicts**: `conda env update -f environment.yaml`
- **R package installation fails**: `conda install -c conda-forge r-essentials`
- **Missing dependencies**: Check both `requirements.txt` and `requirements_r.txt`

## 📊 Key Network Properties

| Property | Leaves | Roots | Impact |
|----------|--------|--------|---------|
| Network Density | 0.354 | 0.192 | Higher leaf density enables rapid stress response |
| Transitivity | 0.740-0.804 | 0.686-0.714 | Better leaf network coordination |
| Modularity | 0.097-0.162 | 0.213-0.288 | Root networks more compartmentalized |
| Components | 6 | 18-21 | Roots show more independent modules |


# 📈 Analysis Pipeline
Our metabolomics data analysis pipeline consists of four major phases, with comprehensive preprocessing steps to ensure data quality and reliability.

## Detailed Data Preprocessing Workflow
```mermaid
   graph TD
    %% Major Steps with darker shades
    A([Raw Data]) --> B[Keep columns with ≥3 reps]
    B --> C[Visualise missing values]
    C --> D[Test for MCAR<br>Little's MCAR test]
    D --> E[Test for MAR<br>Logistic Regression]
    E --> F{Impute missing data}
    F --> |R|G1[Random Forest, PMM]
    F --> |Python|G2[kNN, Median, SVD, GPR, EM]
    G1 & G2 --> H[Evaluate imputation methods]
    H --> H1[EMD]
    H --> H2[Hellinger Distance]
    H --> H3[Calculate richness, Shannon entropy,<br>Simpson's diversity index, & sparsity]
    H --> H4[Visualisations: Q-Q, ECDF, KDE plots]
    H1 & H2 & H3 & H4 --> I[Select best method:<br>Random Forest]
    I --> J{Outlier detection}
    J --> K[Methods: Z-Score, IQR, Isolation Forest,<br>Elliptic Envelope, Mahalanobis, Robust PCA]
    K --> L[Evaluate outlier detection methods]
    L --> L1[PCA and t-SNE visualisations]
    L --> L2[Plots of 30 most impacted variables]
    L --> L3[Number of outliers per method]
    L1 & L2 & L3 --> M[Select method: Isolation Forest]
    M --> N[Remove outliers and<br>impute with Random Forest]
    N --> O{Data Transformation}
    O --> P[Methods: Log, Square Root, Box-Cox,<br>Yeo-Johnson, asinh, glog, Anscombe]
    P --> Q[Evaluate transformations]
    Q --> Q1[Metrics: CV, MA-transform,<br>RSD, rMAD]
    Q --> Q2[Normality tests:<br>Shapiro-Wilk, Anderson-Darling]
    Q --> Q3[Visualise: Density plots]
    Q1 & Q2 & Q3 --> R{Variable Selection}
    R --> S[Exclude variables with rMAD > 30%]
    S --> T([End: Clean Data])

    %% Styling major steps (dark green)
    style A fill:#2e7d32,stroke:#1b5e20,stroke-width:3px,color:#fff
    style F fill:#2e7d32,stroke:#1b5e20,stroke-width:3px,color:#fff
    style J fill:#2e7d32,stroke:#1b5e20,stroke-width:3px,color:#fff
    style O fill:#2e7d32,stroke:#1b5e20,stroke-width:3px,color:#fff
    style R fill:#2e7d32,stroke:#1b5e20,stroke-width:3px,color:#fff

    %% Styling endpoint (lightest green)
    style T fill:#f1f8f1,stroke:#2e7d32,stroke-width:2px

    %% Styling edges
    linkStyle default stroke:#2e7d32,stroke-width:1px
```


## Pipeline Overview

### 1️⃣ Data Preprocessing



### 2️⃣ Network Construction
Building on the clean dataset, we constructed:
- Spearman correlation networks (|r| > 0.7)
- Tissue-specific topology analysis
- Module detection using Louvain algorithm
- Hub identification and analysis

### 3️⃣ Temporal Analysis
Tracking network dynamics through:
- Cross-tissue correlation analysis
- Network stability assessment
- Module preservation analysis
- Pathway-level temporal patterns

### 4️⃣ Validation Framework
Ensuring reliability through:
- Bootstrap validation (n=5,000)
- Permutation testing
- Bayesian network validation
- Multiple null model comparisons

## 🔍 Quality Metrics
- Initial features: 4,255 (negative mode), 3,199 (positive mode)
- Final clean dataset: 2,471 molecular features
- Network validation: P < 0.001 (Bayesian analysis)

## 📚 Citation
If you use this pipeline in your research, please cite:
[Citation information will be added upon publication]

## 📝 License
This project is licensed under the MIT License - see the LICENSE file for details.
