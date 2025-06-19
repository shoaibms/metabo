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
    A["Clean Dataset<br/>2471 features"] --> B("PLS-DA Feature Selection<br/>VIP scores > 1")
    B --> C("High-Confidence Features<br/>964 selected")
    
    C --> D{"Network Analysis Strategy"}
    
    D --> E
    D --> F
    D --> G

    subgraph AnalysisLayers ["Network Analysis Layers"]
        E["LAYER 1: Correlation Networks<br/>Spearman |r| > 0.7, FDR < 0.05<br/>Purpose: Establish initial metabolite associations"]
        F["LAYER 2: Topology Analysis<br/>Density, Transitivity, Modularity, Hubs<br/>Purpose: Quantify network architecture"]
        G["LAYER 3: Bayesian Networks<br/>Structure Learning Hill-Climbing DAG<br/>Purpose: Infer potential dependencies and validate non-randomness"]
    end
    
    E --> H("Integration and Interpretation")
    F --> H
    G --> H
    
    H --> I{"VALIDATION FRAMEWORK"}
    
    I --> I1
    I --> I2
    I --> I3

    subgraph ValidationAspects ["Validation Aspects"]
        I1["Statistical Robustness<br/>Permutation Testing n=5000<br/>Purpose: Assess significance vs random chance"]
        I2["Temporal Stability<br/>Cross-time point consistency<br/>Purpose: Evaluate dynamic network integrity"]
        I3["Biological Significance<br/>Pathway Enrichment and Known Interactions<br/>Purpose: Link network features to function"]
    end
    
    I1 --> J(("Complementary Insights<br/>Each layer provides distinct non-redundant information"))
    I2 --> J
    I3 --> J
    
    J --> K["INTEGRATION RATIONALE<br/>Spearman alone = association without direction<br/>Topology alone = structure without relationships<br/>Bayesian alone = causality without architecture<br/>Combined = Complete network characterisation"]

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
    style K fill:#0d4f3c,stroke:#1b5e20,stroke-width:3px,color:#fff
```

### Core Capabilities
- Analysis of 2,471 molecular features across tissues
- Comprehensive validation framework (permutation + bootstrap + cross-validation)
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

- **Visualisation**:
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
| Modularity | 0.097-0.162 | 0.213-0.288 | Root networks more compartmentalised |
| Components | 6 | 18-21 | Roots show more independent modules |


# 📈 Analysis Pipeline
Our metabolomics data analysis pipeline consists of four major phases, with comprehensive preprocessing steps to ensure data quality and reliability.

## Detailed Data Preprocessing Workflow
```mermaid

flowchart TD
    A["Start: Raw Data"] --> B["Keep columns with at least 3 replicates"]
    B --> C["Visualise missing values"]
    C --> D["Test for MCAR<br>Little's MCAR test"]
    D -->|"Not MCAR"| E["Test for MAR<br>Logistic Regression"]
    E -->|"MAR or MCAR"| F["Impute missing data"]
    
    F --> G1["R: Random Forest, PMM"]
    F --> G2["Python: kNN, Median, SVD, GPR, EM"]
    
    G1 --> H["Evaluate imputation methods"]
    G2 --> H
    
    H --> I1["EMD"]
    H --> I2["Hellinger Distance"]
    H --> I3["Calculated richness, Shannon entropy,<br>Simpson's diversity index, & sparsity"]
    H --> I4["Visualisations: Q-Q, ECDF, KDE plots"]
    
    I1 --> J["Select best method:<br>Random Forest"]
    I2 --> J
    I3 --> J
    I4 --> J
    
    J --> K["Outlier detection"]
    
    K --> L["Methods: Z-Score, IQR, Isolation Forest,<br>Elliptic Envelope, Mahalanobis, Robust PCA"]
    
    L --> M["Evaluate outlier detection methods"]
    
    M --> N1["PCA and t-SNE visualisations"]
    M --> N2["Plots of 30 most impacted variables"]
    M --> N3["Number of outliers per method"]
    
    N1 --> O["Select method: Isolation Forest"]
    N2 --> O
    N3 --> O
    
    O --> P["Remove outliers and<br>impute with Random Forest"]
    
    P --> Q["Data Transformation"]
    
    Q --> R["Methods: Log, Square Root, Box-Cox,<br>Yeo-Johnson, asinh, glog, Anscombe"]
    
    R --> S["Evaluate transformations"]
    
    S --> T1["Metrics: CV, MA-transform,<br>rSD, rMAD"]
    S --> T2["Normality tests:<br>Shapiro-Wilk, Anderson-Darling"]
    S --> T3["Visualise: Density plots"]
    
    T1 --> U["Variable Selection"]
    T2 --> U
    T3 --> U
    
    U --> V["Exclude variables with rMAD > 30%"]
    
    V --> W["End: Clean Data"]

    %% Style definitions - different shades of green based on hierarchy
    classDef mainDecision fill:#5d9b7e,stroke:#000,stroke-width:1.5px,shape:diamond,color:black,font-size:16px
    classDef process fill:#a8e6cf,stroke:#000,stroke-width:1px,color:black,font-size:14px
    classDef evaluateProcess fill:#97d1bc,stroke:#000,stroke-width:1px,color:black,font-size:14px
    classDef methodProcess fill:#c3eadc,stroke:#000,stroke-width:1px,color:black,font-size:14px
    
    %% Apply styles
    class A,F,K,Q,U mainDecision
    class B,C,D,E,J,O,P,V,W process
    class H,M,S evaluateProcess
    class G1,G2,I1,I2,I3,I4,L,N1,N2,N3,R,T1,T2,T3 methodProcess
```


## Pipeline Overview

### 1️⃣ Data Preprocessing



### 2️⃣ Multi-Layer Network Analysis
Our three-layer analytical framework addresses distinct biological questions:
- **Layer 1**: Spearman correlation networks (|r| > 0.7, FDR < 0.05) - metabolite co-regulation
- **Layer 2**: Topology analysis (density, transitivity, modularity) - architectural principles  
- **Layer 3**: Bayesian networks (Hill-climbing DAG) - directional dependencies

### 3️⃣ Temporal Analysis
Tracking network dynamics through:
- Cross-tissue correlation analysis
- Network stability assessment
- Module preservation analysis
- Pathway-level temporal patterns

### 4️⃣ Statistical Validation Framework
Multi-tier validation ensuring robustness:
- **Permutation Testing**: n=1,000-10,000 iterations (optimised by analysis type)
- **Bootstrap Validation**: n=5,000 resamples for confidence intervals
- **Cross-Validation**: Temporal stability and module preservation
- **Significance Threshold**: P < 0.001 across all analyses
  
## 🔍 Quality Metrics
- Initial features: 4,255 (negative mode), 3,199 (positive mode)
- Final clean dataset: 2,471 molecular features
- Network validation: P < 0.001 (Bayesian analysis)

## 📚 Citation
If you use this pipeline in your research, please cite:
[Citation information will be added upon publication]

## 📝 License
This project is licensed under the MIT License - see the LICENSE file for details.
