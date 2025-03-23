# Clinical Trial Analysis Dashboard

A Bayesian approach to analyzing clinical trial outcomes using Quarto, R, and Python.

![Dashboard Preview](images/study_design.svg)

## Overview

This dashboard presents analyses of a hypothetical clinical trial examining the efficacy of Treatment X versus standard care in patients with Condition Y. The analyses incorporate Bayesian methods to quantify uncertainty and provide probabilistic interpretations of treatment effects.

## Features

- **Bayesian analysis** of treatment effects with full posterior distributions
- **Interactive visualizations** of patient outcomes and biomarker data
- **Dynamic reporting** with parameter adjustments
- **Probabilistic decision metrics** for clinical relevance
- **Subgroup analysis** for treatment effect heterogeneity

## Project Structure

```
clinical-trial-analysis/
├── _quarto.yml                  # Quarto configuration file
├── index.qmd                    # Landing page
├── static-report.qmd            # Static analysis report
├── interactive-report-static.qmd # Interactive analysis dashboard
├── data/
│   ├── trial_data.csv           # Clinical trial data
│   └── biomarker_data.csv       # Biomarker data
├── scripts/
│   ├── data_preparation.py      # Python script for data prep
│   ├── bayesian_models.R        # R script with Bayesian models
│   └── stan_models/             # Stan models folder
│       ├── survival_model.stan  # Survival analysis model
│       └── biomarker_model.stan # Longitudinal biomarker model
├── images/                      # Static images
│   └── study_design.svg         # Study design diagram
├── bibliography.bib             # References
└── styles.css                   # Custom CSS styles
```

## Getting Started

### Prerequisites

- [Quarto](https://quarto.org/docs/get-started/) (latest version)
- R (version 4.1.0 or higher)
  - rstanarm
  - tidyverse
  - plotly
  - DT
  - crosstalk
  - bayesplot
- Python (version 3.8 or higher)
  - pandas
  - numpy
  - matplotlib
  - seaborn

### Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/clinical-trial-analysis.git
   cd clinical-trial-analysis
   ```

2. Install R dependencies:
   ```r
   install.packages(c("rstanarm", "tidyverse", "plotly", "DT", "crosstalk", "bayesplot"))
   ```

3. Install Python dependencies:
   ```bash
   pip install pandas numpy matplotlib seaborn
   ```

### Running the Dashboard

To render the entire dashboard:

```bash
quarto render
```

To preview the dashboard with live updates:

```bash
quarto preview
```

## Reports

### Static Report

The static report (`static-report.qmd`) provides:

- Comprehensive methodological details
- Full Bayesian analysis of primary and secondary outcomes
- Biomarker trajectory analysis
- Subgroup analysis and model diagnostics

### Interactive Dashboard

The interactive dashboard (`interactive-report-static.qmd`) offers:

- Interactive filtering and visualization of outcomes
- Exploration of treatment effects across subgroups
- Probability threshold calculators for clinical significance
- Biomarker trajectory visualization

## Data Generation

The project includes scripts to generate synthetic trial data:

- `scripts/data_preparation.py` generates synthetic clinical trial data
- Use the parameter `n_patients` to control the sample size
- Use the parameter `treatment_effect` to control the effect size

## Customization

To customize the dashboard for your own clinical trial:

1. Replace the dummy data with your actual trial data
2. Update the Stan models to match your study design
3. Modify the visualizations for your specific outcomes
4. Update the study design diagram to match your trial

## Bayesian Models

The analysis uses the following Bayesian models:

- **Primary/Secondary Outcomes**: Linear models with treatment effects and covariate adjustment
- **Biomarker Analysis**: Longitudinal mixed-effects models
- **Survival Analysis**: Parametric Weibull survival models (for time-to-event data)

## License

[MIT License](LICENSE)

## References

See `bibliography.bib` for a complete list of references.

## Acknowledgments

- This project template was inspired by the CONSORT reporting guidelines for randomized trials
- Bayesian analysis approaches were adapted from the "Bayesian Data Analysis" textbook by Gelman et al.
