#!/usr/bin/env python
# coding: utf-8

"""
Data preparation script for clinical trial analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Create data directory if it doesn't exist
os.makedirs('../data', exist_ok=True)

# Set random seed for reproducibility
np.random.seed(42)

def generate_trial_data(n_patients=200, treatment_effect=1.5):
    """
    Generate synthetic trial data
    
    Parameters:
    -----------
    n_patients : int
        Number of patients in the trial
    treatment_effect : float
        Effect size for the treatment
        
    Returns:
    --------
    pandas.DataFrame
        Synthetic trial data
    """
    # Generate treatment assignments
    treatment = np.repeat(['Treatment', 'Control'], n_patients // 2)
    
    # Generate demographic variables
    age = np.random.normal(65, 10, n_patients).round().astype(int)
    sex = np.random.choice(['Male', 'Female'], n_patients, p=[0.55, 0.45])
    
    # Generate baseline risk score
    baseline_risk = np.random.normal(5, 1, n_patients)
    
    # Add treatment effect
    treatment_effect_vec = np.where(
        treatment == 'Treatment',
        treatment_effect + np.random.normal(0, 0.5, n_patients // 2),
        np.random.normal(0, 0.5, n_patients // 2)
    )
    treatment_effect_vec = np.concatenate([treatment_effect_vec, treatment_effect_vec])
    
    # Generate outcomes
    primary_outcome = baseline_risk + treatment_effect_vec + \
                     0.05 * (age - 65) + \
                     0.5 * (sex == 'Female') + \
                     np.random.normal(0, 1, n_patients)
    
    secondary_outcome = 0.6 * primary_outcome + \
                        np.random.normal(0, 0.8, n_patients)
    
    # Create DataFrame
    trial_data = pd.DataFrame({
        'patient_id': range(1, n_patients + 1),
        'treatment': treatment,
        'age': age,
        'sex': sex,
        'baseline_risk': baseline_risk,
        'primary_outcome': primary_outcome,
        'secondary_outcome': secondary_outcome
    })
    
    return trial_data

def generate_biomarker_data(trial_data, n_timepoints=3):
    """
    Generate longitudinal biomarker data
    
    Parameters:
    -----------
    trial_data : pandas.DataFrame
        Trial data with patient information
    n_timepoints : int
        Number of timepoints to generate
        
    Returns:
    --------
    pandas.DataFrame
        Longitudinal biomarker data
    """
    patient_ids = []
    timepoints = []
    treatments = []
    biomarker_values = []
    
    timepoint_labels = ['Baseline', 'Week 6', 'Week 12']
    
    # Patient-specific random effects
    n_patients = len(trial_data)
    patient_intercepts = np.random.normal(0, 5, n_patients)
    patient_slopes = np.random.normal(0, 1, n_patients)
    
    for i, row in trial_data.iterrows():
        patient_id = row['patient_id']
        treatment = row['treatment']
        
        # Patient-specific intercept
        baseline = 50 + patient_intercepts[i]
        
        for t in range(n_timepoints):
            # Time effect
            time_effect = 5 * t
            
            # Treatment effect (increases over time)
            if treatment == 'Treatment':
                treatment_effect = 3 * t
            else:
                treatment_effect = 0
                
            # Patient-specific slope
            patient_effect = t * patient_slopes[i]
            
            # Add noise
            noise = np.random.normal(0, 2)
            
            # Calculate biomarker value
            biomarker = baseline + time_effect + treatment_effect + patient_effect + noise
            
            # Append to lists
            patient_ids.append(patient_id)
            timepoints.append(timepoint_labels[t])
            treatments.append(treatment)
            biomarker_values.append(biomarker)
    
    # Create DataFrame
    biomarker_data = pd.DataFrame({
        'patient_id': patient_ids,
        'treatment': treatments,
        'timepoint': timepoints,
        'biomarker_value': biomarker_values
    })
    
    return biomarker_data

if __name__ == "__main__":
    # Generate trial data
    trial_data = generate_trial_data(n_patients=200, treatment_effect=1.7)
    
    # Generate biomarker data
    biomarker_data = generate_biomarker_data(trial_data)
    
    # Save data
    trial_data.to_csv('../data/trial_data.csv', index=False)
    biomarker_data.to_csv('../data/biomarker_data.csv', index=False)
    
    # Print summary statistics
    print("Trial Data Summary:")
    print(trial_data.groupby('treatment')[['primary_outcome', 'secondary_outcome']].agg(['mean', 'std']))
    
    print("\nBiomarker Data Summary:")
    print(biomarker_data.groupby(['treatment', 'timepoint'])['biomarker_value'].agg(['mean', 'std']))
    
    # Create exploratory plot
    plt.figure(figsize=(10, 6))
    
    sns.boxplot(x='treatment', y='primary_outcome', data=trial_data)
    plt.title('Primary Outcome by Treatment')
    plt.tight_layout()
    
    plt.savefig('../images/primary_outcome_boxplot.png', dpi=300)
    plt.close()
    
    # Create biomarker trajectory plot
    plt.figure(figsize=(10, 6))
    
    sns.lineplot(x='timepoint', y='biomarker_value', hue='treatment', data=biomarker_data, marker='o', err_style='bars', ci=95)
    
    plt.title('Biomarker Trajectory')
    plt.tight_layout()
    
    plt.savefig('../images/biomarker_trajectory.png', dpi=300)
    plt.close()
    
    print("Data preparation complete.")
    