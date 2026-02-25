This is code accompanying the paper "A new mixture model for spatiotemporal exceedances with flexible tail dependence", available at <https://arxiv.org/abs/2602.13158> Details of the contents are provided below.

# Data folder

The data folder contains the following files:

1.  `SE_gcm_data.csv`, `SW_gcm_data.csv` - GCM data for the two regions

2.  `SE_obs_data.csv`, `SW_obs_data.csv` - observational data for the two regions

3.  `East_data.csv`, `West_data.csv` - bias correction based on CCA and QM for the two regions for 1951-2014

4.  `East_data_cv.csv`, `West_data_cv.csv` - bias correction based on CCA and QM for the two regions for 2001-2014

# Examples folder

This contains an `RMarkdown` file with a toy example demonstrating the methodology

# Results folder

Summaries that compile the output for the various methodologies, in `RData` format

# Files in the main folder of the repository

1. `00_distribution_shift.R` generates the plots for the exploratory data analysis for nonstationarity

2. `00_plotmap.R` generates Figure 2. (Note: the raw data for the plots are not currently uploaded)

3. `01-01_sim_data_generation.R` generates the data for the simulation study

4. `01-02_sim_correction_marginal.R`, `01-03_sim_correction_cross.R`, `01-04_sim_correction_SPECD.R` are the different density correction methods.
The second and third files correspond to SPECD1 and SPECD2 in the text. Results from the first method (marginal) are not presented in the text.

5. `01-05_compile_results_sim.R`, `02-01_compile_results_sim_competing.R` are used to compile the results (including metrics)
for the simulation studies. "Competing" refers to CCA and QM

6. `03-01_density_correction.R` fits SPECD and calibrates the GCM data for the entire time period of 1951-2014, `03-02_compile_results.R` assesses the calibration

7. `04-01_density-correction_validation.R` fits SPECD using GCM data for 1951-1000 and calibrates GCM data for 2001-2014, `04-02_compile_results_validation.R` assesses the calibration

8. `05-compile_results_competing.R` and `06-compile_results_competing_validation.R` assesses calibration using CCA and QM


