This repository contains all the source code used in the study "The Influence of Intercalated Disc Nanostructure on Local Ionic Currents and Cardiac Conduction."

Repository Structure

Main Branch
Contains all MATLAB scripts used for electrophysiological simulations and data analysis presented in the paper.

mesh_data/ Folder
Includes the complete set of mesh files used in the simulations:

384/: The 384 meshes analyzed in the main text.

25_femdata_random/: 20 randomly generated meshes used for MLP model validation.

Additional FEM data and mesh partitioning files for reproducibility.

NN_sensitivity_analysis/ Folder
Contains the neural network–based sensitivity analysis scripts:

run_nn_regression_final.py: Executes the full sensitivity analysis pipeline, generating the loss table, regression plots, and sensitivity heatmap.

run_nn_regression_ID_heterogeneity.ipynb: Provides the foundational data structures and analysis framework used in the study.
