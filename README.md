# High Dimensional Prediction with Supervised Principal Components


## i) General Information

This project, developed by [Fabian Schmidt](https://github.com/schmidtfabian) and myself, was created for the Research Module in Econometrics & Statistics at the University of Bonn. Our objective is to evaluate analytically and through simulations the performance of the supervised principal components (SPC) method introduced by [Bair et. al (2006)](https://www.tandfonline.com/doi/abs/10.1198/016214505000000628) for prediction with high-dimensional data.

Through our simulations, we demonstrate that while the SPC method can outperform popular alternative techniques such as Principal Component Regression (PCR) and Partial Least Squares (PLS) in terms of achieving a lower test mean squared predition error (Test-MSPE), its performance deteriorates once certain technical assumptions are violated.

Below, you'll find an overview of the repository structure.

---------------------
## ii) Overview
---------------------

    rm_supervised_pcreg
    │   README.md
    │   environment.R
    │
    └───functions
    │   │   functions_pcr.R
    │   │   functions_pls.R
    │   │   MC_simulation_functions.R
    │   
    └───output
    │   │ 
    │   └───data
    │       │   (.RDS data for figures)
    │   └───figures
    │       │   (13 .pdf figures)
    │   └───tables
    │       │   (12 LaTeX tables)
    │
    └───plots
    │   │   boxplots_simulation_1.R
    │   │   boxplots_simulation_2.R
    │   │   boxplots_simulation_3.R
    │   │   heatmap.R
    │  
    └───simulations
    │   │   MC_simulation_1.R
    │   │   MC_simulation_2.R
    │   │   MC_simulation_3.R
    │   │   MC_simulation_Appendix_GF.R
