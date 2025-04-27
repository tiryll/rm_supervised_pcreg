# High Dimensional Prediction with Supervised Principal Components


## i) General Information


This project by [Fabian Schmidt](https://github.com/schmidtfabian) and me has been developed for the *Research Module in Econometrics & Statistics* at the University of Bonn. The goal of this work is to replicate the results of [Bair et. al (2006)](https://www.tandfonline.com/doi/abs/10.1198/016214505000000628), who propose the *supervised principal components*-method for prediction with high dimensional data, and evaluate its performance when compare to popular alternative methods.

In our simulations we show that, while the method can outperform alternative methods such as principal component regression (pcr) or partial least squares (pls) in terms of a lower Test-MSE, it performs worse than pcr and pls once certain technical assumptions are violated. 

Below you can find an overview of the repository structure.

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
