### Data information ###

This repository includes Jupyter notebooks & Python scripts for the analysis of HighResMIP models. Model composites of cold North Atlantic SSTs are examined, which are compared to the ERA5 reanalysis product used as reference. The composite analysis is based on the study by Krüger et al. (2023). https://a.tellusjournals.se/articles/10.16993/tellusa.3235  

The repository includes three Jupyter Notebooks:

* 1) *HighResMIP_North_Atlantic_biases.ipynb* - This notebook calculates the North Atlantic biases of SST and SLHF and absolute bias differences. Significance is calculated based on a sample T-test. The results are shown in Fig. 1,2 and Fig. S2, S3.

* 2) *SLHF_ERA5_OA_Flux_diff.ipynb* - Systematic differences of surface latent heat flux based on data sets from ERA5 reanalysis and OA Flux. Systematic differences for each season are now shown in Fig. S1.

* 3) *HighResMIP_ERA5_cold_NASST_composites.ipynb* - This notebook contains the composite analysis for the HighResMIP models and ERA5. Additionally, the climatological summer mean of Z300 is computed for all models individually. Results are shown in Fig. 3, 4, Fig. S4, S5, S6, S7, S8, S9, S10, S11, S12.


Further, the respository contains five python scripts:

* 1) *HighResMIP_bootstrap_map_composites_SST_T2m.py*  
* 2) *HighResMIP_bootstrap_map_composites_Z300.py*
* 3) *HighResMIP_bootstrap_map_composites_SLHF.py*
* 4) *HighResMIP_bootstrap_map_composites_pr.py*
* 5) *HighResMIP_bootstrap_map_composites_egr.py*

These python scripts are used for the bootstrap analysis used for the significance of the anomalies of the five variables SST&T2m, Z300, SLHF, precipitation, Eady Growth Rate shown in Fig. 3, Fig. 4c,d and Fig. S4, S5, S6, S7, S8, S9, S10, S11.
