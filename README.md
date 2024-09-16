This repository includes python scripts & Jupyter notebooks for the analysis of HighResMIP models. Model composites of cold North Atlantic SSTs are examined, which are compared to the ERA5 reanalysis product used as reference. The composite analysis is based on the study by Krüger et al. (2023). https://a.tellusjournals.se/articles/10.16993/tellusa.3235  

The repository includes four Jupyter Notebooks:

* 1) *HighResMIP_North_Atlantic_biases.ipynb* - This notebook calculates the North Atlantic biases of SST and SLHF and absolute bias differences. Significance is calculated using a sample T-test. The results are shown in Fig. 1,2 and Fig. S1,S2.

* 2) *HighResMIP_ERA5_cold_NASST_composites.ipynb* - This notebook contains the composite analysis for the HighResMIP models and ERA5. Results are shown in Fig. 3, 4, Fig. S3, S4, S5, S6.

* 3) *HighResMIP_ERA5_AMV.ipynb* - This notebook includes the calculation of the AMV of the HighResMIP models and their relationship to the composites of the European T2m anomalies. Results are shown in Fig. 5 and Fig. S7, S8.

* 4) *HadISST_AMV.ipynb* - This notebook calculates the AMV based on the HadISST data set. The results are used for Fig. 5.


Further, the respository contains four python scripts:

* 1) *HighResMIP_bootstrap_map_composites_SST_T2m.py*  
* 2) *HighResMIP_bootstrap_map_composites_Z300.py*
* 3) *HighResMIP_bootstrap_map_composites_SLHF.py*
* 4) *HighResMIP_bootstrap_map_composites_pr.py*

These python scripts are used for the bootstrap analysis used for the significance of the map anomalies of the four variables SST&T2m, Z300, SLHF, pr shown in Fig. 3, Fig. 4c,d and Fig. S4.
