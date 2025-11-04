# MATLAB Codebase for "Assimilative Causal Inference"
This repository contains the MATLAB code for the paper "Assimilative Causal Inference", by [Marios Andreou](https://mariosandreou.short.gy/Homepage), [Nan Chen](https://people.math.wisc.edu/~nchen29/), and [Erik Bollt](https://webspace.clarkson.edu/~ebollt/), currently under review in _Nature Communications_.

The MATLAB codebase was written and tested on version R2024b.

# Script Files
The following MATLAB m-file (script) files are included in this repository (the numbered sections correspond to the ones found in the [arXiv](https://arxiv.org/abs/2505.14825) version of this work):
* dyad_interaction_model.m (Sections 3.1 and SI.4.1 - "A Nonlinear Dyad Model with Intermittent Extreme Events")
* noisy_predator_prey_model.m (Section SI.4.2 - "A Noisy Predator-Prey Model")
* El Niño-Southern Oscillation (ENSO) Case Study (Sections 3.2 and SI.4.3 - "A Stochastic Model Capturing the El Niño-Southern Oscillation (ENSO) Diversity"):
  + ENSO_model_cond_ACI_u_h_W_tau_unobs.m
  + ENSO_model_cond_ACI_u_unobs.m
  + ENSO_model_cond_ACI_h_W_unobs.m
  + ENSO_model_cond_ACI_T_C_unobs.m
  + ENSO_model_cond_ACI_tau_unobs.m
* Auxiliary files:
  + progress_bar.m
  + [simps.m](https://www.mathworks.com/matlabcentral/fileexchange/25754-simpson-s-rule-for-numerical-integration)
  + [legendUnq.m](https://www.mathworks.com/matlabcentral/fileexchange/67646-legendunq)

# Data
The user can use real observational data for the ENSO-related script files instead of using synthetic data as the observations, where the latter are formed by the values generated from the model simulation. These datasets can be found in the ["ENSO_DATA"](https://github.com/marandmath/ACI_code/tree/main/ENSO_DATA) subdirectory, with operational and implementation details for each dataset given in the docstrings of the corresponding m-files. The data are obtained from the following sources:
* [NCEP Global Ocean Data Assimilation System (GODAS) Reanalysis Dataset](https://www.esrl.noaa.gov/psd/data/gridded/data.godas.html): Behringer, D.W., M. Ji, and A. Leetmaa, 1998: An improved coupled model for ENSO prediction and implications for ocean initialization. Part I: The ocean data assimilation system. Mon. Wea. Rev., 126, 1013-1021.
* [NOAA Extended Reconstructed Sea Surface Temperature Version 5 (ERSST.v5) Dataset](https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html): Boyin Huang, Peter W. Thorne, Viva F. Banzon, Tim Boyer, Gennady Chepurin, Jay H. Lawrimore, Matthew J. Menne, Thomas M. Smith, Russell S. Vose, and Huai-Min Zhang (2017): NOAA Extended Reconstructed Sea Surface Temperature (ERSST), Version 5. NOAA National Centers for Environmental Information. doi:10.7289/V5T72FNM.
* [NCEP–NCAR Reanalysis 1 Project](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html): Kalnay et al.,The NCEP/NCAR 40-year reanalysis project, Bull. Amer. Meteor. Soc., 77, 437-470, 1996.

# Citing This Work

[Journal] [[arXiv](https://arxiv.org/abs/2505.14825)]

BibTeX Entry:
```
@article{
}
```

# License
This code is released under the MIT License. See the file [```LICENSE```](https://github.com/marandmath/ACI_code/blob/main/LICENSE) for copying permission.
