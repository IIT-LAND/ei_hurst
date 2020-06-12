# Code and data for EI Hurst paper

This repository has all the code and tidy data for the analyses in Trakoshis, Martínez-Cañada et al., Intrinsic excitation-inhibition imbalance affects medial prefrontal cortex differently in autistic men versus women. https://doi.org/10.1101/2020.01.16.909531

The code directory has all of the code for running the primary analyses. The analyses are split into 4 sections A, B, C, and D, and these are denoted at the beginning of each filename. Section A is the code for running in-silico modeling for the Gao model and the recurrent model. The code for the recurrent network model is located here: https://github.com/pablomc88/EEG_proxy_from_network_point_neurons.  Section B is for running in-vivo DREADD analyses. Section C is for running analyses on human rsfMRI data. Section D is for the gene expression enrichment analyses. Other code that these main scripts depend on are also in this directory.


## Requirements

* **AFNI** (https://afni.nimh.nih.gov/) - for preprocessing

* **FSL** (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) - for preprocessing

* **SPM12** (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) - for convolving simulated LFP to canonical HRF (`spm_conv.m`)


* **Python 2.7** and **3.6** or higher from the **Anaconda** distribution (https://www.anaconda.com/distribution/)

  + **neurodsp** (https://neurodsp-tools.github.io/neurodsp/) - for simulating LFP data based on 1/f and oscillations.

  + **nibabel** (https://nipy.org/nibabel/) - used by `dvars_se.py`.


* **R 3.5.1** or higher and **RStudio** (https://rstudio.com/)

  + **R libraries**

  **```c("here","ggplot2","patchwork","matlabr","nlme","reshape2","readxl","psych","ggseg","heplots","gplots")```**

  + **patchwork** (https://github.com/thomasp85/patchwork) - for putting together multiple ggplots.

  + **ggseg** (https://github.com/LCBC-UiO/ggseg) - for making plots of brains.


* **MATLAB R2018b** or higher (https://www.mathworks.com/)

  + **Brain Wavelet Toolbox** (http://www.brainwavelet.org/) - Used in preprocessing for the wavelet denoising step described by <a href="https://www.sciencedirect.com/science/article/pii/S1053811914001578?via%3Dihub">Patel et al., (2014)</a>.

  + **nonfractal** (https://github.com/wonsang/nonfractal) - Used for computation of H.

  + **EIslope** (https://github.com/voytekresearch/EISlope) - Used for running <a href="https://www.sciencedirect.com/science/article/abs/pii/S1053811917305621?via%3Dihub">Gao et al., (2017)</a> model that simulates LFP data based on manipulations of excitation and inhibition.

  + **PLS toolbox in MATLAB** (https://www.rotman-baycrest.on.ca/index.php?section=84) - Used for the PLS analysis on rsfMRI H data.



### In-silico modeling

* `A_insilico_0_analyze_recurrent_model.py` will run the steps for analyzing the recurrent model data. The data for this step is in the `data/recurrent_model` directory.

* `A_insilico_1_eisim.m` will run the model from <a href="https://www.sciencedirect.com/science/article/abs/pii/S1053811917305621?via%3Dihub">Gao et al., (2017)</a> that manipulates E:I ratio and then compute H on the simulated LFP data. H is computed with the `bfn_mfin_ml` function from `nonfractal`. This function saves data into the `data/gao_model` directory and is run as follows:

  ```
  EI_ratio = [2:0.2:6];
  MAKE_PLOT = 0;
  result = A_insilico_1_eisim(EI_ratio, MAKE_PLOT);
  ```

* `A_insilico_2_neural_ts_sim.py` will run the simulations to create LFP data based on 1/f slope. It can be run simply as shown below. This analysis requires python 3.6 or higher and utilizes the neurodsp library (https://neurodsp-tools.github.io/neurodsp/). This is primarily used to simulate the data that gets used in Supplementary Figure 1D.

  `python A_insilico_2_neural_ts_sim.py`

* `A_insilico_3_boldsim_neuralts_oof.m` will take the simulated LFP data from python in the previous step and will utilize it to compute H. This is primarily used for data going into Supplementary Figure 1D. It needs to be run as follows:

  `result = A_insilico_3_boldsim_neuralts_oof(0, 'oof');``


### In-vivo DREADD mouse rsfMRI analyses

* `B_invivo_1_DREADDpfc_excitation.Rmd` runs in RStudio and will call the `MATLAB` script `B_invivo_1_DREADDpfc_excitation.m` as the main code for running sliding window analyses on the DREADD excitation experiment. The remaining parts of the `B_invivo_1_DREADDpfc_excitation.Rmd` code will run the statistics and make plots. Running this `B_invivo_1_DREADDpfc_excitation.Rmd` will produce the `B_invivo_1_DREADDpfc_excitation.html` report found in the `code` directory.

* `B_invivo_2_DREADDpfc_silencing.Rmd` runs in RStudio and will call the `MATLAB` script `B_invivo_2_DREADDpfc_silencing.m` as the main code for running sliding window analyses on the DREADD silencing experiment. The remaining parts of the `B_invivo_2_DREADDpfc_silencing.Rmd` code will run the statistics and make plots. Running this `B_invivo_2_DREADDpfc_silencing.Rmd` will produce the `B_invivo_2_DREADDpfc_silencing.html` report found in the `code` directory.

### Autism rsfMRI analyses

* `C_1_preproc.sh` is a bash script that runs the preprocessing on the rsfMRI data. The main preprocessing script being called is `speedyppX.py`. This script calls many AFNI functions to do the main preprocessing. Note that `speedyppX.py` was written for `python 2.7` and may not work well in more recent versions of `python`. It also calls functions from the `Brain Wavelet Toolbox` to implement the wavelet denoising procedure described by <a href="https://www.sciencedirect.com/science/article/pii/S1053811914001578?via%3Dihub">Patel et al., (2014)</a>. After the preprocessing framewise displacement and DVARS are computed with `fd.py` and `dvars_se.py`. At the end of this bash script, it also calls a `MATLAB` function called `C_1a_parcEst.m` which will call `C_1b_parcellate.m` to parcellate the data by the HCP-MMP parcellation and then compute H based on those parcels.

* `C_2_AIMS_Hurst_Univariate.Rmd` runs in RStudio and does the main univariate analysis on H data. It runs all the stats for the sex-by-diagnosis interaction effect and other main effects, produces plots, and shows the results tables. It produces the `C_2_AIMS_Hurst_Univariate.html` report that can be found in the `code` directory.

* `C_3_AIMS_Hurst_PLS.Rmd` runs in RStudio and will run `C_3a_AIMS_Hurst_PLS.m` in `MATLAB` as the primary analysis of PLS on the rsfMRI data. The rest of the `C_3_AIMS_Hurst_PLS.Rmd` will make plots and produces the `C_3_AIMS_Hurst_PLS.html` report found in the `code` directory.

### Genomic analyses

* `D_asd_risk_genes_dht_de_overlap.Rmd` runs in RStudio and does the main enrichment analyses between autism-associated genes in different cell types and DHT DE genes.


### Data

Inside the `data` directory are subdirectories with the tidy data needed for the different aspects of the analyses. The names and filenames should be pretty self-explanatory and they get used at various points in the code.
