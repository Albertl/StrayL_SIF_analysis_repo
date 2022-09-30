---
editor_options: 
  markdown: 
    wrap: 72
---

# StrayL_SIF_analysis_repo

## Overview

This repository contains analysis and plotting scripts (some in R some
in Matlab) to reproduce the spectral stray light simulations, FLD
retrieval, 3FLD retrieval, and SVD-based retrieval presented in:

Albert et al. in "Sensitivity of solar-induced fluorescence to spectral
stray light in high resolution imaging spectroscopy"

The scripts are numbered to show the order for running, and include
descriptions in the header. Full information on the methods used in this
study can be found in the paper listed above. In brief, the scripts are:

## 1_make_SLC_Omac_lines.R

R script to make plots showing the effect of performing the stray light
correction to laser lines, and to vegetation and spectralon (reference
target). The vegetation and spectralon spectra are output for use in
other scripts: ./Outputs/OmacData_Rad_W and ./Outputs/OmacData_SLC_Rad_W

## 2_make_stray_light_simulations_v2.R

Script for simulations to examine the effect of stray light on SIF,
including sFLD and 3FLD retrievals. Outputs include:

1.  Differences between with and w/o stray light for 3FL:
    3FLD_SIFlX_stray_no_stray_differences.csv (where X is a level from 1
    to 3, see script header.)

2.  An interpolated version of a \~99% reflectance spectrum for use in
    Matlab scripts: HL21-24_Spec_C1HL21_interpolated.csv

3.  A matrix of wavelengths values corresponding to minima of spectral
    features (for use in later scripts for making table in
    SFMresults_summary.R and make_stray_light_error_figure.R):
    Estimated_window_mins_ind_wl.csv

4.  SIF spectra (with name customized to show the scaling for each of
    the 3 levels defined early in the script): SIF_SIFlX.csv (where X is
    a level from 1 to 3, see script header.)

## 3_SIF_vs_stray_light_magnitude_fig.R

Script to make Figure 1, showing magnitude of spectral stray light on
vegetation and reference spectra. Script also includes a version of the
figure including SIF spectra on same plot.

## 4_SFM_results_summary.R

Script examining and organizing SFM retrieval results (SFM retrieval was
conducted seperately by Luis Alonso and Pablo Morcillo Pallarés).
Produces SFM figures for results section and supplement.

## 5_SVD_Frankenberg (directory)

The Matlab scripts in this directory that conduct an SVD-based SIF
retrieval were written by Christian Frankenberg and his lab group
(fitAllRF.m, getEVs_loren.m, read_hdr.m). The concepts behind this
approach are described in Frankenberg et at (2018) [DOI:
<https://doi.org/10.1016/j.rse.2018.08.032>]. Please cite that paper if
you use the Matlab scripts in this directory. KC Cushman and Loren
Albert made edits so that the code ran on paired reference and target
spectra for the spectral stray light sensitivity analysis (amongst other
minor changes), and wrote the 'master' script that runs the SVD-based
SIF retrieval for the various spectral stray light scenarios. Note that
the Matlab File Exchange script read_envihdr is also required:

Jaroslaw Tuszynski (2022). read_envihdr
(<https://www.mathworks.com/matlabcentral/fileexchange/38500-read_envihdr>),
MATLAB Central File Exchange. Retrieved September 30, 2022.

The R script 5_SVDresults_summary was written by Loren Albert for
plotting and organizing the SVD-based retrieval results. For more
information about SVD-based retrievals, readers are referred to the
following articles:

Frankenberg, C., Köhler, P., Magney, T. S., Geier, S., Lawson, P.,
Schwochert, M., et al. 2018. The Chlorophyll Fluorescence Imaging
Spectrometer (CFIS), mapping far red fluorescence from aircraft. Remote
Sensing of Environment, 217, 523--536.
<https://doi.org/10.1016/j.rse.2018.08.032>

Joiner, J., Guanter, L., Lindstrot, R., Voigt, M., Vasilkov, A. P.,
Middleton, E. M., Huemmrich, K. F., Yoshida, Y., and Frankenberg, C.
2013. Global monitoring of terrestrial chlorophyll fluorescence from
moderate-spectral-resolution near-infrared satellite measurements:
methodology, simulations, and application to GOME-2, Atmos. Meas. Tech.,
6, 2803--2823, <doi:10.5194/amt-> 6-2803-2013.

Köhler, P., Guanter, L., Joiner, J. 2015. A linear method for the
retrieval of sun-induced chlorophyll fluorescence from GOME-2 and
SCIAMACHY data. Atmospheric Measurement Techniques, 8(6), 2589--2608.
<https://doi.org/10.5194/amt-8-2589-2015>

Guanter, L., Rossini, M., Colombo, R., Meroni, M., Frankenberg, C., Lee,
J-L., Joiner, J. 2013. Using field spectroscopy to assess the potential
of statistical approaches for the retrieval of sun-induced chlorophyll
fluorescence from ground and space. Remote Sensing of Environment 133,
52-61.

Guanter, L., Frankenberg, C., Dudhia, A., Lewis, P. E., Gómez-Dans, J.,
Kuze, A., et al. 2012. Remote Sensing of Environment. Remote Sensing of
Environment, 121(C), 236--251.
<https://doi.org/10.1016/j.rse.2012.02.006>

## 6_make_stray_light_error_figure.R

This script makes a figure showing the error attributable to stray light
(in radiance units) for all retrievals examined in the paper associated
with this repository. The retrievals are imported from the outputs ofthe
following scripts:

1.  SFMresults_summary.R for "SFM_stray_no_stray_differences.csv"

2.  SVDresults_summary.R for "SVD_stray_no_stray_differences.csv"

3.  make_stray_light_simulations_v2.R for threeFLD_SIFl1_diffs,
    threeFLD_SIFl2_diffs and threeFLD_SIFl3_diffs

## Stray_light_analysis_functions.R R

Script with all general functions for the instrument used in the case
study (including function for stray light correction to data).

------------------------------------------------------------------------
