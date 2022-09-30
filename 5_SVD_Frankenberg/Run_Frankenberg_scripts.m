% Script to retrieve SIF in the far-red based on SVD.
%
% The Matlab scripts in this directory that conduct an SVD-based SIF 
% retrieval were written by Christian Frankenberg and his lab group 
% (fitAllRF.m, getEVs_loren.m). The concepts behind this approach are 
% described in Frankenberg et at (2018) [DOI: https://doi.org/10.1016/j.rse.2018.08.032]. 
% Please cite that paper if you use the Matlab scripts called by this script.
% This script calls functions that Christian Frankenberg wrote and shared,
% (and KC Cushman and Loren Albert edited, summer/fall 2020.) Loren Albert
% and KC Cushman wrote this script to run the other scripts for the various 
% spectral stray light scenarios.
%
% Original email from Christian Frankenberg:
% I put some quick&dirty fitting routine together and the impact seems marginal. 
% See attached screenshot. What I am fitting is basically SIF in the FR spectral 
% range (even assuming a flat spectral shape at the moment, so this can be improved). 
% I am also attaching some scripts, you can easily get fitted SIF like this:
% SIF_SL2 = fitAllFR("Spectra_for_Frankenberg/Vegetation_SL_Dmat_scenario3_SIFL2.hdr");
%
%
% For more information about SVD-based retrievals, readers are referred to the following articles:
% 
% Frankenberg, C., Köhler, P., Magney, T. S., Geier, S., Lawson, P., Schwochert, M., et al. 2018. The Chlorophyll Fluorescence Imaging Spectrometer (CFIS), mapping far red fluorescence from aircraft. Remote Sensing of Environment, 217, 523–536. https://doi.org/10.1016/j.rse.2018.08.032
% 
% Joiner, J., Guanter, L., Lindstrot, R., Voigt, M., Vasilkov, A. P., Middleton, E. M., Huemmrich, K. F., Yoshida, Y., and Frankenberg, C. 2013. Global monitoring of terrestrial chlorophyll fluorescence from moderate-spectral-resolution near-infrared satellite measurements: methodology, simulations, and application to GOME-2, Atmos. Meas. Tech., 6, 2803–2823, doi:10.5194/amt- 6-2803-2013.
% 
% Köhler, P., Guanter, L., Joiner, J. 2015. A linear method for the retrieval of sun-induced chlorophyll fluorescence from GOME-2 and SCIAMACHY data. Atmospheric Measurement Techniques, 8(6), 2589–2608. https://doi.org/10.5194/amt-8-2589-2015
% 
% Guanter, L., Rossini, M., Colombo, R., Meroni, M., Frankenberg, C., Lee, J-L., Joiner, J. 2013. Using field spectroscopy to assess the potential of statistical approaches for the retrieval of sun-induced chlorophyll fluorescence from ground and space. Remote Sensing of Environment 133, 52-61.
% 
% Guanter, L., Frankenberg, C., Dudhia, A., Lewis, P. E., Gómez-Dans, J., Kuze, A., et al. 2012. Remote Sensing of Environment. Remote Sensing of Environment, 121(C), 236–251. https://doi.org/10.1016/j.rse.2012.02.006
% 


%% Housekeeping
% clear all
% close all

% Loren's previous paths
% addpath('/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/Frankenberg_full physics_retrieval/From_Frankenberg')
%addpath('/Users/lalbert/Documents/Postdoc repos/ImagingSpectroscopy') % laptop
addpath('/Users/lalbert/Documents/Postdoc_RS_research/ImagingSpectroscopy') % IBES computer path
% Change directory
cd('/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/Spectra_output/');
% From Brown desktop computer:
% addpath('/Users/lalbert/Documents/Postdoc_RS_research/ImagingSpectroscopy')


% KC's computer paths
% addpath('C:/Users/cushmank/Documents/ImagingSpectroscopy');
% addpath('C:/Users/cushmank/Documents/ImagingSpectroscopy/Spectra_for_Frankenberg/');

%% Tests of functions and code
% new functions use two input files: a vegetation file and a solar file
% exampleVegFile = 'C:/Users/cushmank/Documents/ImagingSpectroscopy/Spectra_for_Frankenberg/Vegetation_no_SL_SIFL1';
% exampleSolarFile = 'C:/Users/cushmank/Documents/ImagingSpectroscopy/Spectra_for_Frankenberg/Subset_for_Frankenberg';
exampleVegFile = '/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/Spectra_output/Vegetation_SL_Dmat_scenario3_SIFL2';
exampleSolarFile = '/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/Spectra_output/Subset_for_Frankenberg';

SIF_test = fitAllFR(exampleVegFile,exampleSolarFile);

% % A test comparing mean veg spectra SIF and mean of 500 frames of SIF
% % from the vegetation cube:
% % Solar_no_SL_frames paired with Vegetation_no_SL_SIFl1.csv (SIF ~ 1mW @ 760 nm)
% Solar_no_SL_frames = 'Solar_no_SL_frames';
% Vegetation_no_SL_SIFl1 = readtable('Vegetation_no_SL_SIFl1.csv');
% SIF_no_SL_SIFl1 = fitAllFR(Vegetation_no_SL_SIFl1,Solar_no_SL_frames);
% % test with veg cube and reference cube (mean SIF here is the same as SIF using
% mean veg spectra only)
% Vegetation_no_SL_SIFl1 = 'Vegetation_no_SL_SIFl1';
% SIF_no_SL_SIFl1_cube = fitAllFR(Vegetation_no_SL_SIFl1,Solar_no_SL_frames);

%% No spectral stray light scenarios

% Fit vegetation with three levels of SIF
%Solar_no_SL_frames paired with Vegetation_no_SL_SIFl1.csv (SIF ~ 1mW @ 760 nm)
%Solar_no_SL_frames paired with Vegetation_no_SL_SIFl2
%Solar_no_SL_frames paired with Vegetation_no_SL_SIFl3.csv (SIF ~ 3mW @ 760 nm)

Solar_no_SL_frames = 'Solar_no_SL_frames';
Veg_scenario_no_SL = {'Vegetation_no_SL_SIFl1.csv', 'Vegetation_no_SL_SIFl2.csv','Vegetation_no_SL_SIFl3.csv'};
SIF_no_SL = cell(length(Veg_scenario_no_SL),1);
for j = 1:length(Veg_scenario_no_SL)
    veg = readtable(Veg_scenario_no_SL{j});
    SIF_no_SL{j} = fitAllFR(veg,Solar_no_SL_frames);
end

SIF_table_no_SL = table(Veg_scenario_no_SL', SIF_no_SL);

%% Very low level of spectral stray light scenario

% Solar_SL_Dmat_scenario4_frames paired with Vegetation_SL_Dmat_scenario4_SIFl1 (SIF ~ 1mW @ 760 nm)
% Solar_SL_Dmat_scenario4_frames paired with Vegetation_SL_Dmat_scenario4_SIFl2
% Solar_SL_Dmat_scenario4_frames paired with Vegetation_SL_Dmat_scenario4_SIFl3 (SIF ~ 3mW @ 760 nm)

Solar_SL_Dmat_scenario4_frames = 'Solar_SL_Dmat_scenario4_frames';
Veg_Dmat_scenario4 = {'Vegetation_SL_Dmat_scenario4_SIFl1.csv', 'Vegetation_SL_Dmat_scenario4_SIFl2.csv','Vegetation_SL_Dmat_scenario4_SIFl3.csv'};
SIF_Dmat_scenario4 = cell(length(Veg_Dmat_scenario4),1);
for j = 1:length(Veg_Dmat_scenario4)
    veg = readtable(Veg_Dmat_scenario4{j});
    SIF_Dmat_scenario4{j} = fitAllFR(veg,Solar_SL_Dmat_scenario4_frames);
end

SIF_table_Dmat_scenario4 = table(Veg_Dmat_scenario4', SIF_Dmat_scenario4);

%% Low level of spectral stray light scenario

% Solar_SL_Dmat_scenario3_frames paired with Vegetation_SL_Dmat_scenario3_SIFl1 (SIF ~ 1mW @ 760 nm)
% Solar_SL_Dmat_scenario3_frames paired with Vegetation_SL_Dmat_scenario3_SIFl2 
% Solar_SL_Dmat_scenario3_frames paired with Vegetation_SL_Dmat_scenario3_SIFl3 (SIF ~ 3mW @ 760 nm)

Solar_SL_Dmat_scenario3_frames = 'Solar_SL_Dmat_scenario3_frames';
Veg_Dmat_scenario3 = {'Vegetation_SL_Dmat_scenario3_SIFl1.csv', 'Vegetation_SL_Dmat_scenario3_SIFl2.csv','Vegetation_SL_Dmat_scenario3_SIFl3.csv'};
SIF_Dmat_scenario3 = cell(length(Veg_Dmat_scenario3),1);
for j = 1:length(Veg_Dmat_scenario3)
    veg = readtable(Veg_Dmat_scenario3{j});
    SIF_Dmat_scenario3{j} = fitAllFR(veg,Solar_SL_Dmat_scenario3_frames);
end

SIF_table_Dmat_scenario3 = table(Veg_Dmat_scenario3', SIF_Dmat_scenario3);

%% Medium level of spectral stray light scenario

% Solar_SL_Dmat_scenario2_frames paired with Vegetation_SL_Dmat_scenario2_SIFl1 (SIF ~ 1mW @ 760 nm)
% Solar_SL_Dmat_scenario2_frames paired with Vegetation_SL_Dmat_scenario2_SIFl2
% Solar_SL_Dmat_scenario2_frames paired with Vegetation_SL_Dmat_scenario2_SIFl3 (SIF ~ 3mW @ 760 nm)

Solar_SL_Dmat_scenario2_frames = 'Solar_SL_Dmat_scenario2_frames';
Veg_Dmat_scenario2 = {'Vegetation_SL_Dmat_scenario2_SIFl1.csv', 'Vegetation_SL_Dmat_scenario2_SIFl2.csv','Vegetation_SL_Dmat_scenario2_SIFl3.csv'};
SIF_Dmat_scenario2 = cell(length(Veg_Dmat_scenario2),1);
for j = 1:length(Veg_Dmat_scenario2)
    veg = readtable(Veg_Dmat_scenario2{j});
    SIF_Dmat_scenario2{j} = fitAllFR(veg,Solar_SL_Dmat_scenario2_frames);
end

SIF_table_Dmat_scenario2 = table(Veg_Dmat_scenario2', SIF_Dmat_scenario2);

%% Very high level of spectral stray light scenario

% Solar_SL_Dmat_scenario1_frames paired with Vegetation_SL_Dmat_scenario1_SIFl1 (SIF ~ 1mW @ 760 nm)
% Solar_SL_Dmat_scenario1_frames paired with Vegetation_SL_Dmat_scenario1_SIFl2
% Solar_SL_Dmat_scenario1_frames paired with Vegetation_SL_Dmat_scenario1_SIFl3 (SIF ~ 3mW @ 760 nm)

Solar_SL_Dmat_scenario1_frames = 'Solar_SL_Dmat_scenario1_frames';
Veg_Dmat_scenario1 = {'Vegetation_SL_Dmat_scenario1_SIFl1.csv', 'Vegetation_SL_Dmat_scenario1_SIFl2.csv','Vegetation_SL_Dmat_scenario1_SIFl3.csv'};
SIF_Dmat_scenario1 = cell(length(Veg_Dmat_scenario1),1);
for j = 1:length(Veg_Dmat_scenario1)
    veg = readtable(Veg_Dmat_scenario1{j});
    SIF_Dmat_scenario1{j} = fitAllFR(veg,Solar_SL_Dmat_scenario1_frames);
end

SIF_table_Dmat_scenario1 = table(Veg_Dmat_scenario1', SIF_Dmat_scenario1);

%% Combine tables of SIF results and export

SIF_table_no_SL.Properties.VariableNames = {'Scenario','Retrieved_SIF'};
SIF_table_Dmat_scenario4.Properties.VariableNames = {'Scenario','Retrieved_SIF'};
SIF_table_Dmat_scenario3.Properties.VariableNames = {'Scenario','Retrieved_SIF'};
SIF_table_Dmat_scenario2.Properties.VariableNames = {'Scenario','Retrieved_SIF'};
SIF_table_Dmat_scenario1.Properties.VariableNames = {'Scenario','Retrieved_SIF'};
SIF_table = [SIF_table_no_SL; SIF_table_Dmat_scenario4; SIF_table_Dmat_scenario3; SIF_table_Dmat_scenario2; SIF_table_Dmat_scenario1];
writetable(SIF_table, 'SIF_SVD_farRed.csv')


