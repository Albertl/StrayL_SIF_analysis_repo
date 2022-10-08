### Script for simulations to examine the effect of stray light on SIF. Need to define objects under 'housekeeping' section below.
### This second version (v2) has options for different levels of SIF, and can be re-run with different SIF levels. Note that 
### warnings are produced regarding NaNs from log functions and logarithmic plots. This is not unexpected, and results from 
### the dark noise subtraction step of data treatment (when sensors are used in low light, dark noise subtraction results in, 
### essentially, noise minus noise, and some small negative numbers.)
###
### Contact: Loren Albert, lalbert@email.arizona.edu; loren.albert@mail.wvu.edu
###
### Copyright (C) 2022  Loren Albert. See "LICENSE"
### 
### Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
### and associated documentation files (the "Software"), to deal in the Software without restriction, 
### including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
### and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
### subject to the following conditions:
###
### The above copyright notice and this permission notice shall be included in all copies or substantial 
### portions of the Software.
###
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
### NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
### IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
### WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
### SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



######## Housekeeping #########################
## Remove workspace objects if necessary
# rm(list=ls())

## Select the level of SIF emission (for sensitivity of results to different levels of SIF)
# 1 = ~1 mW/m2/nm/sr @760 nm
# 2 = ~2 mW/m2/nm/sr @760 nm  # Main text results
# 3 = ~3 mW/m2/nm/sr @760 nm
SIF_level <- 2

##### INPUTS ##################################
## load library to read SIF spectrometer image
library(caTools)
library(raster)
library(abind)
library(ggplot2)
library(htmlTable)
require(MASS)
library(viridis)

## define path to data files used
path_data <- './Inputs/'

## define path to figure outputs
path_figs <- "./Figures/"

## define path for spectra output and other outputs
path_spectra <- "./Outputs/"
path_outputs <- "./Outputs/"


## READ IN DATA AND FUNCTIONS

## source script with all general functions for this instrument (including function for stray light correction to data)
source("Stray_light_analysis_functions.R")

## C matrix (used to multiply raw data for stray light correction)
Cmatrix <- read.csv(paste(path_data,"SLC_mat.csv",sep=""),header=F)
Cmatrix <- as.matrix(Cmatrix)
## Wavelength values for each spectral pixel
Wavelengths <- read.csv(paste(path_data,"predicted_spectral_pixels_2018.csv",sep=""))
## Laser line data (just for plot explaining stray light)
LineData <- read.csv(paste(path_data,"Bracket_10ms_1000ms_scaled_by_ms_with_correction.csv",sep=""))


## SIF spectrometer (processed by output from our version 1 pipeline (Drone2Radiance.m), followed by
# Subset_Radiance_for_Frankenberg.m, in which 5 pixels on Spectralon
# were averaged, units converted, and spectralon spectra corrected with RELAB scan, while
# maintaining the 500 frames for SVD with temporal variation.
# Thus dimensions are 500 frames x 1 pixel x 2160 bands. Units are already W/m2/sr/nm
SolarData_SLC_frames <- read.ENVI(paste(path_data, "Subset_for_Frankenberg",sep=""))
SolarData_SLC_frames <- drop(SolarData_SLC_frames)

## ASD FluoWat data
OmacData_ASD_SIF <- read.csv(paste(path_data,"OMAC_FluoWatt_SIF.csv",sep=""))
OmacData_ASD_Leaf <- read.csv(paste(path_data,"OMAC_FluoWatt_Leaf.csv",sep=""))
OmacData_ASD_Solar <- read.csv(paste(path_data,"OMAC_FluoWatt_Solar.csv",sep=""))


## Read in RELAB scan for white reference target
col_names <- c("wavelength","mean","stddev")
Refl_99 <- ez.read(paste(path_data, "RELAB_spectralon_scan/HL21-24_Spec/C1HL21.ASC", sep = ""), 
                   skip.rows = c(1,233:247), sep = "", header = FALSE, col.names = col_names)




##### CALIBRATIONS AND CORRECTIONS ###############################
# Note that the input data "SolarData_SLC_frames" already had
# dark subtraction, stray light correction, and is already in radiance units.


# # STRAY LIGHT CORRECTION         

# apply stray light correction to laser line data
LineData_SLC <- SLC_fun_colwise(data_mat = LineData, C_mat = Cmatrix) # takes a few seconds to run


# Correct for the fact that spectralon is not 100% reflective and not perfectly flat.
Refl_99_interp <- approx(x = Refl_99$wavelength, y = Refl_99$mean, xout = Wavelengths$spec_pixel_pred_no710)
Refl_99_interp <- as.numeric(Refl_99_interp$y)
# SolarData_SLC_corrected <- SolarData_SLC / Refl_99_interp
# SolarData_SLC <- SolarData_SLC_corrected


##### SUBSET ASD LEAF REFLECTANCE, RESAMPLE, AND CALCULATE LEAF REFLECTANCE ##################################

## choose x-values
xVals <- Wavelengths[,2]

# Subset to red and infrared regions
SIF_ASD_ind <- which(OmacData_ASD_SIF$Wavelength >= 669 & OmacData_ASD_SIF$Wavelength <= 781)
SIF_ASD_wave_subset <- OmacData_ASD_SIF$Wavelength[SIF_ASD_ind]
SIF_ASD_subset <- OmacData_ASD_SIF$Rad_mean[SIF_ASD_ind]
Leaf_ASD_subset <- OmacData_ASD_Leaf$Rad_mean[SIF_ASD_ind]
Solar_ASD_subset <- OmacData_ASD_Solar$Rad_mean[SIF_ASD_ind]

# Calculate leaf reflectance
ASD_Leaf_R <- Leaf_ASD_subset/Solar_ASD_subset
plot(ASD_Leaf_R, type = "l")
abline(v = c(89,102), col = "red")

# Leaf reflectance has a bump near the oxygen A band. I think this was because FluoWat measurements were done late in the day (rapidly changing solar angle) and
# the leaf was measured before the white reference (so it makes sense that as the sun when down the atmosphere absorbed more light in the oxygen band for the late 
# white reference measurements). Leaf reflectance spectra usually don't have a bump here (e.g. see https://www.researchgate.net/publication/39569434_Using_Hyperspectral_Remote_Sensing_to_map_grape_quality_in_%27Tempranillo%27_Vineyards_affected_by_Iron_deficiency_Chlorosis)
# So linearly interpolate to remove the bump. Resampling to the high resolution spectrometer is done in this same step.
Leaf_R_gap <- ASD_Leaf_R
Leaf_R_gap[89:102] <- NA
Leaf_R_resampled <- approx(x = SIF_ASD_wave_subset, y = Leaf_R_gap, xout = xVals)
plot(Leaf_R_resampled, type = "l")

# Define reflectance shape (earlier version included an if statement for constant shape, but now I think that is unnecessary)
leaf_R <- Leaf_R_resampled$y


## For each frame, calculate leaf reflected radiance with L = r * E for no stray light scenario
L_veg_frames <- apply(SolarData_SLC_frames, 1, function(x) leaf_R*x)

# Calculate mean of leaf reflected radiance (for scaling SIF later)
L_veg_frames_mean <- apply(L_veg_frames, 1, mean)
plot(xVals, L_veg_frames_mean, type = "l")
# Looks similar to Meroni et al 2009 Fig. 1


##### RESAMPLE AND RESCALE SIF SPECTRUM, THEN  ADD SIF TO THE VEGETATION SPECTRA (WITHOUT STRAY LIGHT SPECTRA ADDED) ###################################

# Resample SIF spectrum
SIF_resampled <- approx(x = SIF_ASD_wave_subset, y = SIF_ASD_subset, xout = xVals) 
plot(SIF_resampled$x, SIF_resampled$y)

# Units should be W/m2/sr/nm, but since the lighting conditions were different, and since the Fluowat might reduce light as well, 
# re-scale it so SIF is 2% of signal in the far-red. The left shoulder for FLD is often 758 nm, so use 758 as a reference for scaling the SIF shape
#SIF_scaler <- L_veg[1723] * 0.02

# First make SIF = 1 at 758 nm (to make next steps easier)
SIF_scaled = (1 / SIF_resampled$y[1723]) * (SIF_resampled$y)

# Now make SIF_2perc[1723] / L_veg_frames_mean[1723] = 0.02 (or whatever the level of SIF is set by SIF_level)
if(SIF_level==1){
  SIF_2perc <- SIF_scaled * (L_veg_frames_mean[1723] * 0.007) # SIF_2perc at O2A minimum is SIF_2perc[1772] = 0.001061908
}else if(SIF_level==2){
  SIF_2perc <- SIF_scaled * (L_veg_frames_mean[1723] * 0.014) # SIF_2perc at O2A minimum is SIF_2perc[1772] = 0.00212385
}else if(SIF_level==3){
  SIF_2perc <- SIF_scaled * (L_veg_frames_mean[1723] * 0.02)  # SIF_2perc at O2A minimum is SIF_2perc[1772] = 0.00302789
}



# Add SIF to each frame of vegetation spectra with no stray light
L_veg_SIF_frames <- L_veg_frames + SIF_2perc



##### SIMULATE SEVERAL LEVELS OF STRAY LIGHT ######################################################################################################

# Zong et al 2006 says: "The level of this spectral stray light is of the order of 10^-3 to 10^-5 for the measurements of a monocromatic source or 
# of the order of 10^-1 to 10^-3 for the measurement of a broadband source, depending on the quality of the spectrograph"
# So since sunlight is a broadband source:
# Use 10^-1 for a 'poor quality' instrument measuring a broadband light source
# Use 10^-3 for a 'high quality' instrument measuring a broadband light source
# Use 10^-5 for a high quality instrument that has had stray light corrected (and reduced 2 orders of magnitude)


### FROM D MATRIX
# Try simulating stray light with a matrix
# Want simulated 'D' matrix in Zong et al 2006, and an identity matrix 'I'
# Then simulate measurements in the presence of stray light as Y_meas = (I+D)*Y_IB
# Which is equivalent to Y_meas = A * Y_IB

## Simulate 'D' matrix based on equations in Zong et al 2006 
I_mat <- diag(x = 1, nrow=2160, ncol=2160)
A_Mat <- ginv(Cmatrix)
D_Mat <- A_Mat - I_mat

# Multiply stray light matrix so that stray light will be 10^-1, 10^-2, 10^-3, or 10^-4 compared to the broadband signal
# Original D matrix produces stray light that is 10^-3 across most of the spectrum (but is 10^-2 in O2-A region)
A_Mat_40x <- (D_Mat * 40) + I_mat         # Equivalent to (D_Mat * 10^1.60206) + I_mat
A_Mat_4x <- (D_Mat * 4) + I_mat           # Equivalent to (D_Mat * 10^0.60206) + I_mat
A_Mat_0.4x <- (D_Mat * 0.4) + I_mat       # Equivalent to (D_Mat * 10^-0.39794) + I_mat
A_Mat_0.04x <- (D_Mat * 0.04) + I_mat     # Equivalent to (D_Mat * 10^-1.39794) + I_mat

# Simulate stray light for vegetation spectrum (that includes simulated SIF)
# (Simulation of stray light and the addition to the 'measured' spectrum happen all in one step since this uses matrices)
# Each of the 500 frames
veg_SL_Dmat_1_frames <- simulate_SL_fun_colwise(data_mat = L_veg_SIF_frames, A_mat = A_Mat_40x)
veg_SL_Dmat_2_frames <- simulate_SL_fun_colwise(data_mat = L_veg_SIF_frames, A_mat = A_Mat_4x)
veg_SL_Dmat_3_frames <- simulate_SL_fun_colwise(data_mat = L_veg_SIF_frames, A_mat = A_Mat_0.4x)
veg_SL_Dmat_4_frames <- simulate_SL_fun_colwise(data_mat = L_veg_SIF_frames, A_mat = A_Mat_0.04x)
veg_SL_orig_frames <- simulate_SL_fun_colwise(data_mat = L_veg_SIF_frames, A_mat = A_Mat) # Original (non-scaled A matrix)


# Simulate stray light for solar spectrum.
# (Simulation of stray light and the addition to the 'measured' spectrum happen all in one step since this uses matrices)
# Each of the 500 frames
sol_SL_Dmat_1_frames <- simulate_SL_fun_colwise(data_mat = t(SolarData_SLC_frames), A_mat = A_Mat_40x)
sol_SL_Dmat_2_frames <- simulate_SL_fun_colwise(data_mat = t(SolarData_SLC_frames), A_mat = A_Mat_4x)
sol_SL_Dmat_3_frames <- simulate_SL_fun_colwise(data_mat = t(SolarData_SLC_frames), A_mat = A_Mat_0.4x)
sol_SL_Dmat_4_frames <- simulate_SL_fun_colwise(data_mat = t(SolarData_SLC_frames), A_mat = A_Mat_0.04x)
sol_SL_orig_frames <- simulate_SL_fun_colwise(data_mat = t(SolarData_SLC_frames), A_mat = A_Mat)


### Calculate mean of each cube of 500 frames
L_veg_SIF_frames_mean <- apply(L_veg_SIF_frames, 1, mean)
SolarData_SLC_frames_mean <- apply(t(SolarData_SLC_frames), 1, mean)
sol_SL_Dmat_1_frames_mean <- apply(sol_SL_Dmat_1_frames, 1, mean)
sol_SL_Dmat_2_frames_mean <- apply(sol_SL_Dmat_2_frames, 1, mean)
sol_SL_Dmat_3_frames_mean <- apply(sol_SL_Dmat_3_frames, 1, mean)
sol_SL_Dmat_4_frames_mean <- apply(sol_SL_Dmat_4_frames, 1, mean)
sol_SL_orig_frames_mean <- apply(sol_SL_orig_frames, 1, mean)
veg_SL_Dmat_1_frames_mean <- apply(veg_SL_Dmat_1_frames, 1, mean)
veg_SL_Dmat_2_frames_mean <- apply(veg_SL_Dmat_2_frames, 1, mean)
veg_SL_Dmat_3_frames_mean <- apply(veg_SL_Dmat_3_frames, 1, mean)
veg_SL_Dmat_4_frames_mean <- apply(veg_SL_Dmat_4_frames, 1, mean)
veg_SL_orig_frames_mean <- apply(veg_SL_orig_frames, 1, mean)

## Use mean of frames for the rest of the analysis (these designations mean that 
# the radiometric conversions, etc in R code above with mean spectra, from earlier version
# of make_stray_light_simulations was not necessary and thus removed)
# I added this so that I could use 500 frames for Frankenberg's retrieval
# and their mean for Alonso's retrieval (otherwise the 500 frames result
# and the mean result would come from different R/Matlab pipelines)
L_veg_SIF <- L_veg_SIF_frames_mean
SolarData_SLC <- SolarData_SLC_frames_mean
sol_SL_Dmat_1 <- sol_SL_Dmat_1_frames_mean
sol_SL_Dmat_2 <- sol_SL_Dmat_2_frames_mean
sol_SL_Dmat_3 <- sol_SL_Dmat_3_frames_mean
sol_SL_Dmat_4 <- sol_SL_Dmat_4_frames_mean
sol_SL_orig <- sol_SL_orig_frames_mean
veg_SL_Dmat_1_SIF <- veg_SL_Dmat_1_frames_mean
veg_SL_Dmat_2_SIF <- veg_SL_Dmat_2_frames_mean
veg_SL_Dmat_3_SIF <- veg_SL_Dmat_3_frames_mean
veg_SL_Dmat_4_SIF <- veg_SL_Dmat_4_frames_mean
veg_SL_orig_SIF <- veg_SL_orig_frames_mean

# Note: Simulation of veg_SL_Dmat_1, *_2, *_3, and *_4 all used L_veg_SIF (which already had SIF added) 
# and so the matrix multiplication gives simulated 'measured spectra.' 
# Both SIF and stray light are already components of the veg_SL_Dmat_1_SIF, veg_SL_Dmat_2_SIF, 
# veg_SL_Dmat_3_SIF and veg_SL_Dmat_4_SIF objects. Including 'SIF' in the object names makes this more clear.


# Even though simulation of stray light and its addition to the 'measured spectrum happen all in one step because of the matrix multiplication,
# it is useful to have the stray light quantities as separate objects (e.g. for the lm() later in this script), so define stray light quantities
# separately here:
# Single (mean) frame veg and solar:
veg_SL_magnitude_Dmat_1 <- veg_SL_Dmat_1_SIF - L_veg_SIF
veg_SL_magnitude_Dmat_2 <- veg_SL_Dmat_2_SIF - L_veg_SIF
veg_SL_magnitude_Dmat_3 <- veg_SL_Dmat_3_SIF - L_veg_SIF
veg_SL_magnitude_Dmat_4 <- veg_SL_Dmat_4_SIF - L_veg_SIF
veg_SL_magnitude_orig <- veg_SL_orig_SIF - L_veg_SIF
sol_SL_magnitude_Dmat_1 <- sol_SL_Dmat_1 - SolarData_SLC
sol_SL_magnitude_Dmat_2 <- sol_SL_Dmat_2 - SolarData_SLC
sol_SL_magnitude_Dmat_3 <- sol_SL_Dmat_3 - SolarData_SLC
sol_SL_magnitude_Dmat_4 <- sol_SL_Dmat_4 - SolarData_SLC
sol_SL_magnitude_orig <- sol_SL_orig - SolarData_SLC
# All 500 frames veg and solar:
# [Not needed, since SolarData_SLC now equals SolarData_SLC_frames; and L_veg_SIF now equals L_veg_SIF_frames]







##### PLOTS TO MAKE SURE THAT ALL SIMULATED SPECTRA ARE AS EXPECTED #########################################################################

# quick plot of solar spectrum and simulated vegetation spectrum based on L = E * R
plot(xVals, SolarData_SLC, type = "l")
lines(xVals, L_veg_frames_mean, col = "green")
lines(xVals, SIF_2perc * 50, col = "red")
plot(xVals, SIF_2perc, col = "red")



# Make sure A, I, and D matrices are as expected (based on Zong et al 2006):
# A is nearly the identity matrix: All diagonal components are 1, with adjacent components all 0 
# (􏰉di, J 􏰃= 0 inside the defined set of IB elements). The rest of the components in the matrix 
# are typically 3 orders smaller than the diagonal elements (the values of di, J are typically smaller than 10^-􏰁3). 
plot(raster(A_Mat))
plot(raster(log(A_Mat)))
plot(raster(A_Mat[1:10,1:10], xmn=1,xmx=10,ymn=1,ymx=10))
plot(raster(A_Mat[1:100,1:100], xmn=1,xmx=100,ymn=1,ymx=100))
plot(raster(log(A_Mat[1:100,1:100]), xmn=1,xmx=100,ymn=1,ymx=100))
plot(raster(D_Mat))
plot(raster(log(D_Mat)))
plot(D_Mat[ ,50], type = "l")
plot(D_Mat[ ,50], type = "l", log = "y", xlim = c(50-20,50+20))
plot(A_Mat[ ,50], type = "l", log = "y")
A_Mat[1,1] - I_mat[1,1] #I expected this to be zero, but instead it is a very small number. I attribute this to the ginv() function calculating A. Should be too small to matter, and it is a real empirical SLC matrix.
plot((Cmatrix[ ,50]+0.0001), log = "y", type = "l")

# Tested whether the small non-zero numbers in the IB region affected results by plugging in this 'test' D matrix. It didn't make much difference at all.
# D_Mat_test <- D_Mat
# for(i in 15:2145){
#   ind <- (i-15):(i+15)
#   D_Mat_test[ind,i] <- D_Mat_test[ind,i] / 100
# }
# plot(D_Mat_test[ ,50], type = "l", log = "y", xlim = c(50-20,50+20))
# lines(D_Mat[ ,50], col = "blue", log = "y", xlim = c(50-20,50+20))
# plot(D_Mat_test[ ,2050], type = "l", log = "y", xlim = c(2050-30,2050+30))
# lines(D_Mat[ ,2050], col = "blue", log = "y", xlim = c(2050-30,2050+30))
#D_Mat <- D_Mat_test

# quick plots of stray light simulated by D matrix
plot(veg_SL_Dmat_1_SIF - L_veg_SIF, type = "l")

# quick log scale plots of stray light simulated by D matrix, as a proportion of SolarData_SLC & L_veg_SIF, where D was re-scaled to simulate different quantities of stray light
plot((sol_SL_Dmat_1 - SolarData_SLC)/SolarData_SLC, type = "l", log = "y") #Stray light is 0-1 orders of magnitude of the signal
lines((veg_SL_Dmat_1_SIF - L_veg_SIF)/L_veg_SIF, col = "blue") 
abline(h = c(0.1, 1))
plot((sol_SL_Dmat_2 - SolarData_SLC)/SolarData_SLC, type = "l", log = "y") #Stray light is 2-3 orders of magnitude of the signal
lines((veg_SL_Dmat_2_SIF - L_veg_SIF)/L_veg_SIF, col = "blue")  
abline(h = c(0.1,0.01))
plot((sol_SL_Dmat_3 - SolarData_SLC)/SolarData_SLC, type = "l", log = "y") #Stray light is 2-3 orders of magnitude of the signal
lines((veg_SL_Dmat_3_SIF - L_veg_SIF)/L_veg_SIF, col = "blue")  
abline(h = c(0.001,0.01))
plot((sol_SL_Dmat_4 - SolarData_SLC)/SolarData_SLC, type = "l", log = "y") 
lines((veg_SL_Dmat_4_SIF - L_veg_SIF)/L_veg_SIF, col = "blue")   #Stray light is 4-5 orders of magnitude of the signal
abline(h = c(0.00001, 0.0001)) 



# quick plot of a veg_SL_Dmat object to see if subtracting stray light and SIF reproduces L_veg
test_sim <- veg_SL_Dmat_1_SIF - (veg_SL_Dmat_1_SIF - L_veg_SIF) - SIF_2perc
plot(Wavelengths$spec_pixel_list, test_sim, type = "l")
lines(Wavelengths$spec_pixel_list, L_veg_frames_mean, col = "green")
all.equal(as.numeric(test_sim), as.numeric(L_veg_frames_mean)) # is equal, good.

# quick plots of stray light simulated by D matrix, as a proportion of SolarData_SLC, where D was re-scaled to simulate different quantities of stray light
plot((sol_SL_Dmat_1 - SolarData_SLC)/SolarData_SLC, type = "l")
plot((sol_SL_Dmat_2 - SolarData_SLC)/SolarData_SLC, type = "l")
plot((sol_SL_Dmat_3 - SolarData_SLC)/SolarData_SLC, type = "l")
plot((sol_SL_Dmat_4 - SolarData_SLC)/SolarData_SLC, type = "l")




##### DEFINE FRAUNHOFER/OXYGEN BANDS AND REFERENCE BANDS TO BE USED FOR RETRIEVALS ######################################################

# Define windows around Fraunhofer lines (same as used in CompareWavelengthCalibration.R for FWHM estimation)
CR_breaks <- data.frame(Start=rev(c(777.9, 774, 772.65, 756.6, 755.4, 750.90, 746, 744.5, 682.65, 680.80, 672.40, 671.5)),
                        End=rev(c(778.8, 774.7, 773.45, 757.35, 756.05, 751.7, 747, 745.35, 683.35, 681.5, 673.3, 672.25)))

# Find bottom of the Fraunhofer line within each window (on non-veg target)
SolarData_SLC_df <- cbind(Wavelengths, SolarData_SLC)
window_mins <- rep(NA, times = length(CR_breaks))
window_mins_ind <- rep(NA, times = length(CR_breaks))
for(i in 1:dim(CR_breaks)[1]){
  temp_ind <- SolarData_SLC_df[ ,2] >= CR_breaks$Start[i] & SolarData_SLC_df[ ,2] <= CR_breaks$End[i]
  temp <-SolarData_SLC_df[temp_ind, , ]
  window_mins[i]<- temp$spec_pixel_pred_no710[which(temp[ ,3] == min(temp[ ,3]))]
  window_mins_ind[i] <- which(Wavelengths$spec_pixel_pred_no710==window_mins[i])}


# Now define wavelength interval for finding the maximum (the left and right shoulders, i.e. the reference spectra)
interval <- 1 # (nm)
WI_min <- rep(NA, times = length(CR_breaks))
WI_max <- rep(NA, times = length(CR_breaks))
for(i in 1:length(window_mins)){
  temp1 <- window_mins[i] - interval/2
  WI_min[i] <- Wavelengths$spec_pixel_pred_no710[which(abs(Wavelengths$spec_pixel_pred_no710 - temp1)==min(abs(Wavelengths$spec_pixel_pred_no710 - temp1)))]
  temp2 <- window_mins[i] + interval/2
  WI_max[i] <- Wavelengths$spec_pixel_pred_no710[which(abs(Wavelengths$spec_pixel_pred_no710 - temp2)==min(abs(Wavelengths$spec_pixel_pred_no710 - temp2)))]
}

# Find fraunhofer line shoulders, and plot shoulders, wavelength interval, and mins to make sure code worked
left_shoulder <- matrix(NA, nrow=length(window_mins), ncol=3)
right_shoulder <- matrix(NA, nrow=length(window_mins), ncol=3)
plot(SolarData_SLC_df[,2], SolarData_SLC_df[,3], type = "l")
abline(v = c(window_mins), col = "red")
for(i in 1:length(window_mins)){
  
  # Plot polygon for wavelength interval for finding the shoulders
  polygon(x=c(WI_min[i],WI_max[i],WI_max[i],WI_min[i]),y=c(0,0,max(SolarData_SLC_df),max(SolarData_SLC_df)),
          col=adjustcolor("grey",alpha.f=0.25), border=NA)
  
  # Find and plot the left shoulder (blue)
  left_interval <- which(Wavelengths$spec_pixel_pred_no710 == WI_min[i]):which(Wavelengths$spec_pixel_pred_no710 == window_mins[i])
  left_shoulder_ind <- which(SolarData_SLC_df[left_interval, 3] == max(SolarData_SLC_df[left_interval, 3]))
  if(length(left_shoulder_ind) > 1){
    closest_left_shoulder_ind <- which(abs(left_interval[left_shoulder_ind] - window_mins[i]) == min(abs(left_interval[left_shoulder_ind] - window_mins[i])))
    left_shoulder_ind <- left_shoulder_ind[closest_left_shoulder_ind]
  }
  left_shoulder[i, 1] <- SolarData_SLC_df[left_interval[left_shoulder_ind], 1]
  left_shoulder[i, 2] <- SolarData_SLC_df[left_interval[left_shoulder_ind], 2]
  left_shoulder[i, 3] <- SolarData_SLC_df[left_interval[left_shoulder_ind], 3]
  abline(v = left_shoulder[i,2], col = "blue")
  
  # Find and plot the right shoulder (cyan)
  right_interval <- which(Wavelengths$spec_pixel_pred_no710 == window_mins[i]):which(Wavelengths$spec_pixel_pred_no710 == WI_max[i])
  right_shoulder_ind <- which(SolarData_SLC_df[right_interval, 3] == max(SolarData_SLC_df[right_interval, 3]))
  if(length(right_shoulder_ind) > 1){
    closest_right_shoulder_ind <- which(abs(right_interval[right_shoulder_ind] - window_mins[i]) == min(abs(right_interval[right_shoulder_ind] - window_mins[i])))
    right_shoulder_ind <- right_shoulder_ind[closest_right_shoulder_ind]
  }
  right_shoulder[i, 1] <- SolarData_SLC_df[right_interval[right_shoulder_ind], 1]
  right_shoulder[i, 2] <- SolarData_SLC_df[right_interval[right_shoulder_ind], 2]
  right_shoulder[i, 3] <- SolarData_SLC_df[right_interval[right_shoulder_ind], 3]
  abline(v = right_shoulder[i,2], col = "cyan")
  
}


# # Define reference bands for oxygen bands and then add to matrix of references for Fraunhofer lines
# # Oxygen B band (based on ref: Ni et al 2018 Sensors sensitivity analysis section, and visual examination of spectra for O2B_right)
# O2B_left <- 685
# O2B_min <- 687
# O2B_right <- ? No good candidate 
# Oxygen A band (based on ref: Damm et al 2014 RSE equation 7)
O2A_left_coarse <- 758
O2A_min_coarse <- 760
O2A_right_coarse <- 771
# Using the mean solar spectra, find the exact wavelengths for the band and reference spectra (because the resolution is so high that there are multiple options)
O2A_left_candidates <- which(round(Wavelengths$spec_pixel_pred_no710, digits = 0) == O2A_left_coarse)
O2A_left <- Wavelengths$spec_pixel_pred_no710[which(SolarData_SLC == max(SolarData_SLC[O2A_left_candidates]))]
O2A_min_candidates <- which(round(Wavelengths$spec_pixel_pred_no710, digits = 0) == O2A_min_coarse)
O2A_min <- Wavelengths$spec_pixel_pred_no710[which(SolarData_SLC == min(SolarData_SLC[O2A_min_candidates]))]
O2A_right_candidates <- which(round(Wavelengths$spec_pixel_pred_no710, digits = 0) == O2A_right_coarse)
O2A_right <- Wavelengths$spec_pixel_pred_no710[which(SolarData_SLC == max(SolarData_SLC[O2A_right_candidates]))]


# Plot polygon for O2A wavelength interval 
plot(Wavelengths$spec_pixel_pred_no710[1569:2160], SolarData_SLC[1569:2160], type = "l")
polygon(x=c(O2A_left,O2A_right,O2A_right,O2A_left),y=c(0,0,max(SolarData_SLC),max(SolarData_SLC)),
        col=adjustcolor("grey",alpha.f=0.25), border=NA)
# Plot the left shoulder (blue), min (black), and right shoulder (cyan)
abline(v = O2A_left, col = "blue")
abline(v = O2A_min, col = "black")
abline(v = O2A_right, col = "cyan")


## Add O2A bands to matrices with Fraunhofer line minima and shoulders
# Minima
window_mins <- c(window_mins, O2A_min)
O2A_min_ind <- Wavelengths$spec_pixel_list[which(SolarData_SLC == min(SolarData_SLC[O2A_min_candidates]))]
window_mins_ind <- c(window_mins_ind, O2A_min_ind)
# Left shoulder
O2A_left_ind <- Wavelengths$spec_pixel_list[which(SolarData_SLC == max(SolarData_SLC[O2A_left_candidates]))]
O2A_left_wl <- Wavelengths$spec_pixel_pred_no710[O2A_left_ind]
O2A_left_shoulder <- data.matrix(SolarData_SLC[O2A_left_ind])
colnames(O2A_left_shoulder) <- NULL
rownames(O2A_left_shoulder) <- NULL
left_shoulder <- rbind(left_shoulder, c(O2A_left_ind, O2A_left_wl, O2A_left_shoulder))
# right shoulder
O2A_right_ind <- Wavelengths$spec_pixel_list[which(SolarData_SLC == max(SolarData_SLC[O2A_right_candidates]))]
O2A_right_wl <- Wavelengths$spec_pixel_pred_no710[O2A_right_ind]
O2A_right_shoulder <- data.matrix(SolarData_SLC[O2A_right_ind])
colnames(O2A_right_shoulder) <- NULL
rownames(O2A_right_shoulder) <- NULL
right_shoulder <- rbind(right_shoulder, c(O2A_right_ind, O2A_right_wl, O2A_right_shoulder))

# Make plots zoomed in on red and far-red regions separately
plot(Wavelengths$spec_pixel_pred_no710[1:400], SolarData_SLC[1:400], type = "l")
abline(v = c(window_mins), col = "red")
abline(v = c(left_shoulder[,2]), col = "blue")
abline(v = c(right_shoulder[,2]), col = "cyan")
for(i in 1:length(window_mins)){
  # Plot polygon for wavelength interval for finding the shoulders
  polygon(x=c(WI_min[i],WI_max[i],WI_max[i],WI_min[i]),y=c(0,0,max(SolarData_SLC),max(SolarData_SLC)),
          col=adjustcolor("grey",alpha.f=0.25), border=NA)
}
plot(Wavelengths$spec_pixel_pred_no710[1380:2160], SolarData_SLC[1380:2160], type = "l")
abline(v = c(window_mins), col = "red")
abline(v = c(left_shoulder[,2]), col = "blue")
abline(v = c(right_shoulder[,2]), col = "cyan")
for(i in 1:length(window_mins)){
  # Plot polygon for wavelength interval for finding the shoulders
  polygon(x=c(WI_min[i],WI_max[i],WI_max[i],WI_min[i]),y=c(0,0,max(SolarData_SLC),max(SolarData_SLC)),
          col=adjustcolor("grey",alpha.f=0.25), border=NA)
}

# Combine FLD line minima wavelengths and wavelength indices
window_mins_mat <- cbind(window_mins_ind, window_mins)



##### 3FLD SIF RETRIEVALS FOR THE VARIOUS SPECTRA WITH AND WITHOUT STRAY LIGHT COMPONENTS ####################################################


### NO SL
### 3FLD retrieval for simulated spectra with no stray light
# loop through all retrieval windows

three_FLD_L_veg_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)

for(k in 1:dim(window_mins_mat)[1]){
  
  three_FLD_L_veg_SIF[k] <-
    threeFLD(
      lambdaL = left_shoulder[k, 2],
      lambdaR = right_shoulder[k, 2],
      lambdaIN = window_mins[k],
      LoutL = L_veg_SIF[left_shoulder[k, 1]],
      LoutR = L_veg_SIF[right_shoulder[k, 1]],
      Lin = L_veg_SIF[window_mins_mat[k, 1]],
      EoutL = SolarData_SLC[left_shoulder[k, 1]],
      EoutR = SolarData_SLC[right_shoulder[k, 1]],
      Ein = SolarData_SLC[window_mins_mat[k, 1]]
    )
}

# Quick plot to see if 3FLD retrieval makes sense
plot(window_mins_mat[,2], three_FLD_L_veg_SIF, col = "blue")
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")





# D_MATRIX stray light scenarios
### 3FLD retrieval for simulated spectra with 'D matrix simulated' stray light
# loop through all retrieval windows

three_FLD_veg_SL_Dmat_1_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)
three_FLD_veg_SL_Dmat_2_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)
three_FLD_veg_SL_Dmat_3_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)
three_FLD_veg_SL_Dmat_4_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)

for(k in 1:dim(window_mins_mat)[1]){
  
 # 3FLD for spectra simulated with 'D matrix' and roughly 1 order of magnitude lower than signal
  three_FLD_veg_SL_Dmat_1_SIF[k] <-
    threeFLD(
      lambdaL = left_shoulder[k, 2],
      lambdaR = right_shoulder[k, 2],
      lambdaIN = window_mins[k],
      LoutL = veg_SL_Dmat_1_SIF[left_shoulder[k, 1]],
      LoutR = veg_SL_Dmat_1_SIF[right_shoulder[k, 1]],
      Lin = veg_SL_Dmat_1_SIF[window_mins_mat[k, 1]],
      EoutL = sol_SL_Dmat_1[left_shoulder[k, 1]],
      EoutR = sol_SL_Dmat_1[right_shoulder[k, 1]],
      Ein = sol_SL_Dmat_1[window_mins_mat[k, 1]]
    )
  # 3FLD for spectra simulated with 'D matrix' and roughly 2 orders of magnitude lower than signal
  three_FLD_veg_SL_Dmat_2_SIF[k] <-
    threeFLD(
      lambdaL = left_shoulder[k, 2],
      lambdaR = right_shoulder[k, 2],
      lambdaIN = window_mins[k],
      LoutL = veg_SL_Dmat_2_SIF[left_shoulder[k, 1]],
      LoutR = veg_SL_Dmat_2_SIF[right_shoulder[k, 1]],
      Lin = veg_SL_Dmat_2_SIF[window_mins_mat[k, 1]],
      EoutL = sol_SL_Dmat_2[left_shoulder[k, 1]],
      EoutR = sol_SL_Dmat_2[right_shoulder[k, 1]],
      Ein = sol_SL_Dmat_2[window_mins_mat[k, 1]]
    )
  # 3FLD for spectra simulated with 'D matrix' and roughly 3 orders of magnitude lower than signal
  three_FLD_veg_SL_Dmat_3_SIF[k] <-
    threeFLD(
      lambdaL = left_shoulder[k, 2],
      lambdaR = right_shoulder[k, 2],
      lambdaIN = window_mins[k],
      LoutL = veg_SL_Dmat_3_SIF[left_shoulder[k, 1]],
      LoutR = veg_SL_Dmat_3_SIF[right_shoulder[k, 1]],
      Lin = veg_SL_Dmat_3_SIF[window_mins_mat[k, 1]],
      EoutL = sol_SL_Dmat_3[left_shoulder[k, 1]],
      EoutR = sol_SL_Dmat_3[right_shoulder[k, 1]],
      Ein = sol_SL_Dmat_3[window_mins_mat[k, 1]]
    )
  # 3FLD for spectra wsimulated with 'D matrix' and roughly 4 orders of magnitude lower than signal  
  three_FLD_veg_SL_Dmat_4_SIF[k] <-
    threeFLD(
      lambdaL = left_shoulder[k, 2],
      lambdaR = right_shoulder[k, 2],
      lambdaIN = window_mins[k],
      LoutL = veg_SL_Dmat_4_SIF[left_shoulder[k, 1]],
      LoutR = veg_SL_Dmat_4_SIF[right_shoulder[k, 1]],
      Lin = veg_SL_Dmat_4_SIF[window_mins_mat[k, 1]],
      EoutL = sol_SL_Dmat_4[left_shoulder[k, 1]],
      EoutR = sol_SL_Dmat_4[right_shoulder[k, 1]],
      Ein = sol_SL_Dmat_4[window_mins_mat[k, 1]]
    )
}



### 3FLD retrieval for different sensors scenario: simulated spectra with 'D matrix simulated' 
### stray light, but veg and ref have different stray light levels
# loop through all retrieval windows

three_FLD_veg_SL_Dmat_2to3_SIF_mismatch <- matrix(NA, nrow = length(window_mins), ncol = 1)

for(k in 1:dim(window_mins_mat)[1]){
  
  # 3FLD for spectra simulated with mismatching 'D matrices': roughly 2 order of magnitude 
  # lower than signal for reference and 3 order or magnitude for Veg
  three_FLD_veg_SL_Dmat_2to3_SIF_mismatch[k] <-
    threeFLD(
      lambdaL = left_shoulder[k, 2],
      lambdaR = right_shoulder[k, 2],
      lambdaIN = window_mins[k],
      LoutL = veg_SL_Dmat_3_SIF[left_shoulder[k, 1]],
      LoutR = veg_SL_Dmat_3_SIF[right_shoulder[k, 1]],
      Lin = veg_SL_Dmat_3_SIF[window_mins_mat[k, 1]],
      EoutL = sol_SL_Dmat_2[left_shoulder[k, 1]],
      EoutR = sol_SL_Dmat_2[right_shoulder[k, 1]],
      Ein = sol_SL_Dmat_2[window_mins_mat[k, 1]]
    )
}

three_FLD_veg_SL_Dmat_3to2_SIF_mismatch <- matrix(NA, nrow = length(window_mins), ncol = 1)

for(k in 1:dim(window_mins_mat)[1]){
  
  # 3FLD for spectra simulated with mismatching 'D matrices': roughly 3 order of magnitude 
  # lower than signal for reference and 2 order of magnitude for Veg
  three_FLD_veg_SL_Dmat_3to2_SIF_mismatch[k] <-
    threeFLD(
      lambdaL = left_shoulder[k, 2],
      lambdaR = right_shoulder[k, 2],
      lambdaIN = window_mins[k],
      LoutL = veg_SL_Dmat_2_SIF[left_shoulder[k, 1]],
      LoutR = veg_SL_Dmat_2_SIF[right_shoulder[k, 1]],
      Lin = veg_SL_Dmat_2_SIF[window_mins_mat[k, 1]],
      EoutL = sol_SL_Dmat_3[left_shoulder[k, 1]],
      EoutR = sol_SL_Dmat_3[right_shoulder[k, 1]],
      Ein = sol_SL_Dmat_3[window_mins_mat[k, 1]]
    )
}






##### sFLD SIF RETRIEVALS FOR THE VARIOUS SPECTRA WITH AND WITHOUT STRAY LIGHT COMPONENTS ####################################################


### NO SL
### sFLD retrieval for simulated spectra with no stray light
# loop through all retrieval windows
sFLD_L_veg_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)

for(k in 1:dim(window_mins_mat)[1]){
  
  sFLD_L_veg_SIF[k] <- sFLD(Lout = L_veg_SIF[left_shoulder[k,1]], Lin = L_veg_SIF[window_mins_mat[k, 1]],
                            Eout = SolarData_SLC[left_shoulder[k,1]], Ein = SolarData_SLC[window_mins_mat[k, 1]])
  
}

# Quick plot to see if sFLD retrieval makes sense
plot(window_mins_mat[,2], sFLD_L_veg_SIF, col = "black")
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
# Not as good as 3FlD




### D_MATRIX stray light scenario
### sFLD retrieval for simulated spectra with 'D matrix simulated' stray light
# loop through all retrieval windows

sFLD_veg_SL_Dmat_1_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)
sFLD_veg_SL_Dmat_2_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)
sFLD_veg_SL_Dmat_3_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)
sFLD_veg_SL_Dmat_4_SIF <- matrix(NA, nrow = length(window_mins), ncol = 1)

for(k in 1:dim(window_mins_mat)[1]) {
  sFLD_veg_SL_Dmat_1_SIF[k] <-
    sFLD(
      Lout = veg_SL_Dmat_1_SIF[left_shoulder[k, 1]],
      Lin = veg_SL_Dmat_1_SIF[window_mins_mat[k, 1]],
      Eout = sol_SL_Dmat_1[left_shoulder[k, 1]],
      Ein = sol_SL_Dmat_1[window_mins_mat[k, 1]]
    )
  sFLD_veg_SL_Dmat_2_SIF[k] <-
    sFLD(
      Lout = veg_SL_Dmat_2_SIF[left_shoulder[k, 1]],
      Lin = veg_SL_Dmat_2_SIF[window_mins_mat[k, 1]],
      Eout = sol_SL_Dmat_2[left_shoulder[k, 1]],
      Ein = sol_SL_Dmat_2[window_mins_mat[k, 1]]
    )
  sFLD_veg_SL_Dmat_3_SIF[k] <-
    sFLD(
      Lout = veg_SL_Dmat_3_SIF[left_shoulder[k, 1]],
      Lin = veg_SL_Dmat_3_SIF[window_mins_mat[k, 1]],
      Eout = sol_SL_Dmat_3[left_shoulder[k, 1]],
      Ein = sol_SL_Dmat_3[window_mins_mat[k, 1]]
    )
  sFLD_veg_SL_Dmat_4_SIF[k] <-
    sFLD(
      Lout = veg_SL_Dmat_4_SIF[left_shoulder[k, 1]],
      Lin = veg_SL_Dmat_4_SIF[window_mins_mat[k, 1]],
      Eout = sol_SL_Dmat_4[left_shoulder[k, 1]],
      Ein = sol_SL_Dmat_4[window_mins_mat[k, 1]]
    )
}




#### PLOTS TO EXAMINE HOW MUCH STRAY LIGHT MATTERS ####################################################################################



### D_MATRIX SL
# Quick plots to see how 'D matrix simulated' stray light matters with 3FLD
ymax <- max(three_FLD_veg_SL_Dmat_1_SIF) + 0.01
plot(window_mins_mat[,2], three_FLD_L_veg_SIF, col = "black", ylim = c(0, ymax))
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
points(window_mins_mat[,2], three_FLD_veg_SL_Dmat_1_SIF, col = "#08519c")

plot(window_mins_mat[,2], three_FLD_L_veg_SIF, col = "black", ylim = c(0, ymax))
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
points(window_mins_mat[,2], three_FLD_veg_SL_Dmat_2_SIF, col = "#08519c")

plot(window_mins_mat[,2], three_FLD_L_veg_SIF, col = "black", ylim = c(0, ymax))
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
points(window_mins_mat[,2], three_FLD_veg_SL_Dmat_3_SIF, col = "#08519c")

plot(window_mins_mat[,2], three_FLD_L_veg_SIF, col = "black", ylim = c(0, ymax))
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
points(window_mins_mat[,2], three_FLD_veg_SL_Dmat_4_SIF, col = "#08519c")


# Quick plots to see how 'D matrix simulated' stray light matters with standard FLD
ymax <- max(sFLD_veg_SL_Dmat_1_SIF) + 0.01
plot(window_mins_mat[,2], sFLD_L_veg_SIF, col = "black", ylim = c(0, ymax))
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
points(window_mins_mat[,2], sFLD_veg_SL_Dmat_1_SIF, col = "#08519c")

plot(window_mins_mat[,2], sFLD_L_veg_SIF, col = "black", ylim = c(0, ymax))
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
points(window_mins_mat[,2], sFLD_veg_SL_Dmat_2_SIF, col = "#08519c")

plot(window_mins_mat[,2], sFLD_L_veg_SIF, col = "black", ylim = c(0, ymax))
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
points(window_mins_mat[,2], sFLD_veg_SL_Dmat_3_SIF, col = "#08519c")

plot(window_mins_mat[,2], sFLD_L_veg_SIF, col = "black", ylim = c(0, ymax))
lines(Wavelengths$spec_pixel_pred_no710, SIF_2perc, col = "red")
points(window_mins_mat[,2], sFLD_veg_SL_Dmat_4_SIF, col = "#08519c")



#### CALCULATE DIFFERENCE DUE TO STRAY LIGHT AS A PROPORTION OF TRUE SIF ################################################



### D_MATRIX SL
# 3FLD
three_FLD_veg_SL_Dmat_1_SIF_error <- (three_FLD_veg_SL_Dmat_1_SIF - three_FLD_L_veg_SIF)/SIF_2perc[window_mins_ind]
three_FLD_veg_SL_Dmat_2_SIF_error <- (three_FLD_veg_SL_Dmat_2_SIF - three_FLD_L_veg_SIF)/SIF_2perc[window_mins_ind]
three_FLD_veg_SL_Dmat_3_SIF_error <- (three_FLD_veg_SL_Dmat_3_SIF - three_FLD_L_veg_SIF)/SIF_2perc[window_mins_ind]
three_FLD_veg_SL_Dmat_4_SIF_error <- (three_FLD_veg_SL_Dmat_4_SIF - three_FLD_L_veg_SIF)/SIF_2perc[window_mins_ind]
three_FLD_veg_SL_Dmat_2to3_SIF_mismatch_error <- (three_FLD_veg_SL_Dmat_2to3_SIF_mismatch - three_FLD_L_veg_SIF)/SIF_2perc[window_mins_ind]
three_FLD_veg_SL_Dmat_3to2_SIF_mismatch_error <- (three_FLD_veg_SL_Dmat_3to2_SIF_mismatch - three_FLD_L_veg_SIF)/SIF_2perc[window_mins_ind]

# standard FLD
sFLD_veg_SL_Dmat_1_SIF_error <- (sFLD_veg_SL_Dmat_1_SIF - sFLD_L_veg_SIF)/SIF_2perc[window_mins_ind]
sFLD_veg_SL_Dmat_2_SIF_error <- (sFLD_veg_SL_Dmat_2_SIF - sFLD_L_veg_SIF)/SIF_2perc[window_mins_ind]
sFLD_veg_SL_Dmat_3_SIF_error <- (sFLD_veg_SL_Dmat_3_SIF - sFLD_L_veg_SIF)/SIF_2perc[window_mins_ind]
sFLD_veg_SL_Dmat_4_SIF_error <- (sFLD_veg_SL_Dmat_4_SIF - sFLD_L_veg_SIF)/SIF_2perc[window_mins_ind]



#### TABLE SUMMARIZING HOW MUCH STRAY LIGHT MATTERS ##########################################################################


### Make table for microsoft word
# SNR_wavelengths_ind <- c(200, 330, 1686, 1706, 1754)
# SNR_wavelengths <- SNR[SNR_wavelengths_ind, ]


stray_summary_table_3FLD_Dmat <- t(cbind(three_FLD_veg_SL_Dmat_1_SIF_error,
                                         three_FLD_veg_SL_Dmat_2_SIF_error,
                                         three_FLD_veg_SL_Dmat_3_SIF_error,
                                         three_FLD_veg_SL_Dmat_4_SIF_error))

stray_summary_table_sFLD_Dmat <- t(cbind(sFLD_veg_SL_Dmat_1_SIF_error,
                                         sFLD_veg_SL_Dmat_2_SIF_error,
                                         sFLD_veg_SL_Dmat_3_SIF_error,
                                         sFLD_veg_SL_Dmat_4_SIF_error))

# Summary table of simulations to be reported
stray_summary_table <- rbind(
  stray_summary_table_sFLD_Dmat,
  stray_summary_table_3FLD_Dmat
)

# Convert from proportion to percentage for the table
stray_summary_table_perc <- stray_summary_table * 100

htmlTable(round(stray_summary_table_perc, digits = 2), header =  as.character(round(window_mins, digits = 2)), 
          rnames = paste("Scenario ", as.character(rep(c(1,2,3,4), times = 2)), "orders of magnitude"))
# The html table can be pasted directly into word, or, as a hack, the html table can be pasted into excel then word and can be edited in word.

## Some useful text for the results
# Fraunhofer lines
three_FLD_veg_SL_Dmat_1_SIF_error_perc_mean_Fraunhofer <- mean(three_FLD_veg_SL_Dmat_1_SIF_error[1:12]) * 100
three_FLD_veg_SL_Dmat_2_SIF_error_perc_mean_Fraunhofer <- mean(three_FLD_veg_SL_Dmat_2_SIF_error[1:12]) * 100
three_FLD_veg_SL_Dmat_3_SIF_error_perc_mean_Fraunhofer <- mean(three_FLD_veg_SL_Dmat_3_SIF_error[1:12]) * 100
three_FLD_veg_SL_Dmat_4_SIF_error_perc_mean_Fraunhofer <- mean(three_FLD_veg_SL_Dmat_4_SIF_error[1:12]) * 100
sFLD_veg_SL_Dmat_1_SIF_error_perc_mean_Fraunhofer <- mean(sFLD_veg_SL_Dmat_1_SIF_error[1:12]) * 100
sFLD_veg_SL_Dmat_2_SIF_error_perc_mean_Fraunhofer <- mean(sFLD_veg_SL_Dmat_2_SIF_error[1:12]) * 100
sFLD_veg_SL_Dmat_3_SIF_error_perc_mean_Fraunhofer <- mean(sFLD_veg_SL_Dmat_3_SIF_error[1:12]) * 100
sFLD_veg_SL_Dmat_4_SIF_error_perc_mean_Fraunhofer <- mean(sFLD_veg_SL_Dmat_4_SIF_error[1:12]) * 100
print_mean_perc_error <- function(x){
  print(paste("The mean stray light error as a percentage of SIF was", round(x, digits = 3), "% for", deparse(substitute(x))))
}
print_mean_perc_error(three_FLD_veg_SL_Dmat_1_SIF_error_perc_mean_Fraunhofer)
print_mean_perc_error(three_FLD_veg_SL_Dmat_2_SIF_error_perc_mean_Fraunhofer)
print_mean_perc_error(three_FLD_veg_SL_Dmat_3_SIF_error_perc_mean_Fraunhofer)
print_mean_perc_error(three_FLD_veg_SL_Dmat_4_SIF_error_perc_mean_Fraunhofer)
print_mean_perc_error(sFLD_veg_SL_Dmat_1_SIF_error_perc_mean_Fraunhofer)
print_mean_perc_error(sFLD_veg_SL_Dmat_2_SIF_error_perc_mean_Fraunhofer)
print_mean_perc_error(sFLD_veg_SL_Dmat_3_SIF_error_perc_mean_Fraunhofer)
print_mean_perc_error(sFLD_veg_SL_Dmat_4_SIF_error_perc_mean_Fraunhofer)
# O2_A_band
print_mean_perc_error_O2A <- function(x){
  O2A_perc <- x[13] * 100
  print(paste("The mean stray light error as a percentage of SIF was", round(O2A_perc, digits = 3), "% for", deparse(substitute(x))))
}
print_mean_perc_error_O2A(three_FLD_veg_SL_Dmat_1_SIF_error)
print_mean_perc_error_O2A(three_FLD_veg_SL_Dmat_2_SIF_error)
print_mean_perc_error_O2A(three_FLD_veg_SL_Dmat_3_SIF_error)
print_mean_perc_error_O2A(three_FLD_veg_SL_Dmat_4_SIF_error)
print_mean_perc_error_O2A(sFLD_veg_SL_Dmat_1_SIF_error)
print_mean_perc_error_O2A(sFLD_veg_SL_Dmat_2_SIF_error)
print_mean_perc_error_O2A(sFLD_veg_SL_Dmat_3_SIF_error)
print_mean_perc_error_O2A(sFLD_veg_SL_Dmat_4_SIF_error)
# Context for results:
print(paste('Rascher et al 2015 found that forest had', ((0.868369352 - 0.502946955)/0.502946955)*100,'% higher SIF @760 nm than grass'))
print(paste('Rascher et al 2015 found that corn had', ((0.919449902 - 0.868369352)/0.868369352)*100,'% higher SIF @760 nm than forest'))

# Note: Lee et al 2013 Proceedings of the royal academy B says: "However, over central Amazonia with a relatively moderate dry season of approximately three months
# (region B), fluorescence decreases significantly (approx. 27%) during the normal (2009) dry season compared with the wet season, while LAI increases by 11 per cent during the dry
# season. MODIS EVI, on the other hand, shows similar variations to fluorescence (r2 ¼ 0.52). Table 1 summarizes r2 between different parameters."


#### CODE QUALITY CONTROL #####################################################################################################

# Does everything make sense going into the retrievals?

# # Quick check of vegetation spectrum with and without stray light
# plot(L_veg, type = "l")
# test2 <- veg_SL_const_0 - veg_SL_magnitude_const_0
# lines(test2, col = "red")

# # Quick check of solar spectrum with and without stray light
# plot(SolarData_SLC, type = "l", ylim = c(0, 0.8))
# temp <- sol_SL_const_0 - sol_SL_magnitude_const_0
# lines(temp, col = "red")
# temp0 <- SolarData_SLC - temp
# all.equal(SolarData_SLC, temp)
# # Solar spectrum looks good

# There was a test of the derivation in earlier version (make_stray_light_simulations.R) that showed the 
# FLD derivation in the paper works when its assumptions are met.


test_L_veg_SIF <- L_veg_frames_mean + SIF_2perc

## There is one outlier with a very large negative error for Fraunhofer line at 778.266 nm. Double check a calculation of SIF error for this point.
test_veg_SL_Dmat_1 <- simulate_SL_fun_colwise(data_mat = as.matrix(test_L_veg_SIF), A_mat = A_Mat_40x)
test_sol_SL_Dmat_1 <- simulate_SL_fun_colwise(data_mat = as.matrix(SolarData_SLC), A_mat = A_Mat_40x)
test_veg_SL_Dmat_3 <- simulate_SL_fun_colwise(data_mat = as.matrix(test_L_veg_SIF), A_mat = A_Mat_0.4x)
test_sol_SL_Dmat_3 <- simulate_SL_fun_colwise(data_mat = as.matrix(SolarData_SLC), A_mat = A_Mat_0.4x)
# no stray sFLD for Fraunhofer line at 778.266
test_778_SIF <- sFLD(Lout = test_L_veg_SIF[left_shoulder[12,1]], Lin = test_L_veg_SIF[window_mins_mat[12, 1]],
                                             Eout = SolarData_SLC[left_shoulder[12,1]], Ein = SolarData_SLC[window_mins_mat[12, 1]])
# with stray sFLD for Fraunhofer line at 778.266
test_778_Dmat_1_SIF <- sFLD(Lout = test_veg_SL_Dmat_1[left_shoulder[12,1]], Lin = test_veg_SL_Dmat_1[window_mins_mat[12,1]],
                            Eout = test_sol_SL_Dmat_1[left_shoulder[12,1]], Ein = test_sol_SL_Dmat_1[window_mins_mat[12,1]])
test_778_Dmat_3_SIF <- sFLD(Lout = test_veg_SL_Dmat_3[left_shoulder[12,1]], Lin = test_veg_SL_Dmat_3[window_mins_mat[12,1]],
                            Eout = test_sol_SL_Dmat_3[left_shoulder[12,1]], Ein = test_sol_SL_Dmat_3[window_mins_mat[12,1]])
# Calculate error
test_778_err_1 <- (test_778_Dmat_1_SIF - test_778_SIF)/SIF_2perc[window_mins_ind[12]]
test_778_err_3 <- (test_778_Dmat_3_SIF - test_778_SIF)/SIF_2perc[window_mins_ind[12]]
# Error matches result from main script in sFLD_veg_SL_Dmat_1_df and sFLD_veg_SL_Dmat_3_df. Examine plot:
plot_ind_778 <- (window_mins_mat[12,1] - 100) : (window_mins_mat[12,1] + 39)
plot(test_L_veg_SIF[plot_ind_778], type = "l", ylim = c(0.11, 0.21))
lines(test_veg_SL_Dmat_1[plot_ind_778], type = "l", col = "blue")
abline(v = 101)
plot(SolarData_SLC[plot_ind_778], type = "l", ylim = c(0.30, 0.44))
lines(test_sol_SL_Dmat_1[plot_ind_778], type = "l", col = "blue")
abline(v = 101)
# Stray light seems to make a big difference in the 778  nm region. Examine stray light:
plot(sol_SL_magnitude_Dmat_1[plot_ind_778], type = "l")
abline(v = 101)
plot(veg_SL_magnitude_Dmat_1[plot_ind_778], type = "l")
abline(v = 101)
# plot whole spectrum of stray light
plot(sol_SL_magnitude_Dmat_1, type = "l")
abline(v = window_mins_mat[12,1])
abline(v = left_shoulder[12,1], lty = "dotted")
plot(veg_SL_magnitude_Dmat_1, type = "l")
abline(v = window_mins_mat[12,1])
abline(v = left_shoulder[12,1], lty = "dotted")
# conclusion: stray light adds large error near 778, and true SIF is somewhat lower. Negative sign can result if Ein > Eout. Combined, this makes the error large.

# What about for 3FLD Dmat_1 Ein versus Eout for Dmat_1 stray light level? 
# Is the right shoulder value, or left shoulder value, or shoulder mean less than the absorption feature minimum?
veg_SL_Dmat_1_SIF[left_shoulder[ ,1]] > veg_SL_Dmat_1_SIF[window_mins_mat[ ,1]]
veg_SL_Dmat_1_SIF[right_shoulder[ ,1]] > veg_SL_Dmat_1_SIF[window_mins_mat[ ,1]]
sol_SL_Dmat_1[left_shoulder[ ,1]] > sol_SL_Dmat_1[window_mins_mat[ ,1]]
sol_SL_Dmat_1[right_shoulder[ ,1]] > sol_SL_Dmat_1[window_mins_mat[ ,1]]
# For SIF level 3, all shoulders of veg spectra are greater than the absorption feature min.
# For SIF level 1, all shoulders of veg spectra are greater than the absorption feature min.
# For SIF level 1, all shoulders of the reference spectra are greater than the absorption feature min.
# Conclusion: stray light does not distort the veg spectra or solar spectra to the extent that shoulders are less than valley


# Plot to understand why SIF level 1 has negative bias for retrieval (this plot, and others, will differ based on SIF level setting)
ggplot() +
  geom_line(data = data.frame(x = xVals[1:left_shoulder[1,1]+20], y = L_veg_SIF[1:left_shoulder[1,1]+20]),
            aes(x=x, y=y, colour = "no stray"), alpha = 0.85, size = 0.6) +
  geom_point(data = data.frame (x = xVals[right_shoulder[1,1]], y = L_veg_SIF[right_shoulder[1,1]]),
             aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1) +
  # Highest stray light scenario
  geom_line(data = data.frame(x = xVals[1:left_shoulder[1,1]+20], y = veg_SL_Dmat_1_SIF[1:left_shoulder[1,1]+20]), 
            aes(x=x, y=y, colour = "stray"), alpha = 0.85, size = 0.6) +
  geom_point(data = data.frame(x = xVals[right_shoulder[1, 1]], y = veg_SL_Dmat_1_SIF[right_shoulder[1, 1]]),
             aes(x=x, y=y, colour = "stray"), alpha = 0.85, size = 1) +
  geom_point(data = data.frame(x = xVals[left_shoulder[1, 1]], y = veg_SL_Dmat_1_SIF[left_shoulder[1, 1]]),
             aes(x=x, y=y, colour = "stray"), alpha = 0.85, size = 1) +
  geom_point(data = data.frame(x = xVals[window_mins_mat[1, 1]], y = veg_SL_Dmat_1_SIF[window_mins_mat[1, 1]]),
             aes(x=x, y=y, colour = "stray"), alpha = 0.85, size = 1)
# Not very useful plot  because hard to see relative shapes of the spectra. Normalizing to same scale is more useful.


# Make function to calculate the component of the 3FLD equation that represents the R flux. This is useful in plots below.
# Could move this to script of functions
R_deduced_3FLD <- function(lambdaL, lambdaR, lambdaIN, LoutL, LoutR, EoutL, EoutR, Ein){
  # See Damm et al 2011 RSE p. 1885 for description of this component of 3FLD retrival equation
  w_21 <- (lambdaR - lambdaIN)/(lambdaR-lambdaL)
  w_22 <- (lambdaIN - lambdaL)/(lambdaR - lambdaL)
  R_deduced <- (Ein/(w_21*EoutL+w_22*EoutR)) * (w_21*LoutL + w_22*LoutR)
}

# Make function to calculate the denominator used in the 3FLD equation. This is useful in the plots below.
# Could move this to script of functions
threeFLD_denom <- function(lambdaL, lambdaR, lambdaIN, EoutL, EoutR, Ein){
  # See Damm et al 2011 RSE p. 1885 for description of this component of 3FLD retrival equation
  w_21 <- (lambdaR - lambdaIN)/(lambdaR-lambdaL)
  w_22 <- (lambdaIN - lambdaL)/(lambdaR - lambdaL)
  denom <- 1 - (Ein/(w_21*EoutL+w_22*EoutR))
  return(denom)
}



colsSpecStray <- c("stray"="#DAC228","no stray"="#122B5D","Simulated SIF"="black")

# # Example plot for highest stray light scenario: veg spectrum
# ggplot() +
#   geom_line(data = data.frame(x = xVals[1:left_shoulder[1,1]+20], 
#                               y = L_veg_SIF[1:left_shoulder[1,1]+20]/max(L_veg_SIF[1:left_shoulder[1,1]+20])),
#             aes(x=x, y=y, colour = "no stray"), alpha = 0.85, size = 0.6) +
#   geom_point(data = data.frame (x = xVals[right_shoulder[1,1]], 
#                                 y = L_veg_SIF[right_shoulder[1,1]]/max(L_veg_SIF[1:left_shoulder[1,1]+20])),
#              aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1.5) +
#   geom_point(data = data.frame (x = xVals[left_shoulder[1,1]], 
#                                 y = L_veg_SIF[left_shoulder[1,1]]/max(L_veg_SIF[1:left_shoulder[1,1]+20])),
#              aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1.5) +
#   geom_point(data = data.frame (x = xVals[window_mins_mat[1,1]], 
#                                 y = L_veg_SIF[window_mins_mat[1,1]]/max(L_veg_SIF[1:left_shoulder[1,1]+20])),
#              aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1.5) +
#   geom_point(data = data.frame (x = xVals[window_mins_mat[1,1]], 
#                                 y = R_deduced_no_SL/max(L_veg_SIF[1:left_shoulder[1,1]+20])),
#              aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1.5, shape = 15) +
#   # Highest stray light scenario
#   geom_line(data = data.frame(x = xVals[1:left_shoulder[1,1]+20], 
#                               y = veg_SL_Dmat_1_SIF[1:left_shoulder[1,1]+20]/max(veg_SL_Dmat_1_SIF[1:left_shoulder[1,1]+20])), 
#             aes(x=x, y=y, colour = "stray"), alpha = 0.85, size = 0.6) +
#   geom_point(data = data.frame (x = xVals[window_mins_mat[1,1]], 
#                                 y = veg_SL_Dmat_1_SIF[window_mins_mat[1,1]]/max(veg_SL_Dmat_1_SIF[1:left_shoulder[1,1]+20])),
#              aes(x = x, y = y, colour = "stray"), alpha = 0.85, size = 1.5) +
#   geom_point(data = data.frame (x = xVals[window_mins_mat[1,1]], 
#                                 y = R_deduced_Dmat_1/max(veg_SL_Dmat_1_SIF[1:left_shoulder[1,1]+20])),
#              aes(x = x, y = y, colour = "stray"), alpha = 0.85, size = 1.5, shape = 15) +
#   geom_text(data = data.frame (x = xVals[window_mins_mat[1,1]], 
#                                y = R_deduced_Dmat_1/max(veg_SL_Dmat_1_SIF[1:left_shoulder[1,1]+20])),
#             aes(x = x, y = y, colour = "stray"),
#             label = sprintf("L1 - L[deduced] == '%.3e'", 
#                           (veg_SL_Dmat_1_SIF[window_mins_mat[1,1]] - R_deduced_Dmat_1)/max(veg_SL_Dmat_1_SIF[1:left_shoulder[1,1]+20])),
#             nudge_x = -0.6,
#             size = 3.5,
#             show.legend = FALSE,
#             parse = TRUE) +
#   xlab("Wavelength (nm)") +
#   scale_colour_manual(values=colsSpecStray, name = "")

# Idea is to give diff between L1 and deduced R.   

  


# Function for version of plot normalized to each max within the interval annotated
make_norm_plot_show_bias <- function(x_range, 
                                     y_range_no_SL, y_max_no_SL, 
                                     x_left_shoulder_pt, y_left_shoulder_pt_no_SL, 
                                     x_right_shoulder_pt, y_right_shoulder_pt_no_SL, 
                                     x_win_min, y_win_min_no_SL,
                                     R_ded_no_SL,
                                     denominator_no_SL,
                                     y_range_SL, y_max_SL, 
                                     y_left_shoulder_pt_SL,
                                     y_right_shoulder_pt_SL, 
                                     y_win_min_SL,
                                     R_ded_SL,
                                     denominator_SL){
  ggplot() +
    # No SL scenario
    geom_line(data = data.frame(x = x_range, y = y_range_no_SL/y_max_no_SL),
              aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 0.6) +
    geom_point(data = data.frame (x = x_left_shoulder_pt, y = y_left_shoulder_pt_no_SL/y_max_no_SL),
               aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1) +
    geom_point(data = data.frame (x = x_right_shoulder_pt, y = y_right_shoulder_pt_no_SL/y_max_no_SL),
               aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1) +    
    geom_point(data = data.frame (x = x_win_min, y = y_win_min_no_SL/y_max_no_SL),
               aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1) +
    geom_point(data = data.frame (x = x_win_min, y = R_ded_no_SL/y_max_no_SL),
               aes(x = x, y = y, colour = "no stray"), alpha = 0.85, size = 1.5, shape = 15) +
    geom_text(data = data.frame (x = x_win_min, y = R_ded_no_SL/y_max_no_SL),
              aes(x = x, y = y, colour = "no stray"),
              label = sprintf("L1 - L[ded.] == '%.3e'", 
                              (y_win_min_no_SL - R_ded_no_SL)),
              nudge_x = -0.45,
              nudge_y = -0.00,
              size = 3,
              show.legend = FALSE,
              parse = TRUE) +
    geom_text(data = data.frame (x = x_win_min, y = R_ded_no_SL/y_max_no_SL),
              aes(x = x, y = y, colour = "no stray"),
              label = sprintf("Denominator == '%.3e'",
                              denominator_no_SL),
              nudge_x = 0.5, #.45 works for most panels
              nudge_y = -0.00,
              size = 3,
              show.legend = FALSE,
              parse = TRUE) +
    # Stray light scenario
    geom_line(data = data.frame(x = x_range, y = y_range_SL/y_max_SL), 
              aes(x=x, y=y, colour = "stray"), alpha = 0.85, size = 0.6) +
    geom_point(data = data.frame (x = x_left_shoulder_pt, y = y_left_shoulder_pt_SL/y_max_SL),
               aes(x = x, y = y, colour = "stray"), alpha = 0.85, size = 1) +
    geom_point(data = data.frame (x = x_right_shoulder_pt, y = y_right_shoulder_pt_SL/y_max_SL),
               aes(x = x, y = y, colour = "stray"), alpha = 0.85, size = 1) +  
    geom_point(data = data.frame (x = x_win_min, y = y_win_min_SL/y_max_SL),
               aes(x = x, y = y, colour = "stray"), alpha = 0.85, size = 1.5, shape = 17) +
    geom_point(data = data.frame (x = x_win_min, y = R_ded_SL/y_max_SL),
               aes(x = x, y = y, colour = "stray"), alpha = 0.85, size = 1.5, shape = 15) +
    geom_text(data = data.frame (x = x_win_min, y = R_ded_SL/y_max_SL),
              aes(x = x, y = y, colour = "stray"),
              label = sprintf("L1 - L[ded.] == '%.3e'", 
                              (y_win_min_SL - R_ded_SL)),
              nudge_x = -0.45,
              nudge_y = -0.006,
              size = 3,
              show.legend = FALSE,
              parse = TRUE) +
    geom_text(data = data.frame (x = x_win_min, y = R_ded_SL/y_max_SL),
              aes(x = x, y = y, colour = "stray"),
              label = sprintf("Denominator == '%.3e'",
                              denominator_SL),
              nudge_x = 0.5, #.45 works for most panels
              nudge_y = -0.006,
              check_overlap = TRUE,
              size = 3,
              show.legend = FALSE,
              parse = TRUE) +
    # Other arguments
    theme_classic() +
    xlab("Wavelength (nm)") +
    ylab("Radiance normalized to max within interval") +
    scale_colour_manual(values=colsSpecStray, name = "") +
    theme(legend.position="none")
}



# Highest Stray light scenario: Use function to plot all 12 features for stray light versus non-stray light
for(k in 1:dim(window_mins_mat)[1]){
  
  R_deduced_no_SL <- R_deduced_3FLD(lambdaL = left_shoulder[k, 2],
                                    lambdaR = right_shoulder[k, 2],
                                    lambdaIN = window_mins_mat[k, 2],
                                    LoutL = L_veg_SIF[left_shoulder[k, 1]],
                                    LoutR = L_veg_SIF[right_shoulder[k, 1]],
                                    EoutL = SolarData_SLC[left_shoulder[k, 1]],
                                    EoutR = SolarData_SLC[right_shoulder[k, 1]],
                                    Ein = SolarData_SLC[window_mins_mat[k, 1]])
  
  R_deduced_Dmat_1 <- R_deduced_3FLD(lambdaL = left_shoulder[k, 2],
                                     lambdaR = right_shoulder[k, 2],
                                     lambdaIN = window_mins_mat[k, 2],
                                     LoutL = veg_SL_Dmat_1_SIF[left_shoulder[k, 1]],
                                     LoutR = veg_SL_Dmat_1_SIF[right_shoulder[k, 1]],
                                     EoutL = sol_SL_Dmat_1[left_shoulder[k, 1]],
                                     EoutR = sol_SL_Dmat_1[right_shoulder[k, 1]],
                                     Ein = sol_SL_Dmat_1[window_mins_mat[k, 1]])
  
  denom_no_SL <- threeFLD_denom(lambdaL = left_shoulder[k, 2],
                                lambdaR = right_shoulder[k, 2],
                                lambdaIN = window_mins_mat[k, 2],
                                EoutL = SolarData_SLC[left_shoulder[k, 1]],
                                EoutR = SolarData_SLC[right_shoulder[k, 1]],
                                Ein = SolarData_SLC[window_mins_mat[k, 1]])
  
  denom_Dmat_1 <- threeFLD_denom(lambdaL = left_shoulder[k, 2],
                             lambdaR = right_shoulder[k, 2],
                             lambdaIN = window_mins_mat[k, 2],
                             EoutL = sol_SL_Dmat_1[left_shoulder[k, 1]],
                             EoutR = sol_SL_Dmat_1[right_shoulder[k, 1]],
                             Ein = sol_SL_Dmat_1[window_mins_mat[k, 1]])
  
  veg_win_min_Dmat_1_fig <- make_norm_plot_show_bias(
    x_range = xVals[(left_shoulder[k,1]-12):(right_shoulder[k,1]+12)],
    y_range_no_SL = L_veg_SIF[(left_shoulder[k,1]-12):(right_shoulder[k,1]+12)],
    y_max_no_SL = max(L_veg_SIF[(left_shoulder[k,1]-12):(right_shoulder[k,1]+12)]),
    x_left_shoulder_pt = xVals[left_shoulder[k,1]],
    y_left_shoulder_pt_no_SL = L_veg_SIF[left_shoulder[k,1]],
    x_right_shoulder_pt = xVals[right_shoulder[k,1]],
    y_right_shoulder_pt_no_SL = L_veg_SIF[right_shoulder[k,1]],
    x_win_min = xVals[window_mins_mat[k,1]],
    y_win_min_no_SL = L_veg_SIF[window_mins_mat[k,1]],
    R_ded_no_SL = R_deduced_no_SL,
    denominator_no_SL = denom_no_SL,
    y_range_SL = veg_SL_Dmat_1_SIF[(left_shoulder[k,1]-12):(right_shoulder[k,1]+12)],
    y_max_SL = max(veg_SL_Dmat_1_SIF[(left_shoulder[k,1]-12):(right_shoulder[k,1]+12)]),
    y_left_shoulder_pt_SL = veg_SL_Dmat_1_SIF[left_shoulder[k,1]],
    y_right_shoulder_pt_SL = veg_SL_Dmat_1_SIF[right_shoulder[k, 1]],
    y_win_min_SL = veg_SL_Dmat_1_SIF[window_mins_mat[k, 1]],
    R_ded_SL = R_deduced_Dmat_1,
    denominator_SL = denom_Dmat_1)
  
  ggsave(filename = paste(path_figs, "Fig_veg_win_min_", window_mins_mat[k,2], "_Dmat_1_SIFl",SIF_level,".pdf",sep=""),
         plot = veg_win_min_Dmat_1_fig, width = 3.9, height = 3.9, units = "in", dpi = 600)
  
}




### MAKE A TIDY DATA FRAME OF RESULTS FOR ANALYSIS #################################################################################################

## Make the results of each simulation into a data frame with these columns: retrieved SIF; Error/SIF; Retrieval type; Left shoulder stray light; 
# right shoulder stray light; left shoulder reflectance; right shoulder reflectance


three_FLD_veg_SL_Dmat_1_df <- data.frame(three_FLD_veg_SL_Dmat_1_SIF,
                                         three_FLD_veg_SL_Dmat_1_SIF_error,
                                         rep("3FLD", length.out = length(three_FLD_veg_SL_Dmat_1_SIF)),
                                         rep("Dmat_simulated", length.out = length(three_FLD_veg_SL_Dmat_1_SIF)),
                                         sol_SL_magnitude_Dmat_1[left_shoulder[ ,1]],
                                         sol_SL_magnitude_Dmat_1[right_shoulder[ ,1]],
                                         leaf_R[left_shoulder[ ,1]],
                                         leaf_R[right_shoulder[ ,1]],
                                         SolarData_SLC[window_mins_ind],
                                         SolarData_SLC[left_shoulder[ ,1]],
                                         SolarData_SLC[right_shoulder[ ,1]])

three_FLD_veg_SL_Dmat_2_df <- data.frame(three_FLD_veg_SL_Dmat_2_SIF,
                                         three_FLD_veg_SL_Dmat_2_SIF_error,
                                         rep("3FLD", length.out = length(three_FLD_veg_SL_Dmat_2_SIF)),
                                         rep("Dmat_simulated", length.out = length(three_FLD_veg_SL_Dmat_2_SIF)),
                                         sol_SL_magnitude_Dmat_2[left_shoulder[ ,1]],
                                         sol_SL_magnitude_Dmat_2[right_shoulder[ ,1]],
                                         leaf_R[left_shoulder[ ,1]],
                                         leaf_R[right_shoulder[ ,1]],
                                         SolarData_SLC[window_mins_ind],
                                         SolarData_SLC[left_shoulder[ ,1]],
                                         SolarData_SLC[right_shoulder[ ,1]])

three_FLD_veg_SL_Dmat_3_df <- data.frame(three_FLD_veg_SL_Dmat_3_SIF,
                                         three_FLD_veg_SL_Dmat_3_SIF_error,
                                         rep("3FLD", length.out = length(three_FLD_veg_SL_Dmat_3_SIF)),
                                         rep("Dmat_simulated", length.out = length(three_FLD_veg_SL_Dmat_1_SIF)),
                                         sol_SL_magnitude_Dmat_3[left_shoulder[ ,1]],
                                         sol_SL_magnitude_Dmat_3[right_shoulder[ ,1]],
                                         leaf_R[left_shoulder[ ,1]],
                                         leaf_R[right_shoulder[ ,1]],
                                         SolarData_SLC[window_mins_ind],
                                         SolarData_SLC[left_shoulder[ ,1]],
                                         SolarData_SLC[right_shoulder[ ,1]])

three_FLD_veg_SL_Dmat_4_df <- data.frame(three_FLD_veg_SL_Dmat_4_SIF,
                                         three_FLD_veg_SL_Dmat_4_SIF_error,
                                         rep("3FLD", length.out = length(three_FLD_veg_SL_Dmat_4_SIF)),
                                         rep("Dmat_simulated", length.out = length(three_FLD_veg_SL_Dmat_1_SIF)),
                                         sol_SL_magnitude_Dmat_4[left_shoulder[ ,1]],
                                         sol_SL_magnitude_Dmat_4[right_shoulder[ ,1]],
                                         leaf_R[left_shoulder[ ,1]],
                                         leaf_R[right_shoulder[ ,1]],
                                         SolarData_SLC[window_mins_ind],
                                         SolarData_SLC[left_shoulder[ ,1]],
                                         SolarData_SLC[right_shoulder[ ,1]])
 
 


sFLD_veg_SL_Dmat_1_df <- data.frame(sFLD_veg_SL_Dmat_1_SIF, 
                                    sFLD_veg_SL_Dmat_1_SIF_error, 
                                    rep("sFLD", length.out = length(sFLD_veg_SL_Dmat_1_SIF)),
                                    rep("Dmat_simulated", length.out = length(three_FLD_veg_SL_Dmat_1_SIF)),
                                    sol_SL_magnitude_Dmat_1[left_shoulder[ ,1]],
                                    sol_SL_magnitude_Dmat_1[right_shoulder[ ,1]],
                                    leaf_R[left_shoulder[ ,1]],
                                    leaf_R[right_shoulder[ ,1]],
                                    SolarData_SLC[window_mins_ind],
                                    SolarData_SLC[left_shoulder[ ,1]],
                                    SolarData_SLC[right_shoulder[ ,1]])

sFLD_veg_SL_Dmat_2_df <- data.frame(sFLD_veg_SL_Dmat_2_SIF, 
                                    sFLD_veg_SL_Dmat_2_SIF_error, 
                                    rep("sFLD", length.out = length(sFLD_veg_SL_Dmat_2_SIF)),
                                    rep("Dmat_simulated", length.out = length(three_FLD_veg_SL_Dmat_1_SIF)),
                                    sol_SL_magnitude_Dmat_2[left_shoulder[ ,1]],
                                    sol_SL_magnitude_Dmat_2[right_shoulder[ ,1]],
                                    leaf_R[left_shoulder[ ,1]],
                                    leaf_R[right_shoulder[ ,1]],
                                    SolarData_SLC[window_mins_ind],
                                    SolarData_SLC[left_shoulder[ ,1]],
                                    SolarData_SLC[right_shoulder[ ,1]])

sFLD_veg_SL_Dmat_3_df <- data.frame(sFLD_veg_SL_Dmat_3_SIF, 
                                    sFLD_veg_SL_Dmat_3_SIF_error, 
                                    rep("sFLD", length.out = length(sFLD_veg_SL_Dmat_3_SIF)),
                                    rep("Dmat_simulated", length.out = length(three_FLD_veg_SL_Dmat_1_SIF)),
                                    sol_SL_magnitude_Dmat_3[left_shoulder[ ,1]],
                                    sol_SL_magnitude_Dmat_3[right_shoulder[ ,1]],
                                    leaf_R[left_shoulder[ ,1]],
                                    leaf_R[right_shoulder[ ,1]],
                                    SolarData_SLC[window_mins_ind],
                                    SolarData_SLC[left_shoulder[ ,1]],
                                    SolarData_SLC[right_shoulder[ ,1]])

sFLD_veg_SL_Dmat_4_df <- data.frame(sFLD_veg_SL_Dmat_4_SIF, 
                                    sFLD_veg_SL_Dmat_4_SIF_error, 
                                    rep("sFLD", length.out = length(sFLD_veg_SL_Dmat_4_SIF)),
                                    rep("Dmat_simulated", length.out = length(three_FLD_veg_SL_Dmat_1_SIF)),
                                    sol_SL_magnitude_Dmat_4[left_shoulder[ ,1]],
                                    sol_SL_magnitude_Dmat_4[right_shoulder[ ,1]],
                                    leaf_R[left_shoulder[ ,1]],
                                    leaf_R[right_shoulder[ ,1]],
                                    SolarData_SLC[window_mins_ind],
                                    SolarData_SLC[left_shoulder[ ,1]],
                                    SolarData_SLC[right_shoulder[ ,1]])

# Make column names consistent and then bind
cols <- c("SIF","stray_error_over_SIF","retrieval_type", "stray_light_type", "sol_SL_left", "sol_SL_right", "leaf_R_left", "leaf_R_right", "sol_min", "sol_left", "sol_right")

colnames(three_FLD_veg_SL_Dmat_1_df) <- cols
colnames(three_FLD_veg_SL_Dmat_2_df) <- cols
colnames(three_FLD_veg_SL_Dmat_3_df) <- cols
colnames(three_FLD_veg_SL_Dmat_4_df) <- cols

colnames(sFLD_veg_SL_Dmat_1_df) <- cols
colnames(sFLD_veg_SL_Dmat_2_df) <- cols
colnames(sFLD_veg_SL_Dmat_3_df) <- cols
colnames(sFLD_veg_SL_Dmat_4_df) <- cols


all_FLD_veg_SL_df <- rbind(
                           three_FLD_veg_SL_Dmat_1_df,
                           three_FLD_veg_SL_Dmat_2_df,
                           three_FLD_veg_SL_Dmat_3_df,
                           three_FLD_veg_SL_Dmat_4_df,
                           sFLD_veg_SL_Dmat_1_df,
                           sFLD_veg_SL_Dmat_2_df,
                           sFLD_veg_SL_Dmat_3_df,
                           sFLD_veg_SL_Dmat_4_df)

all_FLD_veg_SL_df$window_mins <- window_mins

all_FLD_veg_SL_df$delta_stray <- all_FLD_veg_SL_df$sol_SL_left - all_FLD_veg_SL_df$sol_SL_right

all_FLD_veg_SL_df$delta_leaf_R <- all_FLD_veg_SL_df$leaf_R_left - all_FLD_veg_SL_df$leaf_R_right

all_FLD_veg_SL_df$line_depth_left <- all_FLD_veg_SL_df$sol_left - all_FLD_veg_SL_df$sol_min

all_FLD_veg_SL_df$Abs_SIF_error_perc <- abs(all_FLD_veg_SL_df$stray_error_over_SIF) * 100




#### FIGURES FOR PUBLICATION ############################################################################################


### Define text for use in graph
radiance_text = expression('Radiance (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')
radiance_text_log = expression('Log Radiance (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')
radiance_text_error = expression('Stray light error (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')
radiance_text_error = expression(atop('Stray light error', paste('(W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')))


### Plot of solar spectrum, SIF spectrum, and simulated vegetation spectrum based on L = E * R
# Useful: https://stackoverflow.com/questions/17148679/construct-a-manual-legend-for-a-complicated-plot
# Can't figure out how to remove the legend, except that 'color' can be defined outside of aes(), then scale_colour_manual line removed
# Also, if theme() is before theme_classic, then the legend commands seem to be overridden

cols <- c("Solar radiance"="black","Vegetation radiance"="#4C1495","SIF * 50"="#EE604A")

Fig_L_E_SIF <- ggplot() +
  # Mean measured solar spectrum
  geom_line(data = data.frame(x = xVals, y = SolarData_SLC), aes(x=x, y=y, colour = "Solar radiance"), alpha = 0.85, size = 0.6) +
  # Simulated L = E * R
  geom_line(data = data.frame(x = xVals, y = L_veg_SIF), aes(x=x, y=y, colour = "Vegetation radiance"), alpha = 0.85, size = 0.6) +
  # SIF spectrum * 50
  geom_line(data = data.frame(x = xVals, y = SIF_2perc * 50), aes(x=x, y=y, colour = "SIF * 50"), alpha = 0.85, size = 0.6) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=cols, name = "") +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        #text = element_text(size=20), #Makes everything text size 20
        legend.text = element_text(size = 14),
        legend.position="top", 
        legend.direction="horizontal") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_L_E_SIF_SIFl",SIF_level,".pdf", sep = ""), 
       plot = Fig_L_E_SIF, width = 8, height = 4, units = "in", dpi = 300)


# ### Plot of native spectral stray light for solar spectrum
# # Also see version of this figure from script called make_SIF_vs_stray_light_magnitude_Fig (that was version used for manuscript)
# 
# cols_native_SL <- c("SL in solar"="black","SL in vegetation"="grey55","SIF"="#EE604A")
# 
# Fig_original_SL_SIF <- ggplot() +
# #  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "SIF"), alpha = 0.85, size = 0.6) +
#   geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_orig), aes(x=x, y=y, colour = "SL in solar"), alpha = 0.85, size = 0.6) +
#   geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_orig), aes(x=x, y=y, colour = "SL in vegetation"), alpha = 0.85, size = 0.6) +
#   theme_classic() +
#   scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
#   xlab("Wavelength (nm)") +
#   ylab(radiance_text) +
#   scale_colour_manual(values=cols_native_SL, name = "") +
#   theme(axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         axis.text.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         legend.position="top", 
#         legend.direction="horizontal") + 
#   geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)
# 
# ggsave(filename = paste(path_figs, "Fig_original_SL_SIF_SIFl",SIF_level,".pdf", sep = ""), 
#        plot = Fig_original_SL_SIF, width = 8, height = 4, units = "in", dpi = 300)
  

### Plot of leaf reflectance spectrum
Fig_Leaf_R <- ggplot() +
  geom_line(data = data.frame(x = xVals, y = Leaf_R_resampled$y), aes(x = x, y = y), color = "black", alpha = 0.9, size = 0.6) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 20, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  xlab("Wavelength (nm)") +
  ylab("Leaf reflectance")

ggsave(filename = paste(path_figs, "Fig_Leaf_R.pdf", sep = ""), plot = Fig_Leaf_R, width = 8, height = 4, units = "in", dpi = 300)



### Plot of D_MATRIX simulated stray light for vegetation spectrum
stray_cols <- c("stray_1"="#0868ac", "stray_2"="#43a2ca","stray_3"="#7bccc4", "stray_4"="#bae4bc")

veg_strayL_Dmat <- ggplot() +
  geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_Dmat_1), aes(x=x, y=y, colour = "stray_1"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_Dmat_2), aes(x=x, y=y, colour = "stray_2"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_Dmat_3), aes(x=x, y=y, colour = "stray_3"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_Dmat_4), aes(x=x, y=y, colour = "stray_4"), alpha = 0.8, size = 0.6) +
  scale_y_log10() +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text_log) +
  scale_colour_manual(values=stray_cols, name = "") +
  geom_vline(xintercept = window_mins, alpha = 0.6, size = 0.4) +
  geom_vline(xintercept = left_shoulder[ ,2], alpha = 0.6, size = 0.4, linetype = "dashed") +
  geom_vline(xintercept = right_shoulder[ ,2], alpha = 0.6, size = 0.4, linetype = "dotted") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(filename = paste(path_outputs, "Fig_veg_strayL_Dmat_SIFl",SIF_level,".pdf", sep = ""),
       plot = veg_strayL_Dmat, width = 8, height = 4, units = "in", dpi = 300)


### Plot of D_MATRIX simulated stray light for solar spectrum

sol_strayL_Dmat <- ggplot() +
  geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_Dmat_1), aes(x=x, y=y, colour = "stray_1"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_Dmat_2), aes(x=x, y=y, colour = "stray_2"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_Dmat_3), aes(x=x, y=y, colour = "stray_3"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_Dmat_4), aes(x=x, y=y, colour = "stray_4"), alpha = 0.8, size = 0.6) +
  scale_y_log10() +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text_log) +
  scale_colour_manual(values=stray_cols, name = "") +
  geom_vline(xintercept = window_mins, alpha = 0.6, size = 0.4) +
  geom_vline(xintercept = left_shoulder[ ,2], alpha = 0.6, size = 0.4, linetype = "dashed") +
  geom_vline(xintercept = right_shoulder[ ,2], alpha = 0.6, size = 0.4, linetype = "dotted") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank())

ggsave(filename = paste(path_figs, "Fig_sol_strayL_Dmat.pdf", sep = ""), plot = sol_strayL_Dmat, width = 8, height = 4, units = "in", dpi = 300)


### Plots of D Matrix simulated stray light as a proportion of solar signal

sol_strayL_Dmat_proportion_sig <- ggplot() +
  geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_Dmat_1/SolarData_SLC), aes(x=x, y=y, colour = "stray_1"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_Dmat_2/SolarData_SLC), aes(x=x, y=y, colour = "stray_2"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_Dmat_3/SolarData_SLC), aes(x=x, y=y, colour = "stray_3"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = sol_SL_magnitude_Dmat_4/SolarData_SLC), aes(x=x, y=y, colour = "stray_4"), alpha = 0.8, size = 0.6) +
  scale_y_log10() +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Wavelength (nm)") +
  ylab("Proportion of signal") +
  scale_colour_manual(values=stray_cols, name = "") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank())

ggsave(filename = paste(path_figs,"Fig_sol_strayL_Dmat_relative_signal.pdf", sep = ""), plot = sol_strayL_Dmat_proportion_sig, width = 8, height = 4, units = "in", dpi = 300)


veg_strayL_Dmat_proportion_sig <- ggplot()+
  geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_Dmat_1/L_veg_SIF), aes(x=x, y=y, colour = "stray_1"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_Dmat_2/L_veg_SIF), aes(x=x, y=y, colour = "stray_2"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_Dmat_3/L_veg_SIF), aes(x=x, y=y, colour = "stray_3"), alpha = 0.8, size = 0.6) +
  geom_line(data = data.frame(x = xVals, y = veg_SL_magnitude_Dmat_4/L_veg_SIF), aes(x=x, y=y, colour = "stray_4"), alpha = 0.8, size = 0.6) +
  scale_y_log10() +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Wavelength (nm)") +
  ylab("Proportion of signal") +
  scale_colour_manual(values=stray_cols, name = "") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank())

ggsave(filename = paste(path_outputs, "Fig_veg_SL_Dmat_relative_signal_SIFl",SIF_level,".pdf",sep=""),
       plot = veg_strayL_Dmat_proportion_sig, width = 8, height = 4, units = "in", dpi = 300)

### Plot of retrieval windows

# Red region
red_region <- c(670, 689)
red_region_y <- c(0.16, 0.41)
rect_df <- data.frame(x1=left_shoulder[1:4 ,2], x2=right_shoulder[1:4 ,2], 
             y1 = rep(red_region_y[1], length.out = 4), y2 = rep(red_region_y[2], length.out = 4))
sol_red_windows <- ggplot() +
  geom_line(data = SolarData_SLC_df, aes(x = spec_pixel_pred_no710, y = SolarData_SLC), color = "black", alpha = 0.85, size = 0.6) +
  scale_x_continuous(limits = red_region, expand = c(0, 0)) +
  scale_y_continuous(limits = red_region_y, expand = c(0,0)) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  geom_rect(data = rect_df, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color=NA, alpha=0.2) +
  geom_vline(xintercept = rect_df$x1, linetype="dashed", color = "black", size=0.5) +
  geom_vline(xintercept = rect_df$x2, linetype="dotted", color = "black", size=0.5) +
  geom_vline(xintercept = 687.1, linetype="solid", color = "#0072B2", size=0.5) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = unit(c(.25, .5, .25, .25), "cm"))   # top, right, bottom, left

suppressWarnings(ggsave(filename = paste(path_figs,"Fig_sol_red_windows.pdf", sep = ""), plot = sol_red_windows, width = 8, height = 4, units = "in", dpi = 300))
# Note, the warning about 1886 rows removed is because of the x and y axis limits selected to focus on red region of spectrum. Safe to ignore warning.


# far-red region up to 758
far_red_region <- c(742, 758)
far_red_region_y <- c(0.28, 0.38)
rect_df <- data.frame(x1=left_shoulder[5:9 ,2], x2=right_shoulder[5:9 ,2], 
                      y1 = rep(far_red_region_y[1], length.out = 5), y2 = rep(far_red_region_y[2], length.out = 5))
FR_rect_df <- data.frame(FR1 = xVals[1530], FR2 = xVals[1700], yFR1 = far_red_region_y[1], yFR2 = far_red_region_y[2])
sol_far_red_windows <- ggplot() +
  geom_line(data = SolarData_SLC_df, aes(x = spec_pixel_pred_no710, y = SolarData_SLC), color = "black", alpha = 0.85, size = 0.6) +
  scale_x_continuous(limits = far_red_region, expand = c(0, 0)) +
  scale_y_continuous(limits = far_red_region_y, expand = c(0,0),
                    breaks = seq(from = far_red_region_y[1], to = far_red_region_y[2], length.out = 6)) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), color=NA, alpha=0.2) +
  geom_rect(data = FR_rect_df, aes(xmin = FR1, xmax = FR2, ymin = yFR1, ymax = yFR2), color = NA, fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = rect_df$x1, linetype="dashed", color = "black", size=0.5) +
  geom_vline(xintercept = rect_df$x2, linetype="dotted", color = "black", size=0.5) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = unit(c(.25, .5, .25, .25), "cm"))   # top, right, bottom, left

suppressWarnings(ggsave(filename = paste(path_figs, "Fig_sol_far_red_windows.pdf", sep = ""), plot = sol_far_red_windows, width = 8, height = 4, units = "in", dpi = 300))
# Note, the warning about 1847 rows removed is because of the x and y limits chosen to focus on the far-red region of the spectrum. Safe to ignore this warning.

# far-red region from 770 to 780
farther_red_region <- c(770, 780)
farther_red_region_y <- c(0.28, 0.36)
rect_df <- data.frame(x1=left_shoulder[10:12 ,2], x2=right_shoulder[10:12 ,2], 
                      y1 = rep(farther_red_region_y[1], length.out = 3), y2 = rep(farther_red_region_y[2], length.out = 3))
sol_farther_red_windows <- ggplot() +
  geom_line(data = SolarData_SLC_df, aes(x = spec_pixel_pred_no710, y = SolarData_SLC), color = "black", alpha = 0.85, size = 0.6) +
  scale_x_continuous(limits = farther_red_region, expand = c(0, 0),
                     breaks = seq(from = farther_red_region[1], to = farther_red_region[2], length.out = 6)) +
  scale_y_continuous(limits = farther_red_region_y, expand = c(0,0)) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  geom_rect(data = rect_df, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color=NA, alpha=0.2) +
  geom_vline(xintercept = rect_df$x1, linetype="dashed", color = "black", size=0.5) +
  geom_vline(xintercept = rect_df$x2, linetype="dotted", color = "black", size=0.5) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = unit(c(.25, .5, .25, .25), "cm"))   # top, right, bottom, left

suppressWarnings(ggsave(filename = paste(path_figs,"Fig_sol_farther_red_windows.pdf", sep = ""), plot = sol_farther_red_windows, width = 8, height = 4, units = "in", dpi = 300))
# Note, the warning about 1965 rows removed is because of the x and y limits chosen to focus on this region of the spectrum. Safe to ignore this warning.

# O2A region
O2A_region <- c(757, 771)
O2A_region_y <- c(0.0, 0.35)
rect_df <- data.frame(x1=left_shoulder[13 ,2], x2=right_shoulder[13 ,2], 
                      y1 = rep(O2A_region_y[1], length.out = 1), y2 = rep(O2A_region_y[2], length.out = 1))
sol_O2A_windows <- ggplot() +
  geom_line(data = SolarData_SLC_df, aes(x = spec_pixel_pred_no710, y = SolarData_SLC), color = "black", alpha = 0.85, size = 0.6) +
  scale_x_continuous(limits = O2A_region, expand = c(0, 0)) +
  scale_y_continuous(limits = O2A_region_y, expand = c(0,0),
                     breaks = seq(from = O2A_region_y[1], to = O2A_region_y[2], length.out = 6)) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  geom_rect(data = rect_df, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color=NA, alpha=0.2) +
  geom_vline(xintercept = rect_df$x1, linetype="dashed", color = "black", size=0.5) +
  geom_vline(xintercept = rect_df$x2, linetype="dotted", color = "black", size=0.5) +
  geom_vline(xintercept = 760.5, linetype="solid", color = "#0072B2", size=0.5) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = unit(c(.25, .5, .25, .25), "cm"))   # top, right, bottom, left

suppressWarnings(ggsave(filename = paste(path_figs,"Fig_sol_O2A_windows.pdf", sep = ""), 
                        plot = sol_O2A_windows, width = 8, height = 4, units = "in", dpi = 300))
# Note, the warning about 1887 rows removed is because of the x and y limits chosen to focus on this region of the spectrum. Safe to ignore this warning.


### Plots of retrieved SIF with and without stray light for 3FLD

sif_ylim_sFLD <- c(-.002, max(sFLD_veg_SL_Dmat_1_SIF) + 0.0005) # This range works for all SIF level 2 sFLD and 3FLD retrievals without warnings (and  maybe other SIF levels)
sif_ylim <- c(-.0001, max(SIF_2perc) + 0.0005) # A minimum of -.0005 includes lower sFLD values

colsRetStray <- c("3FLD with stray"="#DAC228","3FLD no stray"="#122B5D","Simulated SIF"="black")


# 3FLD: Very high stray light scenario
SIF_ret_Dmat_1 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = three_FLD_L_veg_SIF), aes(x=x, y=y, colour = "3FLD no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = three_FLD_veg_SL_Dmat_1_SIF), aes(x=x, y=y, colour = "3FLD with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  # theme(legend.position="none") +
  theme(legend.position = c(0.21, 0.9),
        legend.direction="vertical") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_3FLD_SIF_ret_Dmat_1_SIFl", SIF_level,".pdf", sep = ""),
       plot = SIF_ret_Dmat_1, width = 10, height = 6, units = "cm", dpi = 300)


# 3FLD: High stray light scenario
SIF_ret_Dmat_2 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = three_FLD_L_veg_SIF), aes(x=x, y=y, colour = "3FLD no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = three_FLD_veg_SL_Dmat_2_SIF), aes(x=x, y=y, colour = "3FLD with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_3FLD_SIF_ret_Dmat_2_SIFl",SIF_level,".pdf",sep = ""),
       plot = SIF_ret_Dmat_2, width = 10, height = 6, units = "cm", dpi = 300)
  

# 3FLD: Medium stray light scenario
SIF_ret_Dmat_3 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = three_FLD_L_veg_SIF), aes(x=x, y=y, colour = "3FLD no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = three_FLD_veg_SL_Dmat_3_SIF), aes(x=x, y=y, colour = "3FLD with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_3FLD_SIF_ret_Dmat_3_SIFl",SIF_level,".pdf",sep = ""),
       plot = SIF_ret_Dmat_3, width = 10, height = 6, units = "cm", dpi = 300)


# 3FLD: Low stray light scenario
SIF_ret_Dmat_4 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = three_FLD_L_veg_SIF), aes(x=x, y=y, colour = "3FLD no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = three_FLD_veg_SL_Dmat_4_SIF), aes(x=x, y=y, colour = "3FLD with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_3FLD_SIF_ret_Dmat_4_SIFl",SIF_level,".pdf",sep = ""),
       plot = SIF_ret_Dmat_4, width = 10, height = 6, units = "cm", dpi = 300)



### Plots of retrieved SIF with and without stray light for sFLD

# sFLD: high stray light scenario

sFLD_SIF_ret_Dmat_1 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = sFLD_L_veg_SIF), aes(x=x, y=y, colour = "3FLD no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = sFLD_veg_SL_Dmat_1_SIF), aes(x=x, y=y, colour = "3FLD with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim_sFLD) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_sFLD_SIF_ret_Dmat_1_SIFl",SIF_level,".pdf",sep=""),
       plot = sFLD_SIF_ret_Dmat_1, width = 10, height = 6, units = "cm", dpi = 300)


# sFLD: Medium stray light scenario
sFLD_SIF_ret_Dmat_2 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = sFLD_L_veg_SIF), aes(x=x, y=y, colour = "3FLD no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = sFLD_veg_SL_Dmat_2_SIF), aes(x=x, y=y, colour = "3FLD with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim_sFLD) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_sFLD_SIF_ret_Dmat_2_SIFl",SIF_level,".pdf",sep=""),
       plot = sFLD_SIF_ret_Dmat_2, width = 10, height = 6, units = "cm", dpi = 300)


# sFLD: low stray light scenario
sFLD_SIF_ret_Dmat_3 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = sFLD_L_veg_SIF), aes(x=x, y=y, colour = "3FLD no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = sFLD_veg_SL_Dmat_3_SIF), aes(x=x, y=y, colour = "3FLD with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim_sFLD) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_sFLD_SIF_ret_Dmat_3_SIFl",SIF_level,".pdf",sep=""),
       plot = sFLD_SIF_ret_Dmat_3, width = 10, height = 6, units = "cm", dpi = 300)


# sFLD: very low stray light scenario
sFLD_SIF_ret_Dmat_4 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_2perc), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = sFLD_L_veg_SIF), aes(x=x, y=y, colour = "3FLD no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light
  geom_point(data = data.frame(x = window_mins_mat[,2], y = sFLD_veg_SL_Dmat_4_SIF), aes(x=x, y=y, colour = "3FLD with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim_sFLD) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_sFLD_SIF_ret_Dmat_4_SIFl",SIF_level,".pdf",sep = ""),
       plot = sFLD_SIF_ret_Dmat_4, width = 10, height = 6, units = "cm", dpi = 300)



### Plots of stray light error for sFLD and 3FLD

colsRet <- c("sFLD"="#3C0E66", "3FLD"="#058279")

sif_error_ylim <- c(-0.0011, 0.0025)


# Stray errors for sFLD and 3FLD: High stray light scenario
Stray_error_Dmat_1 <- ggplot() +
  geom_line(data = data.frame(x = window_mins_mat[,2], y = (sFLD_veg_SL_Dmat_1_SIF - sFLD_L_veg_SIF)), aes(x=x, y=y, colour = "sFLD"), 
            alpha = 0.70, size = 1) +
  geom_point(data = data.frame(x = window_mins_mat[,2], y = (sFLD_veg_SL_Dmat_1_SIF - sFLD_L_veg_SIF)), aes(x=x, y=y, colour = "sFLD"), 
             alpha = 0.70, size = 3.3, shape = 15) +
  geom_line(data = data.frame(x = window_mins_mat[,2], y = (three_FLD_veg_SL_Dmat_1_SIF - three_FLD_L_veg_SIF)), aes(x=x, y=y, colour = "3FLD"), 
            alpha = 0.70, size = 1) +
  geom_point(data = data.frame(x = window_mins_mat[,2], y = (three_FLD_veg_SL_Dmat_1_SIF - three_FLD_L_veg_SIF)), aes(x=x, y=y, colour = "3FLD"), 
             alpha = 0.70, size = 4.5, shape = 18) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_error_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text_error) +
  scale_colour_manual(values=colsRet, name = "") +
  theme(legend.position = c(0.79, 0.9), 
        legend.direction="vertical") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_Stray_error_Dmat_1_SIFl",SIF_level,".pdf",sep = ""),
       plot = Stray_error_Dmat_1, width = 10, height = 6, units = "cm", dpi = 300)


# Stray errors for sFLD and 3FLD: Medium stray light scenario
Stray_error_Dmat_2 <- ggplot() +
  geom_line(data = data.frame(x = window_mins_mat[,2], y = (sFLD_veg_SL_Dmat_2_SIF - sFLD_L_veg_SIF)), aes(x=x, y=y, colour = "sFLD"), 
            alpha = 0.70, size = 1) +
  geom_point(data = data.frame(x = window_mins_mat[,2], y = (sFLD_veg_SL_Dmat_2_SIF - sFLD_L_veg_SIF)), aes(x=x, y=y, colour = "sFLD"), 
             alpha = 0.70, size = 3.3, shape = 15) +
  geom_line(data = data.frame(x = window_mins_mat[,2], y = (three_FLD_veg_SL_Dmat_2_SIF - three_FLD_L_veg_SIF)), aes(x=x, y=y, colour = "3FLD"), 
            alpha = 0.70, size = 1) +
  geom_point(data = data.frame(x = window_mins_mat[,2], y = (three_FLD_veg_SL_Dmat_2_SIF - three_FLD_L_veg_SIF)), aes(x=x, y=y, colour = "3FLD"), 
             alpha = 0.70, size = 4.5, shape = 18) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_error_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text_error) +
  scale_colour_manual(values=colsRet, name = "") +
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_Stray_error_Dmat_2_SIFl",SIF_level,".pdf",sep = ""),
       plot = Stray_error_Dmat_2, width = 10, height = 6, units = "cm", dpi = 300)


# Stray light errors for sFLD and 3FLD: Low stray light scenario
Stray_error_Dmat_3 <- ggplot() +
  geom_line(data = data.frame(x = window_mins_mat[,2], y = (sFLD_veg_SL_Dmat_3_SIF - sFLD_L_veg_SIF)), aes(x=x, y=y, colour = "sFLD"), 
            alpha = 0.75, size = 1) +
  geom_point(data = data.frame(x = window_mins_mat[,2], y = (sFLD_veg_SL_Dmat_3_SIF - sFLD_L_veg_SIF)), aes(x=x, y=y, colour = "sFLD"), 
             alpha = 0.75, size = 3.3, shape = 15) +
  geom_line(data = data.frame(x = window_mins_mat[,2], y = (three_FLD_veg_SL_Dmat_3_SIF - three_FLD_L_veg_SIF)), aes(x=x, y=y, colour = "3FLD"), 
            alpha = 0.75, size = 1) +
  geom_point(data = data.frame(x = window_mins_mat[,2], y = (three_FLD_veg_SL_Dmat_3_SIF - three_FLD_L_veg_SIF)), aes(x=x, y=y, colour = "3FLD"), 
             alpha = 0.75, size = 4.5, shape = 18) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_error_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text_error) +
  scale_colour_manual(values=colsRet, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_Stray_error_Dmat_3_SIFl",SIF_level,".pdf",sep=""),
       plot = Stray_error_Dmat_3, width = 10, height = 6, units = "cm", dpi = 300)


# Stray errors for sFLD and 3FLD: Very low stray light scenario
Stray_error_Dmat_4 <- ggplot() +
  geom_line(data = data.frame(x = window_mins_mat[,2], y = (sFLD_veg_SL_Dmat_4_SIF - sFLD_L_veg_SIF)), aes(x=x, y=y, colour = "sFLD"), 
             alpha = 0.75, size = 1) +
  geom_point(data = data.frame(x = window_mins_mat[,2], y = (sFLD_veg_SL_Dmat_4_SIF - sFLD_L_veg_SIF)), aes(x=x, y=y, colour = "sFLD"), 
             alpha = 0.75, size = 3.3, shape = 15) +
  geom_line(data = data.frame(x = window_mins_mat[,2], y = (three_FLD_veg_SL_Dmat_4_SIF - three_FLD_L_veg_SIF)), aes(x=x, y=y, colour = "3FLD"), 
             alpha = 0.75, size = 1) +
  geom_point(data = data.frame(x = window_mins_mat[,2], y = (three_FLD_veg_SL_Dmat_4_SIF - three_FLD_L_veg_SIF)), aes(x=x, y=y, colour = "3FLD"), 
             alpha = 0.75, size = 4.5, shape = 18) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_error_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text_error) +
  scale_colour_manual(values=colsRet, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) 

ggsave(filename = paste(path_figs, "Fig_Stray_error_Dmat_4_SIFl",SIF_level,".pdf",sep=""),
       plot = Stray_error_Dmat_4, width = 10, height = 6, units = "cm", dpi = 300)


### Example laser line to demonstrate what spectral stray light is

## make plots comparing SLC-corrected data and uncorrected data for a laser line
## (to demonstrate what spectral stray light is)


# choose laser lines to plot
# colnames(LineData)[25]
linesToPlot <- c(25)

# define minimum value for y-axis
yAxisMin <- 0.00001

# Normalize by maximum
LineData_norm <- LineData[ ,linesToPlot[1]]/max(LineData[ ,linesToPlot[1]])
LineData_SLC_norm <- LineData_SLC[ ,linesToPlot[1]]/max(LineData_SLC[ ,linesToPlot[1]])


# add small offset to all laser line data for plotting purposes 
Offset <- 0.00002
LineData_Offset <- LineData_norm + Offset
LineData_SLC_Offset <- LineData_SLC_norm + Offset 


# choose colors for laser lines
LL_cols <- data.frame(line=1:3,
                      raw = c("#99d8c9","grey","#a6bddb"),
                      slc = c("#238b45","black","#0868ac"))
LL_cols$raw <- as.character(LL_cols$raw); LL_cols$slc <- as.character(LL_cols$slc)

tiff(width=7, height=5, file=paste(path_figs, "Figure_Example_Laser_line_uncorr.tiff", sep = ""), units="in", res=600)
#par(mfrow=c(1,1), mar=c(3,3,1,1), oma=c(2,2,0,0))
# Plot for uncorrected
plot(x=xVals, y=LineData_Offset, type="l",col=adjustcolor(col=LL_cols$raw[1], alpha.f=0.8), 
     lwd=2,
     cex.axis=1.5,
     cex.lab = 1.5,
     ylim=c(yAxisMin,max(LineData_Offset)),
     xlab="Wavelength (nm)",
     ylab="Relative signal",
     log = 'y')
abline(h= c(0.00001, 0.0001, 0.001, 0.01, 0.1), col="grey80", lty=1)
dev.off()

tiff(width=7, height=5, file=paste(path_figs, "Figure_Example_Laser_line_corr.tiff", sep = ""), units="in", res=600)
#par(mfrow=c(1,1), mar=c(3,3,1,1), oma=c(2,2,0,0))
plot(x =xVals, y = LineData_SLC_Offset, type = "l", col=adjustcolor(col=LL_cols$slc[1],alpha.f=0.8), 
     lwd=2,
     cex.axis=1.5,
     cex.lab = 1.5,
     ylim=c(yAxisMin,max(LineData_Offset)),
     xlab="Wavelength (nm)",
     ylab="Relative signal",
     log = 'y')
abline(h= c(0.00001, 0.0001, 0.001, 0.01, 0.1), col="grey80", lty=1)
dev.off()


### Plots to show D matrix raster and one row (probably for supplemental info)
cols <- rev(magma(100))

tiff(width=7, height=5, units="in", res=600, "Figures/Figure_D_Mat.tiff")
#par(xpd = FALSE, mar = c(0,0,2,5)) # mai command might help reduce white space, e.g. mai = c(1, 0.4, 0.4, 0.4)
D_Mat_rast <- raster(D_Mat, 
                     xmn=1, xmx=2160,
                     ymn=1, ymx=2160)
plot(D_Mat_rast, col = cols, box=F)
dev.off()






#### EXPORT SIMULATIONS ############################################################################################
# Export the reference and simulated vegetation spectra for the four stray light scenarios.
# (These useful to share for trying other retrieval methods)

### Set working directory for outputs
#setwd('../SL_SIF_analysis_repo')


# Export spectra without stray light
# Note that I need to transpose SolarData_SLC_frames because apply functions used to define other cubes transposed them and
# I want all exported arrays to have the same dimensions. I checked that a simple transpositon was needed:
# test_apply <- apply(SolarData_SLC_frames,1,function(x) x+0)
# all.equal(test_apply, t(SolarData_SLC_frames)) #identical
write.table(SolarData_SLC, file = paste(path_spectra,"Solar_no_SL.csv", sep = ""), row.names = FALSE, col.names = c("Solar_radiance"), sep = ",")
write.ENVI(t(SolarData_SLC_frames), paste(path_spectra,"Solar_no_SL_frames", sep = ""))
write.table(L_veg_SIF, file = paste(path_spectra,"Vegetation_no_SL_SIFl",SIF_level,".csv", sep = ""), row.names = FALSE, col.names = c("Veg_radiance"), sep = ",")
write.ENVI(L_veg_SIF_frames, paste(path_spectra,"Vegetation_no_SL_SIFL",SIF_level, sep = ""))

# Export spectra with stray light simulated from D matrix
# solar
write.table(sol_SL_Dmat_1, file = paste(path_spectra, "Solar_SL_Dmat_scenario1.csv", sep = ""), row.names = FALSE, col.names = c("Solar_radiance"), sep = ",")
write.table(sol_SL_Dmat_2, file = paste(path_spectra, "Solar_SL_Dmat_scenario2.csv", sep = ""), row.names = FALSE, col.names = c("Solar_radiance"), sep = ",")
write.table(sol_SL_Dmat_3, file = paste(path_spectra, "Solar_SL_Dmat_scenario3.csv", sep = ""), row.names = FALSE, col.names = c("Solar_radiance"), sep = ",")
write.table(sol_SL_Dmat_4, file = paste(path_spectra, "Solar_SL_Dmat_scenario4.csv", sep = ""), row.names = FALSE, col.names = c("Solar_radiance"), sep = ",")
write.ENVI(sol_SL_Dmat_1_frames, paste(path_spectra, "Solar_SL_Dmat_scenario1_frames", sep = ""))
write.ENVI(sol_SL_Dmat_2_frames, paste(path_spectra, "Solar_SL_Dmat_scenario2_frames", sep = ""))
write.ENVI(sol_SL_Dmat_3_frames, paste(path_spectra, "Solar_SL_Dmat_scenario3_frames", sep = ""))
write.ENVI(sol_SL_Dmat_4_frames, paste(path_spectra, "Solar_SL_Dmat_scenario4_frames", sep = ""))
# vegetation
write.table(veg_SL_Dmat_1_SIF, file = paste(path_spectra, "Vegetation_SL_Dmat_scenario1_SIFl",SIF_level,".csv", sep = ""), row.names = FALSE, col.names = c("Veg_radiance"), sep = ",")
write.table(veg_SL_Dmat_2_SIF, file = paste(path_spectra, "Vegetation_SL_Dmat_scenario2_SIFl",SIF_level,".csv", sep = ""), row.names = FALSE, col.names = c("Veg_radiance"), sep = ",")
write.table(veg_SL_Dmat_3_SIF, file = paste(path_spectra, "Vegetation_SL_Dmat_scenario3_SIFl",SIF_level,".csv", sep = ""), row.names = FALSE, col.names = c("Veg_radiance"), sep = ",")
write.table(veg_SL_Dmat_4_SIF, file = paste(path_spectra, "Vegetation_SL_Dmat_scenario4_SIFl",SIF_level,".csv", sep = ""), row.names = FALSE, col.names = c("Veg_radiance"), sep = ",")
write.ENVI(veg_SL_Dmat_1_frames, paste(path_spectra, "Vegetation_SL_Dmat_scenario1_SIFl",SIF_level, sep = ""))
write.ENVI(veg_SL_Dmat_2_frames, paste(path_spectra, "Vegetation_SL_Dmat_scenario2_SIFl",SIF_level, sep = ""))
write.ENVI(veg_SL_Dmat_3_frames, paste(path_spectra, "Vegetation_SL_Dmat_scenario3_SIFl",SIF_level, sep = ""))
write.ENVI(veg_SL_Dmat_4_frames, paste(path_spectra, "Vegetation_SL_Dmat_scenario4_SIFl",SIF_level, sep = ""))

### EXPORT OBJECTS FOR USE ELSEWHERE ############################################################################################

# Export differences between with and w/o stray light for 3FLD
three_FLD_diffs <- cbind(window_mins_mat[ ,2],
                         (three_FLD_veg_SL_Dmat_1_SIF - three_FLD_L_veg_SIF),
                         (three_FLD_veg_SL_Dmat_2_SIF - three_FLD_L_veg_SIF),
                         (three_FLD_veg_SL_Dmat_3_SIF - three_FLD_L_veg_SIF),
                         (three_FLD_veg_SL_Dmat_4_SIF - three_FLD_L_veg_SIF))
three_FLD_diffs_df <- data.frame(three_FLD_diffs)
colnames(three_FLD_diffs_df) <- c("wavelength",
                                  paste("three_FLD_SL_Dmat_1_", "SIFl", SIF_level, "_diff", sep = ""),
                                  paste("three_FLD_SL_Dmat_2_", "SIFl", SIF_level, "_diff", sep = ""),
                                  paste("three_FLD_SL_Dmat_3_", "SIFl", SIF_level, "_diff", sep = ""),
                                  paste("three_FLD_SL_Dmat_4_", "SIFl", SIF_level, "_diff", sep = ""))
         
write.table(three_FLD_diffs_df, file = paste(path_spectra, "3FLD_", "SIFl", SIF_level, "_stray_no_stray_differences.csv", sep = ""), row.names = FALSE, sep = ",")


# Export Refl_99_interp for use in Matlab
write.table(Refl_99_interp, file = paste(path_spectra,"HL21-24_Spec_C1HL21_interpolated.csv", sep = ""), row.names = FALSE, col.names = c("Reflectance"), sep = ",")

# Export window_mins_mat for use in making table in SFMresults_summary.R and make_stray_light_error_figure.R
write.table(window_mins_mat, file = paste(path_spectra, "Estimated_window_mins_ind_wl.csv", sep = ""), row.names = FALSE, col.names = c("Index","Wavelength(nm)"), sep = ",")

# Export SIF spectra (with name customized to show the scaling)
write.table(SIF_2perc, file = paste(path_spectra, "SIF_SIFl", SIF_level, ".csv", sep = ""), row.names = FALSE, sep = ",")

