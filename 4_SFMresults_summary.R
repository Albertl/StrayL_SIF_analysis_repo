# Script for examining and organizing SFM retrieval results from Luis Alonso. SFM
# was applied to the stray light scenarios described in the main text of the paper.
# Script by KC Cushman and Loren Albert
#
# Contact: Loren Albert, lalbert@email.arizona.edu; loren.albert@mail.wvu.edu
# Here were the scenarios used for SFM retrieval:
#
# Scenarios for SFM:
# No spectral stray light scenario
# 1) Solar_no_SL.csv paired with Vegetation_no_SL_SIFl1.csv (SIF ~ 1mW @ 760 nm)
# 2) Solar_no_SL.csv paired with Vegetation_no_SL_SIFl2.csv (SIF ~ 2mW @ 760 nm)
# 3) Solar_no_SL.csv paired with Vegetation_no_SL_SIFl3.csv (SIF ~ 3mW @ 760 nm)
# 
# Low spectral stray light scenario
# 4) Solar_SL_Dmat_scenario4.csv paired with Vegetation_SL_Dmat_scenario4_SIFl1.csv (SIF ~ 1mW @ 760 nm)
# 5) Solar_SL_Dmat_scenario4.csv paired with Vegetation_SL_Dmat_scenario4_SIFl2.csv (SIF ~ 2mW @ 760 nm)
# 6) Solar_SL_Dmat_scenario4.csv paired with Vegetation_SL_Dmat_scenario4_SIFl3.csv (SIF ~ 3mW @ 760 nm)
# 
# Medium spectral stray light scenario
# 7) Solar_SL_Dmat_scenario3.csv paired with Vegetation_SL_Dmat_scenario3.csv_SIFl1.csv (SIF ~ 1mW @ 760 nm)
# 8) Solar_SL_Dmat_scenario3.csv paired with Vegetation_SL_Dmat_scenario3.csv_SIFl2.csv (SIF ~ 2mW @ 760 nm)
# 9) Solar_SL_Dmat_scenario3.csv paired with Vegetation_SL_Dmat_scenario3.csv_SIFl3.csv (SIF ~ 3mW @ 760 nm)
# 
# High spectral stray light scenario
# 10) Solar_SL_Dmat_scenario2.csv paired with Vegetation_SL_Dmat_scenario2.csv_SIFl1.csv (SIF ~ 1mW @ 760 nm)
# 11) Solar_SL_Dmat_scenario2.csv paired with Vegetation_SL_Dmat_scenario2.csv_SIFl2.csv (SIF ~ 2mW @ 760 nm)
# 12) Solar_SL_Dmat_scenario2.csv paired with Vegetation_SL_Dmat_scenario2.csv_SIFl3.csv (SIF ~ 3mW @ 760 nm)
# 
# Very High spectral stray light scenario
# 13) Solar_SL_Dmat_scenario1.csv paired with Vegetation_SL_Dmat_scenario1_SIFl1.csv (SIF ~ 1mW @ 760 nm)
# 14) Solar_SL_Dmat_scenario1.csv paired with Vegetation_SL_Dmat_scenario1_SIFl2.csv (SIF ~ 2mW @ 760 nm)
# 15) Solar_SL_Dmat_scenario1.csv paired with Vegetation_SL_Dmat_scenario1_SIFl3.csv (SIF ~ 3mW @ 760 nm)

# # Luis' instructions about the columns of the data (shows that they correspond to the scenario numbers above)
# scenario1 correspond to VERY HIGH SL (13-15),
# scenario2 correspond to HIGH SL (10-12),
# scenario3 correspond to MEDIUM SL (7-9),
# and scenario4 to LOW SL (4-6)
#
#
# Copyright (C) 2022  Loren Albert.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.



## Housekeeping--------------------------------------------------------------
## Remove workspace objects if necessary
# rm(list=ls())

require(ncdf4)
require(ggplot2)
require(viridis)
library(htmlTable)

## Change working directory if necessary
# setwd('/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/SL_SIF_analysis_repo/')


## define paths
path_figs <- "./Figures/"
path_spectra <- "./Outputs/"
path_data_LL <- "./Inputs/"

## Import fluoresence shape (F_true) for stray light error calculation and plots
SIF_true_SIFl3 <- read.csv(paste(path_spectra, "SIF_SIFl3.csv", sep = ""))
SIF_true_SIFl3 <- SIF_true_SIFl3[,1]
SIF_true_SIFl2 <- read.csv(paste(path_spectra, "SIF_SIFl2.csv", sep = ""))
SIF_true_SIFl2 <- SIF_true_SIFl2[,1]
SIF_true_SIFl1 <- read.csv(paste(path_spectra, "SIF_SIFl1.csv", sep = ""))
SIF_true_SIFl1 <- SIF_true_SIFl1[,1]

## Import Wavelength values for each spectral pixel
Wavelengths <- read.csv(paste(path_data_LL,"predicted_spectral_pixels_2018.csv",sep=""))

## Import table with index and wavelength of Fraunhofer line and O2A band minima (exported from make_stray_light_simulations_v2.R)
window_mins_mat <- read.csv(paste(path_spectra, "Estimated_window_mins_ind_wl.csv", sep = ""))
window_mins_ind <- window_mins_mat[,1]

## Import and open file from Luis
# Open using the "ncdf4" package:

#SFM_data <- ncdf4::nc_open("SFM_on_SL_data_using_3FLD_inicialization.nc")
SFM_data <- ncdf4::nc_open("./Inputs/SFM_on_SL-SIF_scenariosl.nc")

# You can enter the name of this list to see information about what's included:

SFM_data

# There are six variables, we can get each one by name:

Fluorescence_A <- ncdf4::ncvar_get(SFM_data, varid = "Fluorescence_A")
# Oxygen-A band fluorescence
# array of 588 (wavelengths) by 15 (stray light scenarios)

R_A <- ncdf4::ncvar_get(SFM_data, varid = "R_A")
# Oxygen-A band radiance?
# array of 588 (wavelengths) by 15 (stray light scenarios)

wvlA <- ncdf4::ncvar_get(SFM_data, varid = "wvlA")
# Oxygen-A band wavelengths
# array of 588 (wavelengths) by 15 (stray light scenarios)

Fluorescence_B <- ncdf4::ncvar_get(SFM_data, varid = "Fluorescence_B")
# Oxygen-B band fluorescence
# array of 314 (wavelengths) by 15 (stray light scenarios)

R_B <- ncdf4::ncvar_get(SFM_data, varid = "R_B")
# Oxygen-B band radiance
# array of 314 (wavelengths) by 15 (stray light scenarios)

wvlB <- ncdf4::ncvar_get(SFM_data, varid = "wvlB")
# Oxygen-B band wavelengths
# array of 314 (wavelengths) by 15 (stray light scenarios)




### Plots ---------------------------------------------------------------------

## Set working directory
#setwd(path_figs)

## choose x-values for entire instrument range
xVals <- Wavelengths[,2]

## choose y limits (good to coordinate with ylimits for 3FLD in make_stray_light_simulations_v2.R)
# sif_ylim <- c(-.002, max(SIF_true_SIFl3) + 0.0005)
sif_ylim <- c(-.0001, (max(SIF_true_SIFl2) + 0.0005)) # This range works well for SIF level 1 and 2

## Quick plot of all O2_A scenarios
viridis_cols <- viridis(15)
plot(wvlA[,15], Fluorescence_A[,15], type = "l")
for(i in 1:15){
  plot(wvlA[,i], Fluorescence_A[,i], col = viridis_cols[i], ylim = c(0,0.0045))
}

## Quick plot of all O2_B scenarios
plot(wvlB[,15], Fluorescence_B[,15], type = "l")
for(i in 1:15){
  lines(wvlB[,i], Fluorescence_B[,i], col = viridis_cols[i])
}

## Define color palette
colsRetStray <- c("SFM with stray"="#DAC228","SFM no stray"="#122B5D","Simulated SIF"="black")

## Define text for use in graphs
radiance_text = expression('Radiance (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')
radiance_text_log = expression('Log Radiance (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')
radiance_text_error = expression('Stray light error (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')
radiance_text_error = expression(atop('Stray light error', paste('(W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')))

## Define indices for 687.1nm and 760.5nm bands (within wavelength subsets defined by Luis)
O2B_band_ind <- which(round(wvlB[ ,1], digits = 2) == 687.1)
O2A_band_ind <- which(abs(wvlA[ ,1] - 760.5) == min(abs(wvlA[ ,1] - 760.5)))

# SFM: Very high stray light scenario (Low SIF (SIFl1))
# SFM: Very high stray light scenario (Low SIF (SIFl1))
SIFl1_ret_SFM_Dmat_1 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl1), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,1], y = Fluorescence_B[O2B_band_ind,1]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,13], y = Fluorescence_B[O2B_band_ind,13]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,1], y = Fluorescence_A[O2A_band_ind,1]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,13], y = Fluorescence_A[O2A_band_ind,13]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_1_SIFl1.pdf",sep = ""),
       plot = SIFl1_ret_SFM_Dmat_1, width = 10, height = 6, units = "cm", dpi = 300)

# SFM: Very high stray light scenario (Mid SIF (SIFl2))
SIFl2_ret_SFM_Dmat_1 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl2), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,2], y = Fluorescence_B[O2B_band_ind,2]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,14], y = Fluorescence_B[O2B_band_ind,14]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,2], y = Fluorescence_A[O2A_band_ind,2]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,14], y = Fluorescence_A[O2A_band_ind,14]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_1_SIFl2.pdf",sep = ""),
       plot = SIFl2_ret_SFM_Dmat_1, width = 10, height = 6, units = "cm", dpi = 300)

# SFM: Very high stray light scenario (High SIF (SIFl3))
SIFl3_ret_SFM_Dmat_1 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl3), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,3], y = Fluorescence_B[O2B_band_ind,3]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,15], y = Fluorescence_B[O2B_band_ind,15]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,3], y = Fluorescence_A[O2A_band_ind,3]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,15], y = Fluorescence_A[O2A_band_ind,15]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_1_SIFl3.pdf",sep = ""),
       plot = SIFl3_ret_SFM_Dmat_1, width = 10, height = 6, units = "cm", dpi = 300)



# SFM: High stray light scenario (Low SIF (SIFl1))
SIFl1_ret_SFM_Dmat_2 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl1), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,1], y = Fluorescence_B[O2B_band_ind,1]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,10], y = Fluorescence_B[O2B_band_ind,10]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,1], y = Fluorescence_A[O2A_band_ind,1]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,10], y = Fluorescence_A[O2A_band_ind,10]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_2_SIFl1.pdf",sep = ""),
       plot = SIFl1_ret_SFM_Dmat_2, width = 10, height = 6, units = "cm", dpi = 300)

# SFM: High stray light scenario (Mid SIF (SIFl2))
SIFl2_ret_SFM_Dmat_2 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl2), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,2], y = Fluorescence_B[O2B_band_ind,2]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,11], y = Fluorescence_B[O2B_band_ind,11]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,2], y = Fluorescence_A[O2A_band_ind,2]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,11], y = Fluorescence_A[O2A_band_ind,11]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_2_SIFl2.pdf",sep = ""),
       plot = SIFl2_ret_SFM_Dmat_2, width = 10, height = 6, units = "cm", dpi = 300)

# SFM: High stray light scenario (High SIF (SIFl3))
SIFl3_ret_SFM_Dmat_2 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl3), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,3], y = Fluorescence_B[O2B_band_ind,3]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,12], y = Fluorescence_B[O2B_band_ind,12]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,3], y = Fluorescence_A[O2A_band_ind,3]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,12], y = Fluorescence_A[O2A_band_ind,12]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_2_SIFl3.pdf",sep = ""),
       plot = SIFl3_ret_SFM_Dmat_2, width = 10, height = 6, units = "cm", dpi = 300)



# SFM: Medium stray light scenario (Low SIF (SIFl1))
SIFl1_ret_SFM_Dmat_3 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl1), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,1], y = Fluorescence_B[O2B_band_ind,1]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,7], y = Fluorescence_B[O2B_band_ind,7]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,1], y = Fluorescence_A[O2A_band_ind,1]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,7], y = Fluorescence_A[O2A_band_ind,7]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_3_SIFl1.pdf",sep = ""),
       plot = SIFl1_ret_SFM_Dmat_3, width = 10, height = 6, units = "cm", dpi = 300)

# SFM: Medium stray light scenario (Mid SIF (SIFl2))
SIFl2_ret_SFM_Dmat_3 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl2), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,2], y = Fluorescence_B[O2B_band_ind,2]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,8], y = Fluorescence_B[O2B_band_ind,8]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,2], y = Fluorescence_A[O2A_band_ind,2]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,8], y = Fluorescence_A[O2A_band_ind,8]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_3_SIFl2.pdf",sep = ""),
       plot = SIFl2_ret_SFM_Dmat_3, width = 10, height = 6, units = "cm", dpi = 300)

# SFM: Medium stray light scenario (High SIF (SIFl3))
SIFl3_ret_SFM_Dmat_3 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl3), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,3], y = Fluorescence_B[O2B_band_ind,3]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,9], y = Fluorescence_B[O2B_band_ind,9]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,3], y = Fluorescence_A[O2A_band_ind,3]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,9], y = Fluorescence_A[O2A_band_ind,9]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_3_SIFl3.pdf",sep = ""),
       plot = SIFl3_ret_SFM_Dmat_3, width = 10, height = 6, units = "cm", dpi = 300)




# SFM: Low stray light scenario (Mid SIF (SIFl1))
SIFl1_ret_SFM_Dmat_4 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl1), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,1], y = Fluorescence_B[O2B_band_ind,1]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,4], y = Fluorescence_B[O2B_band_ind,4]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,1], y = Fluorescence_A[O2A_band_ind,1]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,4], y = Fluorescence_A[O2A_band_ind,4]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_4_SIFl1.pdf",sep = ""),
       plot = SIFl1_ret_SFM_Dmat_4, width = 10, height = 6, units = "cm", dpi = 300)

# SFM: Low stray light scenario (Mid SIF (SIFl2))
SIFl2_ret_SFM_Dmat_4 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl2), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,2], y = Fluorescence_B[O2B_band_ind,2]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,5], y = Fluorescence_B[O2B_band_ind,5]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,2], y = Fluorescence_A[O2A_band_ind,2]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,5], y = Fluorescence_A[O2A_band_ind,5]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_4_SIFl2.pdf",sep = ""),
       plot = SIFl2_ret_SFM_Dmat_4, width = 10, height = 6, units = "cm", dpi = 300)

# SFM: Low stray light scenario (High SIF (SIFl3))
SIFl3_ret_SFM_Dmat_4 <- ggplot() +
  # SIF spectrum
  geom_line(data = data.frame(x = xVals, y = SIF_true_SIFl3), aes(x=x, y=y, colour = "Simulated SIF"), alpha = 0.45, size = 1) +
  # Retrieved SIF without stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,3], y = Fluorescence_B[O2B_band_ind,3]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2B region
  geom_point(data = data.frame(x = wvlB[O2B_band_ind,6], y = Fluorescence_B[O2B_band_ind,6]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  # Retrieved SIF without stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,3], y = Fluorescence_A[O2A_band_ind,3]), aes(x=x, y=y, colour = "SFM no stray"), 
             alpha = 0.75, size = 3.3) +
  # Retrieved SIF with stray light in O2A region
  geom_point(data = data.frame(x = wvlA[O2A_band_ind,6], y = Fluorescence_A[O2A_band_ind,6]), aes(x=x, y=y, colour = "SFM with stray"), 
             alpha = 0.70, size = 3.3, shape = 17) +
  theme_classic() +
  scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
  scale_y_continuous(limits = sif_ylim) +
  xlab("Wavelength (nm)") +
  ylab(radiance_text) +
  scale_colour_manual(values=colsRetStray, name = "") +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

ggsave(filename = paste(path_figs,"Fig_SFM_SIF_ret_Dmat_4_SIFl3.pdf",sep = ""),
       plot = SIFl3_ret_SFM_Dmat_4, width = 10, height = 6, units = "cm", dpi = 300)



### Calculate differences between scenarios with and w/o stray light ---------------------------

# Calculate difference between corresponding scenarios
# Oxygen B region fit
SFM_SL_Dmat_1_SIFl1_diff_B <- (Fluorescence_B[,13] - Fluorescence_B[,1])
SFM_SL_Dmat_1_SIFl2_diff_B <- (Fluorescence_B[,14] - Fluorescence_B[,2])
SFM_SL_Dmat_1_SIFl3_diff_B <- (Fluorescence_B[,15] - Fluorescence_B[,3])

SFM_SL_Dmat_2_SIFl1_diff_B <- (Fluorescence_B[,10] - Fluorescence_B[,1])
SFM_SL_Dmat_2_SIFl2_diff_B <- (Fluorescence_B[,11] - Fluorescence_B[,2])
SFM_SL_Dmat_2_SIFl3_diff_B <- (Fluorescence_B[,12] - Fluorescence_B[,3])

SFM_SL_Dmat_3_SIFl1_diff_B <- (Fluorescence_B[,7] - Fluorescence_B[,1])
SFM_SL_Dmat_3_SIFl2_diff_B <- (Fluorescence_B[,8] - Fluorescence_B[,2])
SFM_SL_Dmat_3_SIFl3_diff_B <- (Fluorescence_B[,9] - Fluorescence_B[,3])

SFM_SL_Dmat_4_SIFl1_diff_B <- (Fluorescence_B[,4] - Fluorescence_B[,1])
SFM_SL_Dmat_4_SIFl2_diff_B <- (Fluorescence_B[,5] - Fluorescence_B[,2])
SFM_SL_Dmat_4_SIFl3_diff_B <- (Fluorescence_B[,6] - Fluorescence_B[,3])

# Oxygen A region fit
SFM_SL_Dmat_1_SIFl1_diff_A <- (Fluorescence_A[,13] - Fluorescence_A[,1])
SFM_SL_Dmat_1_SIFl2_diff_A <- (Fluorescence_A[,14] - Fluorescence_A[,2])
SFM_SL_Dmat_1_SIFl3_diff_A <- (Fluorescence_A[,15] - Fluorescence_A[,3])

SFM_SL_Dmat_2_SIFl1_diff_A <- (Fluorescence_A[,10] - Fluorescence_A[,1])
SFM_SL_Dmat_2_SIFl2_diff_A <- (Fluorescence_A[,11] - Fluorescence_A[,2])
SFM_SL_Dmat_2_SIFl3_diff_A <- (Fluorescence_A[,12] - Fluorescence_A[,3])

SFM_SL_Dmat_3_SIFl1_diff_A <- (Fluorescence_A[,7] - Fluorescence_A[,1])
SFM_SL_Dmat_3_SIFl2_diff_A <- (Fluorescence_A[,8] - Fluorescence_A[,2])
SFM_SL_Dmat_3_SIFl3_diff_A <- (Fluorescence_A[,9] - Fluorescence_A[,3])

SFM_SL_Dmat_4_SIFl1_diff_A <- (Fluorescence_A[,4] - Fluorescence_A[,1])
SFM_SL_Dmat_4_SIFl2_diff_A <- (Fluorescence_A[,5] - Fluorescence_A[,2])
SFM_SL_Dmat_4_SIFl3_diff_A <- (Fluorescence_A[,6] - Fluorescence_A[,3])

# Make objects by concatenating O2-A and O2-B regions and filling in NAs for full 2160 length
Filler1 <- as.numeric(rep(NA,276))
Filler2 <- as.numeric(rep(NA,977))
Filler3 <- as.numeric(rep(NA,5))
SFM_SL_Dmat_1_SIFl1_diff <- c(Filler1, SFM_SL_Dmat_1_SIFl1_diff_B, Filler2,
                              SFM_SL_Dmat_1_SIFl1_diff_A, Filler3)
SFM_SL_Dmat_1_SIFl2_diff <- c(Filler1, SFM_SL_Dmat_1_SIFl2_diff_B, Filler2,
                              SFM_SL_Dmat_1_SIFl2_diff_A, Filler3)
SFM_SL_Dmat_1_SIFl3_diff <- c(Filler1, SFM_SL_Dmat_1_SIFl3_diff_B, Filler2,
                              SFM_SL_Dmat_1_SIFl3_diff_A, Filler3)
SFM_SL_Dmat_2_SIFl1_diff <- c(Filler1, SFM_SL_Dmat_2_SIFl1_diff_B, Filler2,
                              SFM_SL_Dmat_2_SIFl1_diff_A, Filler3)
SFM_SL_Dmat_2_SIFl2_diff <- c(Filler1, SFM_SL_Dmat_2_SIFl2_diff_B, Filler2,
                              SFM_SL_Dmat_2_SIFl2_diff_A, Filler3)
SFM_SL_Dmat_2_SIFl3_diff <- c(Filler1, SFM_SL_Dmat_2_SIFl3_diff_B, Filler2,
                              SFM_SL_Dmat_2_SIFl3_diff_A, Filler3)
SFM_SL_Dmat_3_SIFl1_diff <- c(Filler1, SFM_SL_Dmat_3_SIFl1_diff_B, Filler2,
                              SFM_SL_Dmat_3_SIFl1_diff_A, Filler3)
SFM_SL_Dmat_3_SIFl2_diff <- c(Filler1, SFM_SL_Dmat_3_SIFl2_diff_B, Filler2,
                              SFM_SL_Dmat_3_SIFl2_diff_A, Filler3)
SFM_SL_Dmat_3_SIFl3_diff <- c(Filler1, SFM_SL_Dmat_3_SIFl3_diff_B, Filler2,
                              SFM_SL_Dmat_3_SIFl3_diff_A, Filler3)
SFM_SL_Dmat_4_SIFl1_diff <- c(Filler1, SFM_SL_Dmat_4_SIFl1_diff_B, Filler2,
                              SFM_SL_Dmat_4_SIFl1_diff_A, Filler3)
SFM_SL_Dmat_4_SIFl2_diff <- c(Filler1, SFM_SL_Dmat_4_SIFl2_diff_B, Filler2,
                              SFM_SL_Dmat_4_SIFl2_diff_A, Filler3)
SFM_SL_Dmat_4_SIFl3_diff <- c(Filler1, SFM_SL_Dmat_4_SIFl3_diff_B, Filler2,
                              SFM_SL_Dmat_4_SIFl3_diff_A, Filler3)

# Combine the differences between scenarios with and without stray light into a dataframe
SFM_diffs <- data.frame(cbind(Wavelengths[ ,2],
                              SFM_SL_Dmat_1_SIFl1_diff,
                              SFM_SL_Dmat_1_SIFl2_diff,
                              SFM_SL_Dmat_1_SIFl3_diff,
                              SFM_SL_Dmat_2_SIFl1_diff,
                              SFM_SL_Dmat_2_SIFl2_diff,
                              SFM_SL_Dmat_2_SIFl3_diff,
                              SFM_SL_Dmat_3_SIFl1_diff,
                              SFM_SL_Dmat_3_SIFl2_diff,
                              SFM_SL_Dmat_3_SIFl3_diff,
                              SFM_SL_Dmat_4_SIFl1_diff,
                              SFM_SL_Dmat_4_SIFl2_diff,
                              SFM_SL_Dmat_4_SIFl3_diff))



### Plot differences between scenarios with and w/o stray light ------------------------------------

# All of the oxygen B band retrievals show a positive difference between scenarios with and w/o stray light
min(SFM_diffs[1:600, ], na.rm = T)

# Plots to explore negative difference in oxygen A bands for some scenarios

compare_O2A_scenarios <- function(w_stray, no_stray){
  ggplot() +
    geom_point(data = data.frame(x = wvlA[ ,1],
                                 y = w_stray), aes(x = x, y = y, colour = "SFM with stray"), alpha = 0.8, size = 0.6) +
    geom_point(data = data.frame(x = wvlA[ ,1],
                                 y = no_stray), aes(x = x, y = y, colour = "SFM no stray"), alpha = 0.8, size = 0.6) +
    xlab("Wavelength (nm)") +
    ylab("Radiance") +
    scale_colour_manual(values=colsRetStray, name = "") +
    theme_bw()
}
  
stray_no_stray_plot_Dmat_1_SIFl1 <- compare_O2A_scenarios(Fluorescence_A[,13], Fluorescence_A[,1])
stray_no_stray_plot_Dmat_1_SIFl2 <- compare_O2A_scenarios(Fluorescence_A[,14], Fluorescence_A[,2])
stray_no_stray_plot_Dmat_1_SIFl3 <- compare_O2A_scenarios(Fluorescence_A[,15], Fluorescence_A[,3])

stray_no_stray_plot_Dmat_2_SIFl1 <- compare_O2A_scenarios(Fluorescence_A[,10], Fluorescence_A[,1])
stray_no_stray_plot_Dmat_2_SIFl2 <- compare_O2A_scenarios(Fluorescence_A[,11], Fluorescence_A[,2])
stray_no_stray_plot_Dmat_2_SIFl3 <- compare_O2A_scenarios(Fluorescence_A[,12], Fluorescence_A[,3])

stray_no_stray_plot_Dmat_3_SIFl1 <- compare_O2A_scenarios(Fluorescence_A[,7], Fluorescence_A[,1])
stray_no_stray_plot_Dmat_3_SIFl2 <- compare_O2A_scenarios(Fluorescence_A[,8], Fluorescence_A[,2])
stray_no_stray_plot_Dmat_3_SIFl3 <- compare_O2A_scenarios(Fluorescence_A[,9], Fluorescence_A[,3])

stray_no_stray_plot_Dmat_4_SIFl1 <- compare_O2A_scenarios(Fluorescence_A[,4], Fluorescence_A[,1])
stray_no_stray_plot_Dmat_4_SIFl2 <- compare_O2A_scenarios(Fluorescence_A[,5], Fluorescence_A[,2])
stray_no_stray_plot_Dmat_4_SIFl3 <- compare_O2A_scenarios(Fluorescence_A[,6], Fluorescence_A[,3])



### Calculate error as a proportion of SIF --------------------------------------------------------------

# Get indices corresponding with window minima in Fluorescence_B and Fluorescence_A objects
# Looks like O2-B wavelengths start at Wavelengths$spec_pixel_pred_no710[277]
SIF_B_ind <- Wavelengths[277:590, ]
SIF_A_ind <- Wavelengths[1568:2155, ]

# Calculate error as a proportion of SIF_true (called 'SIF_2perc' in make_stray_light_simulations_v2.R)
# Oxygen B region fit
SFM_SL_Dmat_1_SIFl1_error_B <- (Fluorescence_B[,13] - Fluorescence_B[,1])/SIF_true_SIFl1[SIF_B_ind$spec_pixel_list]
SFM_SL_Dmat_1_SIFl2_error_B <- (Fluorescence_B[,14] - Fluorescence_B[,2])/SIF_true_SIFl2[SIF_B_ind$spec_pixel_list]
SFM_SL_Dmat_1_SIFl3_error_B <- (Fluorescence_B[,15] - Fluorescence_B[,3])/SIF_true_SIFl3[SIF_B_ind$spec_pixel_list]

SFM_SL_Dmat_2_SIFl1_error_B <- (Fluorescence_B[,10] - Fluorescence_B[,1])/SIF_true_SIFl1[SIF_B_ind$spec_pixel_list]
SFM_SL_Dmat_2_SIFl2_error_B <- (Fluorescence_B[,11] - Fluorescence_B[,2])/SIF_true_SIFl2[SIF_B_ind$spec_pixel_list]
SFM_SL_Dmat_2_SIFl3_error_B <- (Fluorescence_B[,12] - Fluorescence_B[,3])/SIF_true_SIFl3[SIF_B_ind$spec_pixel_list]

SFM_SL_Dmat_3_SIFl1_error_B <- (Fluorescence_B[,7] - Fluorescence_B[,1])/SIF_true_SIFl1[SIF_B_ind$spec_pixel_list]
SFM_SL_Dmat_3_SIFl2_error_B <- (Fluorescence_B[,8] - Fluorescence_B[,2])/SIF_true_SIFl2[SIF_B_ind$spec_pixel_list]
SFM_SL_Dmat_3_SIFl3_error_B <- (Fluorescence_B[,9] - Fluorescence_B[,3])/SIF_true_SIFl3[SIF_B_ind$spec_pixel_list]

SFM_SL_Dmat_4_SIFl1_error_B <- (Fluorescence_B[,4] - Fluorescence_B[,1])/SIF_true_SIFl1[SIF_B_ind$spec_pixel_list]
SFM_SL_Dmat_4_SIFl2_error_B <- (Fluorescence_B[,5] - Fluorescence_B[,2])/SIF_true_SIFl2[SIF_B_ind$spec_pixel_list]
SFM_SL_Dmat_4_SIFl3_error_B <- (Fluorescence_B[,6] - Fluorescence_B[,3])/SIF_true_SIFl3[SIF_B_ind$spec_pixel_list]

# Oxygen A region fit
SFM_SL_Dmat_1_SIFl1_error_A <- (Fluorescence_A[,13] - Fluorescence_A[,1])/SIF_true_SIFl1[SIF_A_ind$spec_pixel_list]
SFM_SL_Dmat_1_SIFl2_error_A <- (Fluorescence_A[,14] - Fluorescence_A[,2])/SIF_true_SIFl2[SIF_A_ind$spec_pixel_list]
SFM_SL_Dmat_1_SIFl3_error_A <- (Fluorescence_A[,15] - Fluorescence_A[,3])/SIF_true_SIFl3[SIF_A_ind$spec_pixel_list]

SFM_SL_Dmat_2_SIFl1_error_A <- (Fluorescence_A[,10] - Fluorescence_A[,1])/SIF_true_SIFl1[SIF_A_ind$spec_pixel_list]
SFM_SL_Dmat_2_SIFl2_error_A <- (Fluorescence_A[,11] - Fluorescence_A[,2])/SIF_true_SIFl2[SIF_A_ind$spec_pixel_list]
SFM_SL_Dmat_2_SIFl3_error_A <- (Fluorescence_A[,12] - Fluorescence_A[,3])/SIF_true_SIFl3[SIF_A_ind$spec_pixel_list]

SFM_SL_Dmat_3_SIFl1_error_A <- (Fluorescence_A[,7] - Fluorescence_A[,1])/SIF_true_SIFl1[SIF_A_ind$spec_pixel_list]
SFM_SL_Dmat_3_SIFl2_error_A <- (Fluorescence_A[,8] - Fluorescence_A[,2])/SIF_true_SIFl2[SIF_A_ind$spec_pixel_list]
SFM_SL_Dmat_3_SIFl3_error_A <- (Fluorescence_A[,9] - Fluorescence_A[,3])/SIF_true_SIFl3[SIF_A_ind$spec_pixel_list]

SFM_SL_Dmat_4_SIFl1_error_A <- (Fluorescence_A[,4] - Fluorescence_A[,1])/SIF_true_SIFl1[SIF_A_ind$spec_pixel_list]
SFM_SL_Dmat_4_SIFl2_error_A <- (Fluorescence_A[,5] - Fluorescence_A[,2])/SIF_true_SIFl2[SIF_A_ind$spec_pixel_list]
SFM_SL_Dmat_4_SIFl3_error_A <- (Fluorescence_A[,6] - Fluorescence_A[,3])/SIF_true_SIFl3[SIF_A_ind$spec_pixel_list]

# Make error objects by concatenating O2-A and O2-B regions and filling in NAs for full 2160 length
Filler1 <- as.numeric(rep(NA,276))
Filler2 <- as.numeric(rep(NA,977))
Filler3 <- as.numeric(rep(NA,5))
SFM_SL_Dmat_1_SIFl1_error <- c(Filler1, SFM_SL_Dmat_1_SIFl1_error_B, Filler2,
                               SFM_SL_Dmat_1_SIFl1_error_A, Filler3)
SFM_SL_Dmat_1_SIFl2_error <- c(Filler1, SFM_SL_Dmat_1_SIFl2_error_B, Filler2,
                               SFM_SL_Dmat_1_SIFl2_error_A, Filler3)
SFM_SL_Dmat_1_SIFl3_error <- c(Filler1, SFM_SL_Dmat_1_SIFl3_error_B, Filler2,
                               SFM_SL_Dmat_1_SIFl3_error_A, Filler3)
SFM_SL_Dmat_2_SIFl1_error <- c(Filler1, SFM_SL_Dmat_2_SIFl1_error_B, Filler2,
                               SFM_SL_Dmat_2_SIFl1_error_A, Filler3)
SFM_SL_Dmat_2_SIFl2_error <- c(Filler1, SFM_SL_Dmat_2_SIFl2_error_B, Filler2,
                               SFM_SL_Dmat_2_SIFl2_error_A, Filler3)
SFM_SL_Dmat_2_SIFl3_error <- c(Filler1, SFM_SL_Dmat_2_SIFl3_error_B, Filler2,
                               SFM_SL_Dmat_2_SIFl3_error_A, Filler3)
SFM_SL_Dmat_3_SIFl1_error <- c(Filler1, SFM_SL_Dmat_3_SIFl1_error_B, Filler2,
                               SFM_SL_Dmat_3_SIFl1_error_A, Filler3)
SFM_SL_Dmat_3_SIFl2_error <- c(Filler1, SFM_SL_Dmat_3_SIFl2_error_B, Filler2,
                               SFM_SL_Dmat_3_SIFl2_error_A, Filler3)
SFM_SL_Dmat_3_SIFl3_error <- c(Filler1, SFM_SL_Dmat_3_SIFl3_error_B, Filler2,
                               SFM_SL_Dmat_3_SIFl3_error_A, Filler3)
SFM_SL_Dmat_4_SIFl1_error <- c(Filler1, SFM_SL_Dmat_4_SIFl1_error_B, Filler2,
                               SFM_SL_Dmat_4_SIFl1_error_A, Filler3)
SFM_SL_Dmat_4_SIFl2_error <- c(Filler1, SFM_SL_Dmat_4_SIFl2_error_B, Filler2,
                               SFM_SL_Dmat_4_SIFl2_error_A, Filler3)
SFM_SL_Dmat_4_SIFl3_error <- c(Filler1, SFM_SL_Dmat_4_SIFl3_error_B, Filler2,
                               SFM_SL_Dmat_4_SIFl3_error_A, Filler3)

### Make tables --------------------------------------------------------------
# Combine the stray light errors as a proportion of SIF into a dataframe
SFM_error_prop_to_SIF <- data.frame(cbind(Wavelengths[ ,2],
                                          SFM_SL_Dmat_1_SIFl1_error,
                                          SFM_SL_Dmat_1_SIFl2_error,
                                          SFM_SL_Dmat_1_SIFl3_error,
                                          SFM_SL_Dmat_2_SIFl1_error,
                                          SFM_SL_Dmat_2_SIFl2_error,
                                          SFM_SL_Dmat_2_SIFl3_error,
                                          SFM_SL_Dmat_3_SIFl1_error,
                                          SFM_SL_Dmat_3_SIFl2_error,
                                          SFM_SL_Dmat_3_SIFl3_error,
                                          SFM_SL_Dmat_4_SIFl1_error,
                                          SFM_SL_Dmat_4_SIFl2_error,
                                          SFM_SL_Dmat_4_SIFl3_error))

## Define indices for 687.1nm and 760.5nm bands within full 2160 bands
O2B_band_ind_full <- which(round(Wavelengths[ ,2], digits = 2) == 687.1)
O2A_band_ind_full <- which(abs(Wavelengths[ ,2] - 760.5) == min(abs(Wavelengths[ ,2] - 760.5)))

# Define row indices for SIFl2
SIFl2_ind <- c(1, 3, 6, 9, 12)
# Define row indices for SIFl1 (for supplement)
SIFl1_ind <- c(1, 2, 5, 8, 11)
# Define row indices for SIFl3 (for supplement)
SIFl3_ind <- c(1, 4, 7, 10, 13)

# Take the SIFl2 and index by the O2A and O2B bands to make results table for main text
SFM_stray_summary_table <- rbind(SFM_error_prop_to_SIF[O2B_band_ind_full, SIFl2_ind], 
                                 SFM_error_prop_to_SIF[O2A_band_ind_full, SIFl2_ind])
# Take the SIFl1, SIFl2, and SIFl3 and index by the O2A and O2B bands to make larger results table for supplement
SFM_stray_summary_table_supplement_SIFl1 <- rbind(SFM_error_prop_to_SIF[O2B_band_ind_full, SIFl1_ind], 
                                            SFM_error_prop_to_SIF[O2A_band_ind_full, SIFl1_ind])
SFM_stray_summary_table_supplement_SIFl3 <- rbind(SFM_error_prop_to_SIF[O2B_band_ind_full, SIFl3_ind], 
                                            SFM_error_prop_to_SIF[O2A_band_ind_full, SIFl3_ind])

# Convert from proportion to percentage for the tables
SFM_stray_summary_table_perc <- SFM_stray_summary_table * 100
SFM_stray_summary_table_supplement_SIFl1_perc <- SFM_stray_summary_table_supplement_SIFl1 * 100
SFM_stray_summary_table_supplement_SIFl3_perc <- SFM_stray_summary_table_supplement_SIFl3 * 100

# Transpose and organize
SFM_stray_summary_table_perc_cols <- t(SFM_stray_summary_table_perc[ , 2:ncol(SFM_stray_summary_table_perc)])
SFM_stray_summary_table_supplement_SIFl1_perc_cols <- t(SFM_stray_summary_table_supplement_SIFl1_perc[ , 2:ncol(SFM_stray_summary_table_supplement_SIFl1_perc)])
SFM_stray_summary_table_supplement_SIFl3_perc_cols <- t(SFM_stray_summary_table_supplement_SIFl3_perc[ , 2:ncol(SFM_stray_summary_table_supplement_SIFl3_perc)])
SFM_stray_summary_table_supplement_perc_cols <- rbind(SFM_stray_summary_table_supplement_SIFl1_perc_cols,
                                                      SFM_stray_summary_table_perc_cols,
                                                      SFM_stray_summary_table_supplement_SIFl3_perc_cols)

# Make tables for results
# Main text
htmlTable(round(SFM_stray_summary_table_perc_cols, digits = 2), 
          header =  as.character(round(SFM_stray_summary_table[ ,1], digits = 2)),
          rnames = paste("Scenario ", as.character(rep(c(1,2,3,4), times = 2)), "orders of magnitude"))
# Supplementary table with all three SIF levels
htmlTable(round(SFM_stray_summary_table_supplement_perc_cols, digits = 2),
          header =  as.character(round(SFM_stray_summary_table[ ,1], digits = 2)),
          rnames = paste("Scenario ", as.character(rep(c(1,2,3,4), times = 3)), "orders of magnitude"))


 ### Make plots to explore behavior of the errors proportional to SIF ---------------------------------------
plot(SFM_SL_Dmat_1_SIFl2_error, SFM_SL_Dmat_2_SIFl2_error, type = 'l')
plot(SFM_SL_Dmat_1_SIFl2_error, SFM_SL_Dmat_2_SIFl2_error, type = 'l', ylim = c(0,0.0075), xlim = c(0,0.0075))
plot(SFM_SL_Dmat_2_SIFl2_error, SFM_SL_Dmat_3_SIFl2_error, type = 'l')
plot(SFM_SL_Dmat_2_SIFl2_error, SFM_SL_Dmat_3_SIFl2_error, type = 'l', ylim = c(0,0.0075), xlim = c(0,0.0075))
plot(SFM_SL_Dmat_3_SIFl2_error, SFM_SL_Dmat_4_SIFl2_error, type = 'l')
plot(SFM_SL_Dmat_3_SIFl2_error, SFM_SL_Dmat_4_SIFl2_error, type = 'l', ylim = c(0,0.0075), xlim = c(0,0.0075))
all_SFM_errors <- ggplot() +
  geom_point(data = data.frame(x = SFM_SL_Dmat_1_SIFl2_error,
                              y = SFM_SL_Dmat_2_SIFl2_error), aes(x=x, y=y), alpha = 0.8, size = 0.6) +
  #scale_x_continuous(expand = c(0, 0)) +
  xlab("Stray light error as a proportion of SIF for highest stray light scenario") +
  ylab("Stray light error as a proportion of SIF for second highest stray light scenario") +
  theme_bw()
SFM_errors_zoom <- ggplot() +
  geom_point(data = data.frame(x = SFM_SL_Dmat_1_SIFl2_error,
                                y = SFM_SL_Dmat_2_SIFl2_error), aes(x=x, y=y), alpha = 0.8, size = 0.6) +
  scale_x_continuous(limits = c(-0.1, 0.0075)) +
  xlab("Stray light error as a proportion of SIF for highest stray light scenario") +
  ylab("Stray light error as a proportion of SIF for second highest stray light scenario") +
  theme_bw()

hist(SFM_SL_Dmat_1_SIFl1_error_B)
hist(SFM_SL_Dmat_1_SIFl1_error_A)




### Export differences and errors as proportion of SIF (for synthesis with other retrieval methods in a new plot) --------------------------------------------------

write.table(SFM_diffs, file = paste(path_spectra, "SFM_stray_no_stray_differences.csv", sep = ""), row.names = FALSE, sep = ",")

write.table(SFM_error_prop_to_SIF, file = paste(path_spectra, "SFM_errors_prop_to_SIF.csv", sep = ""), row.names = FALSE, sep = ",")

