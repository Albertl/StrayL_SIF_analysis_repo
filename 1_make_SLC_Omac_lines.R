# R script to make plots showing the effect of performing the stray light correction to laser lines, and
# to vegetation and spectralon (reference target).  The vegetation and spectralon spectra are output for use in other scripts. 
# (Based on make_SLC_plots.R script by Loren Albert and KC Cushman). Produces figure for stray light correction
# of laser lines, and radiometrically calibrated vegetation / spectralon image with and without spectral
# stray light correction.
#
# Contact: Loren Albert, lalbert@email.arizona.edu; loren.albert@mail.wvu.edu
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


######## Housekeeping #######################
## Remove workspace objects if necessary
# rm(list=ls())

## Change working directory if necessary
# setwd('/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/SL_SIF_analysis_repo/')


##### INPUTS ##################################
## load library to read SIF spectrometer image
library(caTools)
library(raster)


## Code to define directories 
  ## define path for inputs
  path_data <- './Inputs/'
  
  ## define path to figure outputs
  path_figs <- "./Figures/"
  
  ## define path to other outputs
  path_output <- "./Outputs"


  ## read in data
  ## C matrix (used to multiply raw data for stray light correction)
  Cmatrix <- read.csv(paste(path_data,"SLC_mat.csv",sep=""),header=F)
  Cmatrix <- as.matrix(Cmatrix)
  ## Wavelength values for each spectral pixel
  Wavelengths <- read.csv(paste(path_data,"predicted_spectral_pixels_2018.csv",sep=""))
  ## Laser line data
  LineData <- read.csv(paste(path_data,"Bracket_10ms_1000ms_scaled_by_ms_with_correction.csv",sep=""))

  ## Test data from OMAC roof 
  ## SIF spectrometer
  OmacData <-  read.ENVI(paste(path_data,"omac_BoP_point_stare_2018_04_23_16_48_29/raw_0_mean_dark_subtracted",sep=""))
  ## SIF spectrometer (final spectralon only-- no vegetation. May not ultimately be needed)
  OmacDataSpec <-  read.ENVI(paste(path_data,"omac_spectralon_only_point_stare_2018_04_23_17_38_17/raw_0_mean_dark_subtracted",sep=""))
  ## ASD FluoWatt data
  OmacData_ASD_SIF <- read.csv(paste(path_data,"OMAC_FluoWatt_SIF.csv",sep=""))
  OmacData_ASD_Leaf <- read.csv(paste(path_data,"OMAC_FluoWatt_Leaf.csv",sep=""))
  OmacData_ASD_Solar <- read.csv(paste(path_data,"OMAC_FluoWatt_Solar.csv",sep=""))
  
  ## Spatial pixel info:
  ## 1. Vegetation: 600-650
  ## 2. Spectralon
  ## 1st least reflective: 860-880
  ## 2nd least reflective: 910-930
  ## 3rd least reflective: 960-980
  ## 4th least (most) reflective: 1015-1035
  
  
  ## radiance calibration outputs: Coefficients from segmented lm (spectral x spatial x coefficients) with NAs interpolated
  Rad100_coef2 <- read.ENVI(paste(path_data,"Radiance_calib_figs/RC_DN_100_coefficients2_interp",sep=''))
  
  ## source script with all general functions for this instrument (including function for stray light correction to data)
  source("Stray_light_analysis_functions.R")
  

        

##### CALCULATIONS ###############################

  # STRAY LIGHT          
    # apply stray light correction to laser line data
        LineData_SLC <- SLC_fun_colwise(data_mat = LineData, C_mat = Cmatrix) # takes a few seconds to run
    # # apply stray light correction to validation source data
    #     ValiData_SLC <- SLC_fun_colwise(data_mat = ValiData, C_mat = Cmatrix)
    # apply stray light correction to OMAC solar data
        OmacData_SLC <- SLC_fun_colwise(data_mat = t(OmacData), C_mat = Cmatrix)
        OmacDataSpec_SLC <- SLC_fun_colwise(data_mat = t(OmacDataSpec), C_mat = Cmatrix)
        
        
  # RADIANCE CALIBRATION
    # first, divide DN by integration time
        OmacData_SLC_DNMS <- OmacData_SLC/100
        OmacData_DNMS <- t(OmacData)/100
        
        OmacDataSpec_SLC_DNMS <- OmacDataSpec_SLC/100
        OmacDataSpec_DNMS <- t(OmacDataSpec)/100
        
        
        # For segmented lm: First four coefficients are intercept, and three slope variables. Fifth and sixth coefficients are the breakpoints, (and 7th is r.squared and 8th is sigma)
        coef2incpt <- Rad100_coef2[,,1]
        coef2slp1 <- Rad100_coef2[,,2]
        coef2slp2 <- Rad100_coef2[,,3]
        coef2slp3 <- Rad100_coef2[,,4]
        coef2break1 <- Rad100_coef2[,,5]
        coef2break2 <- Rad100_coef2[,,6]
        
        # Apply segmented lm radiance calibration to OMAC data
        # Use dummy matrix fn with segmented lm coefficients from radiance calibration
        nameU1 <- c("U1.x", "U2.x")
        nameV1 <- c("psi1.x", "psi2.x")
        
        OmacData_SLC_Rad <- matrix(NA, 2160, 1600)
        OmacData_Rad <- matrix(NA, 2160, 1600)
        OmacDataSpec_SLC_Rad <- matrix(NA, 2160, 1600)
        OmacDataSpec_Rad <- matrix(NA, 2160, 1600)
        
        start_time <- Sys.time()
        
        for (i in 1:dim(OmacData_SLC_DNMS)[1]) { # Spectral
          for (j in 1:dim(OmacData_SLC_DNMS)[2]) { # Spatial
            
            # Define the x variable as one pixel for each DN/ms dataset
            tempOD_SLC_DNMS <- OmacData_SLC_DNMS[i,j]
            tempOD_DNMS <- OmacData_DNMS[i,j]
            tempSpec_SLC_DNMS <- OmacDataSpec_SLC_DNMS[i,j]
            tempSpec_DNMS <- OmacDataSpec_DNMS[i,j]
            
            # Make vector coefficients for that pixel
            temp_coef <- c(coef2incpt[i,j], coef2slp1[i,j], coef2slp2[i,j], coef2slp3[i,j], 0, 0) 
            
            # Make design matrix for each dataset
            design_OD_SLC_DNMS <- dummy_matrix(x.values = tempOD_SLC_DNMS, x_names = "x", psi.est = TRUE, 
                                        nameU = nameU1, nameV = nameV1, diffSlope = c(coef2slp2[i,j], coef2slp3[i,j]), 
                                        est.psi = c(coef2break1[i,j], coef2break2[i,j]))
            design_OD_DNMS <- dummy_matrix(x.values = tempOD_DNMS, x_names = "x", psi.est = TRUE, 
                                        nameU = nameU1, nameV = nameV1, diffSlope = c(coef2slp2[i,j], coef2slp3[i,j]), 
                                        est.psi = c(coef2break1[i,j], coef2break2[i,j]))
            design_Spec_SLC_DNMS <- dummy_matrix(x.values = tempSpec_SLC_DNMS, x_names = "x", psi.est = TRUE, 
                                        nameU = nameU1, nameV = nameV1, diffSlope = c(coef2slp2[i,j], coef2slp3[i,j]), 
                                        est.psi = c(coef2break1[i,j], coef2break2[i,j]))
            design_Spec_DNMS <- dummy_matrix(x.values = tempSpec_DNMS, x_names = "x", psi.est = TRUE, 
                                        nameU = nameU1, nameV = nameV1, diffSlope = c(coef2slp2[i,j], coef2slp3[i,j]), 
                                        est.psi = c(coef2break1[i,j], coef2break2[i,j]))

            
            # Estimate the radiance for one pixel
            temp_rad_OD_SLC_DNMS <- temp_coef %*% t(design_OD_SLC_DNMS)
            temp_rad_OD_DNMS <- temp_coef %*% t(design_OD_DNMS)
            temp_rad_Spec_SLC_DNMS <- temp_coef %*% t(design_Spec_SLC_DNMS)
            temp_rad_Spec_DNMS <- temp_coef %*% t(design_Spec_DNMS)
            
            # Compile radiance into a matrix
            OmacData_SLC_Rad[i,j] <- temp_rad_OD_SLC_DNMS
            OmacData_Rad[i,j] <- temp_rad_OD_DNMS
            OmacDataSpec_SLC_Rad[i,j] <- temp_rad_Spec_SLC_DNMS
            OmacDataSpec_Rad[i,j] <- temp_rad_Spec_DNMS
            
          }
        }
        
        end_time <- Sys.time()
        
        end_time - start_time
        
        print(paste("there are ", sum(is.na(OmacData_SLC_Rad)), "NA values")) # If there are NA values, check radiance coefficient inputs
        
        
    # then, convert from mW/cm2/sr/um to W/m2/sr/nm (W/1000mW * 10000cm2/m2 * 1um/1000nm = 0.01 so divide by 100)
        OmacData_SLC_Rad_W <- OmacData_SLC_Rad/100
        OmacData_Rad_W <- OmacData_Rad/100
        
        OmacDataSpec_SLC_Rad_W <- OmacDataSpec_SLC_Rad/100
        OmacDataSpec_Rad_W <- OmacDataSpec_Rad/100

        
    # Average across spatial pixels for spectralon (solar) and vegetation
    # (Note, spatial pixels 1023:1027 were same pixels as used in the subset for the SVD-based retrieval, processed in Subset_Radiance_for_Frankenberg.m)
        SolarData_subset <- OmacData_Rad_W[,1023:1027] # spatial pixel 1025 is in middle of most reflective spectralon level
        SolarData <- apply(SolarData_subset, 1, mean)
        SolarData_SLC_subset <- OmacData_SLC_Rad_W[,1023:1027]
        SolarData_SLC <- apply(SolarData_SLC_subset, 1, mean)
        
    # Just use pixel ## for vegetation spectrum
        VegData <- OmacData_Rad_W[,625] # spatial pixel 625 is in middle of plant
        VegData_SLC <- OmacData_SLC_Rad_W[,625]
    
          
          
          
 
    
##### PLOTS ###############################


## choose x-values
    xVals <- Wavelengths[,2]
    
    
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
    
    tiff(width=7, height=5, file= paste(path_figs, "Figure_Example_Laser_line_corr.tiff", sep = ""), units="in", res=600)
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
    
    
    

### Make plot showing stray light magnitude (for diffuser and vegetation)
    
    ## Define uncorrected color
    uncorr_col <- "#fc9272" # "grey70" #"#4292c6"
    
    # Define text for use in graph
    radiance_text = expression('Radiance (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')
    
    cols_native_SL <- c("SL in solar"="black","SL in vegetation"="grey55","SIF"="#EE604A")
    
    Fig_original_SL <- ggplot() +
      geom_line(data = data.frame(x = xVals, y = (SolarData - SolarData_SLC)), aes(x=x, y=y, colour = "SL in solar"), alpha = 0.85, size = 0.6) +
      geom_line(data = data.frame(x = xVals, y = (VegData - VegData_SLC)), aes(x=x, y=y, colour = "SL in vegetation"), alpha = 0.85, size = 0.6) +
      theme_classic() +
      scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
      xlab("Wavelength (nm)") +
      ylab(radiance_text) +
      scale_colour_manual(values=cols_native_SL, name = "") +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.position="top", 
            legend.direction="horizontal") + 
      geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)
    
    ggsave(filename = paste(path_figs, "Fig_original_SL.pdf", sep = ""), 
           plot = Fig_original_SL, width = 8, height = 4, units = "in", dpi = 300)
    
    

#### EXPORTS #################################

write.ENVI(OmacData_Rad_W, file = "./Outputs/OmacData_Rad_W")
write.ENVI(OmacData_SLC_Rad_W, file = "./Outputs/OmacData_SLC_Rad_W")

