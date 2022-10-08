# Loren Albert
# Feb 2021
# This figure is for the stray light simulations paper (Fig. 1), and is based on a figure I had made in make_SLC_plots.R
# This figure shows the magnitude of spectral stray light on vegetation and reference spectra versus the magnitude of
# SIF.  It uses the SIF levels exported from make_stray_light_simulations_v2.R.  There is a similar figure (but with
# vegetation spectra simulated instead of measured) created by make_stray_light_simulations_v2.R
#
# Copyright (C) 2022  Loren Albert
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
# and associated documentation files (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

### Housekeeping --------------------------

## Remove workspace objects if necessary
# rm(list=ls())

### Load packages
library(ggplot2)
library(caTools)

## define path to figure outputs
path_figs <- "/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/Figures"

## define path for spectra input
path_spectra <- "/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/Spectra_output/"

## define path for wavelengths file
path_data_wave <- "/Volumes/Research/IBES_KellnerLab/Shared/NIST_Jan_2018_organized/"

# Import SIF spectra (output from make_stray_light_simulations_v2)
SIFl1 <- read.csv(file = paste(path_spectra, "SIF_SIFl1.csv", sep = ""))
SIFl2 <- read.csv(file = paste(path_spectra, "SIF_SIFl2.csv", sep = ""))
SIFl3 <- read.csv(file = paste(path_spectra, "SIF_SIFl3.csv", sep = ""))

# Import stray light spectra (output from make_SLC_plots.R)
OmacData_Rad_W <- read.ENVI(file = paste(path_spectra, "OmacData_Rad_W", sep = ""))
OmacData_SLC_Rad_W <- read.ENVI(file = paste(path_spectra, "OmacData_SLC_Rad_W", sep = ""))

# Import wavelengths values for each spectral pixel
Wavelengths <- read.csv(paste(path_data_wave,"predicted_spectral_pixels_2018.csv",sep=""))



#### Processing ------------------------------

# Drop extra dimension in SIF spectra
SIFl1 <- SIFl1[,1]
SIFl2 <- SIFl2[,1]
SIFl3 <- SIFl3[,1]

# Average across spatial pixels for spectralon (solar) and vegetation
# (Note, spatial pixels 1023:1027 were same pixels as used in the subset for the SVD-based retrieval, processed in Subset_Radiance_for_Frankenberg.m)
SolarData_subset <- OmacData_Rad_W[,1023:1027] # spatial pixel 1025 is in middle of most reflective spectralon level
SolarData <- apply(SolarData_subset, 1, mean)
SolarData_SLC_subset <- OmacData_SLC_Rad_W[,1023:1027]
SolarData_SLC <- apply(SolarData_SLC_subset, 1, mean)

# Just use pixel ## for vegetation spectrum
VegData <- OmacData_Rad_W[,625] # spatial pixel 625 is in middle of plant
VegData_SLC <- OmacData_SLC_Rad_W[,625]




### Plot -----------------------------------

### Make plot showing SIF magnitude relative to stray light magnitude (for diffuser and vegetation)

### Set working directory
setwd(path_figs)

## choose x-values
xVals <- Wavelengths[,2]

## Define uncorrected color
uncorr_col <- "#fc9272" # "grey70" #"#4292c6"

# Define text for use in graph
radiance_text = expression('Radiance (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')')

tiff(width=6, height=10, file="Figure_Solar and veg spectra_v3.tiff", units = "in", res=600)  
par(mfrow=c(4,1), mar=c(4,5,1,1), oma=c(1,2,0,0), bty = 'n')

## Solar spectrum
plot(x=xVals, y=SolarData,type='l',
     xlab =NA,
     ylab =radiance_text,
     ylim =c(0,0.45),
     xlim = c(669, 781),
     xaxs = "i",
     xaxt = 'n',
     cex.lab = 1.25,
     # main = "Solar spectrum",
     col = uncorr_col)
Axis(side=1, labels=FALSE)
lines(x=xVals, y=SolarData_SLC, col="black")
legend(x=720,y=0.25,
       c("Uncorrected","Corrected"),
       col=c(uncorr_col,"black"),
       lty=1,
       bty='n')
# text("A", x=670, y= 0.43)

## Vegetation spectrum
plot(x=xVals, y=VegData, type='l',
     xlab=NA,
     ylab=radiance_text,
     ylim=c(0,0.16),
     xlim = c(669, 781),
     xaxs="i",
     xaxt='n',
     cex.lab = 1.25,
     # main = "Vegetation spectrum",
     col=uncorr_col)
Axis(side=1, labels=FALSE)
lines(x=(xVals), y=VegData_SLC, col="black")
# text("B",x=670,y=0.14)

# SIF correction effect in radiance units
SIF_color <- "#cb181d"
veg_color <- "#1a9641"
solar_color <- "#5e4fa2"
plot(x=xVals, y=(VegData - VegData_SLC), type="l",
     col = veg_color,
     ylim = c(0, max(SIFl3)),
     xlab = "Wavelength (nm)",
     ylab = radiance_text,
     xlim = c(669, 781),
     cex.lab = 1.25,
     xaxs ="i")
Axis(side=1, labels=FALSE)
# text("C", x = 670, y = max(SIFl3) - 0.1e-3)
legend(x = 670, y = max(SIFl3) - 0.05e-3, c(as.expression(bquote("SIF"[757]~"="~ .(round(SIFl3[1772], digits = 4)))),
                                          as.expression(bquote("SIF"[757]~"="~.(round(SIFl2[1772], digits = 4)))),
                                          as.expression(bquote("SIF"[757]~"="~.(round(SIFl1[1772], digits = 4)))),
                                          as.expression(bquote("Stray light in veg. spectrum")),
                                          as.expression(bquote("Stray light in solar spectrum"))),
       col=c(adjustcolor(col=SIF_color,alpha.f = 1),
             adjustcolor(col=SIF_color,alpha.f = 0.7),
             adjustcolor(col=SIF_color,alpha.f = 0.4),
             veg_color,
             solar_color),
       lty=1,bty='n')
lines(x = xVals, y = SIFl1, col=adjustcolor(col=SIF_color, alpha.f = 0.4))
lines(x = xVals, y = SIFl2, col=adjustcolor(col=SIF_color, alpha.f = 0.7))
lines(x = xVals, y = SIFl3, col=adjustcolor(col=SIF_color, alpha.f = 1))
lines(x = xVals, y = (SolarData - SolarData_SLC), col = solar_color)

# mtext("Wavelength (nm)", side=1, outer=T)
# mtext(radiance_text, side=2, outer=T, at = .65)


dev.off()

### Make plot just showing stray light magnitude (for diffuser and vegetation)

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

ggsave(filename = paste("Fig_original_SL.pdf", sep = ""), 
       plot = Fig_original_SL, width = 8, height = 4, units = "in", dpi = 300)
