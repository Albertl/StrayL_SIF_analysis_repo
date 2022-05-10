# This figure is for the stray light simulations paper, and shows the error attributable
# to stray light in Radiance units for all retrievals examined. The retrievals are
# exported from the following scripts:
# SFMresults_summary.R for "SFM_stray_no_stray_differences.csv"
# SVDresults_summary.R for "SVD_stray_no_stray_differences.csv"
# make_stray_light_simulations_v2.R for threeFLD_SIFl1_diffs, threeFLD_SIFl2_diffs and threeFLD_SIFl3_diffs
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


### Housekeeping --------------------------

## Remove workspace objects if necessary
# rm(list=ls())

### Load packages
library(ggplot2)
library(viridis)
library(caTools)

## Change working directory if necessary
# setwd('/Users/lalbert/Dropbox (Brown)/CFL_002_calib paper/Stray light simulations paper/SL_SIF_analysis_repo/')

## define paths
path_figs <- "./Figures/"
path_spectra <- "./Outputs/"


# Import
SFM_diffs <- read.csv(file = paste(path_spectra, "SFM_stray_no_stray_differences.csv", sep = ""))
SVD_diffs <- read.csv(file = paste(path_spectra, "SVD_stray_no_stray_differences.csv", sep = ""))
threeFLD_SIFl1_diffs <- read.csv(file = paste(path_spectra, "3FLD_SIFl1_stray_no_stray_differences.csv", sep = ""))
threeFLD_SIFl2_diffs <- read.csv(file = paste(path_spectra, "3FLD_SIFl2_stray_no_stray_differences.csv", sep = ""))
threeFLD_SIFl3_diffs <- read.csv(file = paste(path_spectra, "3FLD_SIFl3_stray_no_stray_differences.csv", sep = ""))


## Define indices for 687.1nm and 760.5nm bands within full 2160 bands
Wavelengths <- SFM_diffs$V1
O2B_band_ind_full <- which(round(Wavelengths, digits = 2) == 687.1)
O2A_band_ind_full <- which(abs(Wavelengths - 760.5) == min(abs(Wavelengths - 760.5)))
SFM_O2_inds <- c(O2B_band_ind_full, O2A_band_ind_full)

### Plots of stray light error for sFLD and 3FLD ----------------------------------------

#colsRet <- c("sFLD"="#3C0E66", "3FLD"="#058279")
colsRet <- c("3FLD" = plasma(4)[1], "SFM" = plasma(4)[2], "SVD" = plasma(4)[3])

#sif_error_ylim <- c(-0.001, 0.0025)
#sif_error_ylim <- c(NA, NA)
#sif_error_ylim <- c(-0.0001, 1.105e-03)
#sif_error_ylim <- c(-0.0002, 0.00115)
sif_error_ylim <- c(-0.00005, 0.00115)

# Define text for use in plot
radiance_text_error = expression('Stray light error (W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')       ') # spaces avoid y label cut-off
radiance_text_error = expression(atop('Stray light error', paste('(W m'^"-2"*' sr'^"-1"*' nm'^"-1"*')'))) # Units on line below label

make_SL_diffs_fig <- function(threeFLD_wl_vec, threeFLD_diff_vec, SFM_wl_vec, SFM_diff_vec, 
                              SVD_wl_vec, SVD_diff_vec, sif_ylim, ylab_text, colVals){
  ggplot() +
    geom_line(data = data.frame(x = threeFLD_wl_vec, y = threeFLD_diff_vec), aes(x=x, y=y, colour = "3FLD"), 
              alpha = 0.70, size = 1) +
    geom_point(data = data.frame(x = threeFLD_wl_vec, y = threeFLD_diff_vec), aes(x=x, y=y, colour = "3FLD"), 
               alpha = 0.70, size = 4.5, shape = 18) +
    geom_point(data = data.frame(x = SFM_wl_vec[SFM_O2_inds], y = SFM_diff_vec[SFM_O2_inds]), aes(x=x, y=y, colour = "SFM"), 
               alpha = 0.70, size = 3.3, shape = 15) +
    geom_line(data = data.frame(x = SVD_wl_vec, y = SVD_diff_vec), aes(x=x, y=y, colour = "SVD"),
              alpha = 0.70, size = 1.5) +
    # geom_point(data = data.frame(x = SVD_wl_vec, y = SVD_diff_vec), aes(x=x, y=y, colour = "SVD"), 
    #            alpha = 0.70, size = 4.5, shape = 18) +
    theme_classic() +
    scale_x_continuous(limits = c(669, 781), breaks=seq(670,780,20), expand = c(0, 0)) +
    scale_y_continuous(limits = sif_ylim) +
    xlab("Wavelength (nm)") +
    ylab(ylab_text) +
    scale_colour_manual(values = colVals, name = "",
                        guide = guide_legend(override.aes = list(linetype = c("blank", "blank", "dashed"), 
                                                                 shape = c(18, 15, NA)))) +
    theme(legend.position = c(0.79, 0.9), 
          legend.direction="vertical") + 
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)} 

# Very high stray light scenario & SIFl1: differences betw stray and no stray for Dmat_1_scenario with SIFl1
Plot_diffs_Dmat_1_SIFl1 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl1_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl1_diffs$three_FLD_SL_Dmat_1_SIFl1_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_1_SIFl1_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_1_SIFl1_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# Very high stray light scenario & SIFl2: differences betw stray and no stray for Dmat_1_scenario with SIFl2
Plot_diffs_Dmat_1_SIFl2 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl2_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl2_diffs$three_FLD_SL_Dmat_1_SIFl2_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_1_SIFl2_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_1_SIFl2_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# Very high stray light scenario & SIFl3: differences betw stray and no stray for Dmat_1_scenario with SIFl3
Plot_diffs_Dmat_1_SIFl3 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl3_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl3_diffs$three_FLD_SL_Dmat_1_SIFl3_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_1_SIFl3_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_1_SIFl3_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# High stray light scenario & SIFl1: differences betw stray and no stray for Dmat_2_scenario with SIFl1
Plot_diffs_Dmat_2_SIFl1 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl1_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl1_diffs$three_FLD_SL_Dmat_2_SIFl1_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_2_SIFl1_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_2_SIFl1_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# High stray light scenario & SIFl2: differences betw stray and no stray for Dmat_2_scenario with SIFl2
Plot_diffs_Dmat_2_SIFl2 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl2_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl2_diffs$three_FLD_SL_Dmat_2_SIFl2_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_2_SIFl2_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_2_SIFl2_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# High stray light scenario & SIFl3: differences betw stray and no stray for Dmat_2_scenario with SIFl3
Plot_diffs_Dmat_2_SIFl3 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl3_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl3_diffs$three_FLD_SL_Dmat_2_SIFl3_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_2_SIFl3_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_2_SIFl3_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# Medium stray light scenario & SIFl1: differences betw stray and no stray for Dmat_3_scenario with SIFl1
Plot_diffs_Dmat_3_SIFl1 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl1_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl1_diffs$three_FLD_SL_Dmat_3_SIFl1_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_3_SIFl1_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_3_SIFl1_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# Medium stray light scenario & SIFl2: differences betw stray and no stray for Dmat_3_scenario with SIFl2
Plot_diffs_Dmat_3_SIFl2 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl2_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl2_diffs$three_FLD_SL_Dmat_3_SIFl2_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_3_SIFl2_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_3_SIFl2_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# Medium stray light scenario & SIFl3: differences betw stray and no stray for Dmat_3_scenario with SIFl3
Plot_diffs_Dmat_3_SIFl3 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl3_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl3_diffs$three_FLD_SL_Dmat_3_SIFl3_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_3_SIFl3_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_3_SIFl3_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# Low stray light scenario & SIFl1: differences betw stray and no stray for Dmat_4_scenario with SIFl1
Plot_diffs_Dmat_4_SIFl1 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl1_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl1_diffs$three_FLD_SL_Dmat_4_SIFl1_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_4_SIFl1_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_4_SIFl1_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# Low stray light scenario & SIFl2: differences betw stray and no stray for Dmat_4_scenario with SIFl2
Plot_diffs_Dmat_4_SIFl2 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl2_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl2_diffs$three_FLD_SL_Dmat_4_SIFl2_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_4_SIFl2_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_4_SIFl2_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

# Low stray light scenario & SIFl3: differences betw stray and no stray for Dmat_4_scenario with SIFl3
Plot_diffs_Dmat_4_SIFl3 <- make_SL_diffs_fig(threeFLD_wl_vec = threeFLD_SIFl3_diffs$wavelength, 
                                             threeFLD_diff_vec = threeFLD_SIFl3_diffs$three_FLD_SL_Dmat_4_SIFl3_diff,
                                             SFM_wl_vec = SFM_diffs$V1,
                                             SFM_diff_vec = SFM_diffs$SFM_SL_Dmat_4_SIFl3_diff,
                                             SVD_wl_vec = SVD_diffs$V1,
                                             SVD_diff_vec = SVD_diffs$SVD_SL_Dmat_4_SIFl3_diff,
                                             sif_ylim = sif_error_ylim,
                                             ylab_text = radiance_text_error,
                                             colVals = colsRet)

### Save plots ---------------------------------------------------------------------------

setwd(path_figs)

ggsave(filename = "Plot_diffs_Dmat_1_SIFl1.pdf",
       plot = Plot_diffs_Dmat_1_SIFl1, width = 10, height = 6, units = "cm", dpi = 400)
# If the lower y axis limit is < -0.00015 or so then this plot won't show some values (hence the warning for removed rows)
# For the main text multi-panel showing only SIF level 2 panels, this doesn't matter.

ggsave(filename = "Plot_diffs_Dmat_1_SIFl2.pdf",
       plot = Plot_diffs_Dmat_1_SIFl2, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_1_SIFl3.pdf",
       plot = Plot_diffs_Dmat_1_SIFl3, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_2_SIFl1.pdf",
       plot = Plot_diffs_Dmat_2_SIFl1, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_2_SIFl2.pdf",
       plot = Plot_diffs_Dmat_2_SIFl2, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_2_SIFl3.pdf",
       plot = Plot_diffs_Dmat_2_SIFl3, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_3_SIFl1.pdf",
       plot = Plot_diffs_Dmat_3_SIFl1, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_3_SIFl2.pdf",
       plot = Plot_diffs_Dmat_3_SIFl2, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_3_SIFl3.pdf",
       plot = Plot_diffs_Dmat_3_SIFl3, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_4_SIFl1.pdf",
       plot = Plot_diffs_Dmat_4_SIFl1, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_4_SIFl2.pdf",
       plot = Plot_diffs_Dmat_4_SIFl2, width = 10, height = 6, units = "cm", dpi = 400)

ggsave(filename = "Plot_diffs_Dmat_4_SIFl3.pdf",
       plot = Plot_diffs_Dmat_4_SIFl3, width = 10, height = 6, units = "cm", dpi = 400)
