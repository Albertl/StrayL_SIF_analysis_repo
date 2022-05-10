# Functions used in spectral stray light analysis (taken from The_CFL002_Fcns.R)

# Contents:
# SLC_fun_colwise: Function for stray light correction (columns of vectors of spectral pixels)
# dummy_matrix: Function for making matrix of dummy variables from segmented.lm outputs (with two break points)
# function for reading in ascii files that have metadata at the bottom
# simulate_SL_fun_colwise: function to simulate stray light for multiple columns given the 'A' and 'Y_IB' matrices from Zong et al (2006)
# ez.read: function for reading in ascii files that have metadata at the bottom
# threeFLD: Function for 3FLD SIF retrieval
# sFLD: Function for original (standard) FLD retrieval for SIF


#### Define function for spectral stray light correction ------------------------------------------------------------------------------------------------
# Version that will loop through columns
SLC_fun_colwise <-function(data_mat, C_mat){
  YIB <- matrix(nrow=dim(data_mat)[1],ncol=dim(data_mat)[2],data=NA)
  for (i in 1:dim(data_mat)[2]){
    Y_meas <- data_mat[,i]
    YIB[,i] <-C_mat %*% Y_meas
  }
  return(YIB)
}


#### Define function for making matrix of dummy variables from segmented.lm outputs -----------------------------------------------------------------------------------------
# This is based on code from predict.segmented())
# Note that this works for segmented lm with TWO breakpoints and would need to be changed/adapted for different number of breakpoints
# See the script "Test_segmented_lm_prediction.R" for toy examples with this function

dummy_matrix <- function(x.values, x_names, psi.est = TRUE, nameU, nameV, diffSlope, est.psi) {
  # This function creates a model matrix with dummy variables for a segmented lm with two breakpoints.
  # Inputs:
  # x.values: the x values of the segmented lm
  # x_names: the name of the column of x values
  # psi.est: this is legacy from the predict.segmented function, leave it set to 'TRUE'
  # obj: the segmented lm object
  # nameU: names (class character) of 3rd and 4th coef, which are "U1.x" "U2.x" for lm with two breaks. Example: names(c(obj$coef[3], obj$coef[4]))
  # nameV: names (class character) of 5th and 6th coef, which are "psi1.x" "psi2.x" for lm with two breaks. Example: names(c(obj$coef[5], obj$coef[6]))
  # diffSlope: the coefficients (class numeric) with the slope differences; called U1.x and U2.x for lm with two breaks. Example: c(o$coef[3], o$coef[4])
  # est.psi: the estimated break points (class numeric); these are the estimated breakpoints from segmented.lm. Example: c(obj$psi[1,2], obj$psi[2,2])
  #
  n <- length(x.values)
  k <- length(est.psi)
  PSI <- matrix(rep(est.psi, rep(n, k)), ncol = k)
  newZ <- matrix(x.values, nrow = n, ncol = k, byrow = FALSE)
  dummy1 <- pmax(newZ - PSI, 0)
  if (psi.est) {
    V <- ifelse(newZ > PSI, -1, 0)
    dummy2 <- if (k == 1) 
      V * diffSlope
    else V %*% diag(diffSlope)
    newd <- cbind(x.values, dummy1, dummy2)
    colnames(newd) <- c(x_names, nameU, nameV)
  } else {
    newd <- cbind(x.values, dummy1)
    colnames(newd) <- c(x_names, nameU)
  }
  # if (!x_names %in% names(coef(obj.seg))) 
  #   newd <- newd[, -1, drop = FALSE]
  
  col1 <- matrix(1, nrow = dim(newd)[1])
  newd2 <- cbind(col1, newd) # Now test is the same as model.matrix(o)
  return(newd2)
}


### Define function for reading in ascii files that have metadata at the bottom -----------------------------------------------------------------------------
# https://stackoverflow.com/questions/39110755/skip-specific-rows-using-read-csv-in-r

#' read csv table, wrapper of \code{\link{read.csv}}
#' @description read csv table, wrapper of \code{\link{read.csv}}
#' @param tolower whether to convert all column names to lower case
#' @param skip.rows rows to skip (1 based) before read in, eg 1:3
#' @return returns a data frame
#' @export
ez.read = function(file, ..., skip.rows=NULL, tolower=FALSE){
  if (!is.null(skip.rows)) {
    tmp = readLines(file)
    tmp = tmp[-(skip.rows)]
    tmpFile = tempfile()
    on.exit(unlink(tmpFile))
    writeLines(tmp,tmpFile)
    file = tmpFile
  }
  result = read.csv(file, ...)
  if (tolower) names(result) = tolower(names(result))
  return(result)
}


### Define function to simulate stray light for multiple columns given the 'A' and 'Y_IB' matrices from Zong et al (2006) ------------------------------------
# Want simulated 'D' matrix in Zong et al 2006, and an identity matrix 'I'
# Then simulate measurements in the presence of stray light as Y_meas = (I+D)*Y_IB
# Which is equivalent to Y_meas = A * Y_IB

simulate_SL_fun_colwise <- function(data_mat, A_mat){
  Y_sim <- matrix(nrow=dim(data_mat)[1],ncol=dim(data_mat)[2],data=NA)
  for (i in 1:dim(data_mat)[2]){
    Y_IB <- data_mat[,i]
    Y_sim[ ,i] <-A_mat %*% Y_IB
  }
  return(Y_sim)
}


### Define function for 3FLD-----------------------------------------------------------------------------------------------------------------


threeFLD <- function(lambdaL, lambdaR, lambdaIN, LoutL, LoutR, Lin, EoutL, EoutR, Ein){
  # Function inputs are:
  # lambdaL: wavelength position outside (left) of the band
  # LambdaR: wavelength position outside (right) of the band
  # LambdaIN: wavelength position inside the band
  # LoutL: canopy radiance outside (left) of the band
  # LoutR: canopy radiance outside (right) of the band
  # Lin: canopy radiance inside the band
  # EoutL: ground radiance (or irradiance, but then would need to add conversion by pi to this function) (left) of the band
  # EoutR: ground radiance (or irradiance, but then would need to add conversion by pi to this function) (right) of the band
  # Ein: ground radiance (or irradiance, but then would need to add conversion by pi to this function) inside the band
  # Refs: Damm et al (2011) RSE, Meroni et al (2009) RSE
  
  w_21 <- (lambdaR - lambdaIN)/(lambdaR-lambdaL)
  w_22 <- (lambdaIN - lambdaL)/(lambdaR - lambdaL)
  denom <- 1 - (Ein/(w_21*EoutL+w_22*EoutR))
  numer <- Lin - (Ein/(w_21*EoutL+w_22*EoutR)) * (w_21*LoutL + w_22*LoutR)
  SIF <- numer/denom
  if(Ein > EoutL){warning(paste('check inputs to threeFLD(): Ein as ', Ein, 'and EoutL as', EoutL, 'because Ein > EoutL'))}
  if(Ein > EoutR){warning(paste('check inputs to threeFLD(): Ein as ', Ein, 'and EoutR as', EoutR, 'because Ein > EoutR'))}
  if(Lin > LoutL){warning(paste('check inputs to threeFLD(): Lin as ', Lin, 'and LoutL as', LoutL, 'because Lin > LoutL'))}
  if(Lin > LoutR){warning(paste('check inputs to threeFLD(): Lin as ', Lin, 'and LoutR as', LoutR, 'because Lin > LoutR'))}
  return(SIF)
}


### Define function for original (standard) FLD --------------------------------------------------------------------------------------------------------


sFLD <- function(Eout, Ein, Lout, Lin){
  # Function inputs are:
  # Eout: ground radiance outside the band (assumes E, typically solar irradiance, is now solar radiance because we use the reflected flux from spectralon to represent E)
  # Ein: ground radiance inside the band (assumes E, typically solar irradiance, is now solar radiance because we use the reflected flux from spectralon to represent E)
  # Lout: target radiance outside the band
  # Lin: target radiance inside the band
  # Refs: Damm et al (2011) RSE, Meroni et al (2009) RSE
  
  numer <- Lin - (Ein/Eout * Lout)
  denom <- 1 - (Ein/Eout)
  SIF <- numer/denom
  return(SIF)
  
}


