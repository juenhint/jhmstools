
#'Get correlation matrix between selected features
#'
#' @param exprs a numeric data.frame of matrix. A dataframe of peak areas features\*samples (rows\*columns)
#' @param names vector of character strings. The names of the features (feature IDs)
#' @param ... all other params are passed to `cor()`
#'
#' @return a square matrix of correlations
#'
#' @export
#' @importFrom stats cor
get_cors <- function(exprs, names, use="complete.obs", ...) {
  return(cor(t(exprs[names,]), use=use, ...))
}

#' Plot histogram of peak areas
#'
#' @param exprs a numeric data.frame of matrix. A dataframe of peak areas features\*samples (rows\*columns)
#' @param var character string. The name of the feature (feature ID) to plot
#' @param log logical. Whether to log-transform
#' @param ... all other params are passed to `hist()`
#'
#' @export
#' @importFrom graphics hist
histo <- function(exprs, var, log=F, ...) {
  if (!log) hist(t(exprs[var,]), ...)
  else hist(log(t(exprs[var,])), ...)
}

#' Find featured by compound names
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param name character string. The compound name to search
#' @param column_to_search character string. The column containing annotations
#' @param ignore.case logical. Whether to ignore case in search
#'
#' @return a vector of strings. The feature IDs of the hits
#' @export
find_by_compound <- function(data, name, column_to_search="Metabolite.name", ignore.case=T) {
  inds <- grep(name, data[,column_to_search], ignore.case = ignore.case)
  return(rownames(data)[inds])
}

#' Get MS2 spectrum from feature data
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param var character string. The rowname of the feature (feature_ID) to extract
#' @param lim numeric. Relative abundance cutoff of MS2 peaks
#' @param itlim numeric. Relative abundance cutoff of isotope peaks
#' @param mz_col character string. Name of the m/z column
#' @param ms2_col character string. Name of the MS2 spectrm column
#' @param adduct_col character string. Name of the adduct column
#' @param isot_col character string. Name of the isotope spectrm column
#'
#' @return A list of 5 objects:
#'  precursor = m/z value of the precursor ion
#'  spectrum = dataframe of the MS2 spectru,
#'  adduct = type of adduct
#'  isotope_spectrum = dataframe of the isotope spectrum
#'  feature = the variable name
#' @export
#'
get_spectrum <- function(data, var, lim=0.01, itlim=0, mz_col="Average.Mz", ms2_col="MS.MS.spectrum", adduct_col="Adduct.type", isot_col="MS1.isotopic.spectrum") {
  precursor <- data[var,mz_col]
  spectrum <- data[var,ms2_col]
  adduct <- data[var,adduct_col]
  itspectrum <- get_isotope_spectrum(data, var, lim=itlim, mz_col=mz_col, adduct_col=adduct_col, isot_col=isot_col)[["spectrum"]]
  id <- var
  spectrum <- strsplit(spectrum,split = " ")[[1]]
  spectrum <- mapply(spectrum, FUN=strsplit, MoreArgs = list(split=":"))
  spectrum <- as.data.frame(apply(as.data.frame(spectrum), MARGIN=1, FUN=as.numeric))
  spectrum[[3]] <- spectrum[,2]/max(spectrum[,2])
  colnames(spectrum) <- c("mz","ab","rel")
  ret <- list(precursor=precursor, spectrum = spectrum[spectrum$rel>lim,], adduct=adduct,
              isotope_spectrum=itspectrum, feature=id)
  plot_spectrum(ret)
  return(ret)
}

#' Get isotope spectrum from feature data
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param var character string. The rowname of the feature (feature_ID) to extract
#' @param lim numeric. Relative abundance cutoff of isotope peaks
#' @param mz_col character string. Name of the m/z column
#' @param isot_col character string. Name of the isotope spectrm column
#' @param adduct_col character string. Name of the adduct column
#'
#' @return A list of 4 objects:
#'  precursor = m/z value of the precursor ion
#'  isotope_spectrum = dataframe of the isotope spectrum
#'  adduct = type of adduct
#'  feature = the variable name
#'
#' @export
#'
get_isotope_spectrum <- function(data, var, lim=0, mz_col="Average.Mz", isot_col="MS1.isotopic.spectrum", adduct_col="Adduct.type") {
  precursor <- data[var,mz_col]
  spectrum <- data[var,isot_col]
  adduct <- data[var,adduct_col]
  id <- var
  spectrum <- strsplit(spectrum,split = " ")[[1]]
  spectrum <- mapply(spectrum, FUN=strsplit, MoreArgs = list(split=":"))
  spectrum <- as.data.frame(apply(as.data.frame(spectrum), MARGIN=1, FUN=as.numeric))
  spectrum[[3]] <- spectrum[,2]/max(spectrum[,2])
  colnames(spectrum) <- c("mz","ab","rel")
  ret <- list(precursor=precursor, isotope_spectrum = spectrum[spectrum$rel>lim,], adduct=adduct,
              feature=id)
  plot_isotope(ret)
  return(ret)
}


#' Plot ms2 spectrum
#'
#' @param sp list. The spectrum returned by `get_spectrum`
#' @param alpha numeric. Transparency of bars
#' @param ... all other args are passed to `geom_col`
#'
#' @return ggplot object. The mass spec plot.
#' @export
#' @import ggplot2
plot_spectrum <- function(sp, alpha=0.6, fill="blue3", width = width, ...) {
  requireNamespace("ggplot2")
  p <- ggplot(sp[["spectrum"]], aes(x=mz, y=rel)) + geom_col(fill = fill, alpha=alpha, width = 1, ...) +
    geom_hline(yintercept = 0) + geom_text(aes(label=mz),  hjust = 0, size=4, alpha=alpha) +
    labs(caption=paste0("Precursor m/z ",sp[["precursor"]]))
  print(p)
  return(p)
}


#' Plot isotope spectrum
#'
#' @param sp list. The spectrum returned by `get_spectrum` or `get_isotope_spectrum`
#' @param alpha numeric. Transparency of bars
#' @param ... all other args are passed to `geom_col`
#'
#' @return ggplot object. The mass spec plot.
#' @export
#' @import ggplot2
plot_isotope <- function(sp, alpha=0.6, fill="blue3", width = width, ...) {
  requireNamespace("ggplot2")
  p <- ggplot(sp[["isotope_spectrum"]], aes(x=mz, y=rel)) + geom_col(fill=fill, alpha=alpha, width=1, ...) +
    geom_hline(yintercept = 0) + geom_text(aes(label=mz),  hjust = 0, size=4, alpha=alpha) +
    labs(caption=paste0("Precursor m/z ",sp[["precursor"]]))
  print(p)
  return(p)
}


#' Find adducts for molecular mass
#'
#' @param mass numeric. Molecular mass of the precursor
#' @param mode character string. Which polarization to use.
#'
#' @return data.frame. The m/z values calculated for different adducts
#' @export
#'
find_adducts <- function(mass, mode="pos") {
  #load("Data/adducts_pos.Rda")
  #load("Data/adducts_neg.Rda")
  addcts <- list(adducts.neg,adducts.pos)[[grepl("pos", mode, ignore.case = T)+1]]
  vals <- data.frame(mz=as.numeric(addcts$Mass + mass*addcts$Mult), row.names=addcts$Ion.name)
  return(vals)
}

#' Find molecular weights for adduct
#'
#' @param feature_data a data.frame object. The feature data of the MS peaks
#' @param var character string. Name of the feature to extract
#' @param mode character string. Optional. Which polarization to use. Extracts from feature name by default.
#' @param mz_col character string. Name of the m/z column
#'
#' @return data.frame. The molecular weights derived from different assumed adducts
#' @export
find_mws <- function(feature_data, var, mode=var, mz_col="Average.Mz") {
  mz = feature_data[var,mz_col]
  #load("R_files/adducts_pos.Rda")
  #load("R_files/adducts_neg.Rda")
  addcts <- list(adducts.neg,adducts.pos)[[grepl("pos", mode, ignore.case = T)+1]]
  vals <- data.frame(calculated_mass=as.numeric((mz-addcts$Mass)/addcts$Mult), row.names=addcts$Ion.name)
  return(vals)
}

#' Find molecular weights for adduct m/z
#'
#' @param mz numeric. m/z of the precursor
#' @param mode character string. Which polarization to use.
#'
#' @return data.frame. The molecular weights derived from different assumed adducts
#' @export
find_mws2 <- function(mz, mode) {
  #load("R_files/adducts_pos.Rda")
  #load("R_files/adducts_neg.Rda")
  addcts <- list(adducts.neg,adducts.pos)[[grepl("pos", mode, ignore.case = T)+1]]
  vals <- data.frame(calculated_mass=as.numeric((mz-addcts$Mass)/addcts$Mult), row.names=addcts$Ion.name)
  return(vals)
}


#' Find shared molecular weights for two adducts
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param var1 character string. Name of the 1st feature to extract
#' @param var2 character string. Name of the 2nd feature to extract
#' @param mode1 character string. Optional. Polarization of 1st adduct. Extracts from feature name by default.
#' @param mode2 character string. Optional. Polarization of 2nd adduct. Extracts from feature name by default.
#' @param tol Δ m/z tolerance for adducts
#' @param use.mass2adduct whether to also get adduct data from the package mass2adduct
#'
#' @export
#'
find_common_mw <- function(data, var1, var2, mode1=var1, mode2=var2, tol=0.005, use.mass2adduct=T) {
  mws1 <- find_mws(data, var1, mode1)
  mws2 <- find_mws(data, var2, mode2)
  name_m <- mapply(rownames(mws1), FUN=function(x) {paste0("var1:",x," var2:", rownames(mws2))})
  diffm <- mapply(mws1$calculated_mass, FUN=function(x) {abs(mws2$calculated_mass - x)})
  hits <- cbind(diffm[diffm<=tol],name_m[diffm<=tol])
  print( hits[order(hits[,1]),] )

  if (use.mass2adduct) {
    if (!requireNamespace("mass2adduct",quietly = T))
    {
      stop("Package \"mass2adduct\" is required for this function")
    }
    diffs <- mass2adduct::adducts
    diffs[["diff"]] <- abs(diffs$mass - abs(data[var1,"Average.Mz"] - data[var2,"Average.Mz"]))
    diffs <- subset(diffs, diff <= tol)
    print(diffs)
  }
}

#'Find features by presence of MS2 peaks
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param peaks numeric or vector of numeric values. The m/z values of the MS2 peaks to search
#' @param operator c(all, any) operator to use with multiple given peaks. Defaults to all.
#' @param tol numeric. m/z tolerance for peak matching
#' @param lim numeric. Relative abundance cutoff of peaks
#' @param ms2_col character string. Name of the MS2 spectrum column
#'
#' @return a vector of hits (feature IDs)
#'
#' @export
find_by_MS2peaks <- function(data, peaks, operator=all, tol=0.005, lim=0.03, ms2_col="MS.MS.spectrum") {
  hits <- c()
  for (r in rownames(data)) {
    spectrum <- data[r,ms2_col]
    if (is.na(spectrum)) next
    spectrum <- strsplit(spectrum,split = " ")[[1]]
    spectrum <- mapply(spectrum, FUN=strsplit, MoreArgs = list(split=":"))
    spectrum <- as.data.frame(apply(as.data.frame(spectrum), MARGIN=1, FUN=as.numeric))
    spectrum[[3]] <- spectrum[,2]/max(spectrum[,2])
    colnames(spectrum) <- c("mz","ab","rel")
    spectrum <- spectrum[spectrum$rel>lim,]
    peaks_hits <- mapply(peaks, FUN=function(p) any(abs(p - spectrum$mz) < tol))
    if (operator(peaks_hits)) hits <- c(hits, r)
  }
  return(hits)
}
