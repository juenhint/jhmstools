
#' Read data from an MS-Dial format excel file (feature data \* abundances \* sample data)
#'
#' @param file string. The name of the .xlsx-file to read.
#' @param feature_IDs_present logical. Whether the first column contains unique feature-IDs for the variables
#'
#' @return A list with 3 elements (data frames):
#'  Expressions or peak integrals
#'  Metadata for the features
#'  Metadata for the samples
#'
#' @export
#'
#' @import openxlsx
read_tables_from_MSDIALfile <- function(file, feature_IDs_present=TRUE) {
  test <- openxlsx::read.xlsx(file, colNames = FALSE, rowNames = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  range.ecols <- 1:(table(is.na(test[1,]))[["TRUE"]] + 1)
  range.dcols <- table(is.na(test[1,]))[["TRUE"]]+1:dim(test)[[2]]
  range.erows <- 1:(table(is.na(test[,1]))[["TRUE"]] + 1)
  range.drows <- table(is.na(test[,1]))[["TRUE"]]+1:dim(test)[[1]]
  feature_data <- openxlsx::read.xlsx(file, colNames = TRUE, rowNames = feature_IDs_present, rows = range.drows, cols=range.ecols)
  exprs <- openxlsx::read.xlsx(file, colNames = TRUE, rowNames = FALSE, rows = range.drows, cols=range.dcols[-1])
  rownames(exprs) <- rownames(feature_data)
  pheno_data <- openxlsx::read.xlsx(file, colNames = FALSE, rowNames = TRUE, rows = range.erows, cols=range.dcols)
  colnames(pheno_data) <- pheno_data[range.drows[1],]
  pheno_data <- as.data.frame(t(pheno_data))
  if (!all(colnames(exprs) == rownames(pheno_data))) print("SAMPLE NAMES DO NOT MATCH")
  if (!all(rownames(exprs) == rownames(feature_data))) print("FEATURE NAMES DO NOT MATCH")
  return(list(exprs=exprs, feature_data=feature_data, pheno_data=pheno_data))
}

#'Get correlation matrix between selected features
#'
#' @param exprs a numeric data.frame of matrix. A dataframe of peak areas features\*samples (rows\*columns)
#' @param names vector of character strings. The names of the features (feature IDs)
#' @param r_lim  numeric. The abs(R) threshold for displaying the value in the lower triangle
#' @param ... all other params are passed to `cor()`
#'
#' @return a square matrix of correlations, lower triangle only contains the values with abs(R) higher than `r_lim`
#'
#' @export
#' @importFrom stats cor
get_cors <- function(exprs, names, r_lim = 0.8, ...) {
  ARGS <- replace(alist(x = t(exprs[names,]), use="complete.obs"), names(list(...)), values = list(...))
  cors <- do.call(cor, ARGS)
  for (r in 1:dim(cors)[1]) {
    for (c in (0+r):dim(cors)[2]) {
      if (abs(cors[r,c]) < r_lim) cors[c,r] = NA
    }
  }
  return(cors)
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
histo <- function(exprs, var, log=FALSE, ...) {
  if (!log) hist(t(exprs[var,]), ...)
  else hist(log10(t(exprs[var,])), ...)
}

#' Find features by annotated names
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param name regular expression. The compound name to search
#' @param column_to_search regular expression. The column containing annotations
#' @param ignore.case logical. Whether to ignore case in search
#'
#' @return a vector of strings. The feature IDs of the hits
#' @export
find_by_compound <- function(data, name, column_to_search=".*Metabolite.*name", ignore.case=TRUE) {
  column_to_search<-grep(column_to_search, colnames(data))[1]
  inds <- grep(name, data[,column_to_search], ignore.case = ignore.case)
  return(rownames(data)[inds])
}

#' Get MS2 spectrum from feature data
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param var character string. The rowname of the feature (feature_ID) to extract
#' @param lim numeric. Relative abundance cutoff of MS2 peaks
#' @param itlim numeric. Relative abundance cutoff of isotope peaks
#' @param mz_col regular expression. Name of the m/z column
#' @param ms2_col regular expression. Name of the MS2 spectrm column
#' @param adduct_col regular expression. Name of the adduct column
#' @param isot_col regular expression. Name of the isotope spectrm column
#'
#' @return A list of 5 objects:
#'  precursor = m/z value of the precursor ion
#'  spectrum = dataframe of the MS2 spectru,
#'  adduct = type of adduct
#'  isotope_spectrum = dataframe of the isotope spectrum
#'  feature = the variable name
#' @export
#'
get_spectrum <- function(data, var, lim=0.01, itlim=0, mz_col=".*Mz.*", ms2_col=".*MS.*MS.*spectrum", adduct_col=".*Adduct.*", isot_col=".*isotopic.*spectrum") {
  mz_col <- grep(mz_col, colnames(data), ignore.case = T, value = T)[1]
  ms2_col <- grep(ms2_col, colnames(data), ignore.case = T, value = T)[1]
  adduct_col <- grep(adduct_col, colnames(data), ignore.case = T, value = T)[1]
  isot_col <- grep(isot_col, colnames(data), ignore.case = T, value = T)[1]
  precursor <- data[var,mz_col]
  spectrum <- data[var,ms2_col]
  adduct <- data[var,adduct_col]
  itspectrum <- get_isotope_spectrum(data, var, lim=itlim, mz_col=mz_col, adduct_col=adduct_col, isot_col=isot_col)[["isotope_spectrum"]]
  id <- var

  spectrum <- strsplit(spectrum,split = " ")[[1]]
  spectrum <- mapply(spectrum, FUN=strsplit, MoreArgs = list(split=":"))
  spectrum <- matrix(apply(as.data.frame(spectrum), MARGIN=1, FUN=as.numeric), ncol=2)
  spectrum <- as.data.frame(cbind(spectrum, spectrum[,2]/max(spectrum[,2])))

  colnames(spectrum) <- c("mz","ab","rel")
  ret <- list(precursor=precursor, spectrum = spectrum[spectrum$rel>lim,], adduct=adduct,
              isotope_spectrum=itspectrum, feature=id)
  try(plot_spectrum(ret))
  return(ret)
}

#' Get isotope spectrum from feature data
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param var character string. The rowname of the feature (feature_ID) to extract
#' @param lim numeric. Relative abundance cutoff of isotope peaks
#' @param mz_col regular expression. Name of the m/z column
#' @param isot_col regular expression. Name of the isotope spectrm column
#' @param adduct_col regular expression. Name of the adduct column
#'
#' @return A list of 4 objects:
#'  precursor = m/z value of the precursor ion
#'  isotope_spectrum = dataframe of the isotope spectrum
#'  adduct = type of adduct
#'  feature = the variable name
#'
#' @export
#'
get_isotope_spectrum <- function(data, var, lim=0, mz_col=".*Mz", isot_col=".*isotopic.*spectrum", adduct_col=".*Adduct.*") {
  mz_col <- grep(mz_col, colnames(data), ignore.case = T, value = T)[1]
  adduct_col <- grep(adduct_col, colnames(data), ignore.case = T, value = T)[1]
  isot_col <- grep(isot_col, colnames(data), ignore.case = T, value = T)[1]
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
  try(plot_isotope(ret))
  return(ret)
}


#' Plot ms2 spectrum
#'
#' @param sp list. The spectrum returned by `get_spectrum`
#' @param alpha numeric. Transparency of bars
#' @param fill string or RGB value. Color of bars
#' @param width numeric. Width of bars
#' @param ... all other args are passed to `geom_col`
#'
#' @return ggplot object. The mass spec plot.
#' @export
#' @import ggplot2
plot_spectrum <- function(sp, alpha=0.6, fill = "blue3", width = 1, ...) {
  requireNamespace("ggplot2")
  p <- ggplot(sp[["spectrum"]], aes(x=mz, y=rel)) + geom_col(fill=fill, alpha=alpha, width=width, ...) +
    geom_hline(yintercept = 0) + geom_text(aes(label=mz),  hjust = 0, size=4, alpha=alpha) +
    labs(caption=paste0("Precursor m/z ",sp[["precursor"]]))
  print(p)
  return(p)
}


#' Plot isotope spectrum
#'
#' @param sp list. The spectrum returned by `get_spectrum` or `get_isotope_spectrum`
#' @param alpha numeric. Transparency of bars
#' @param fill string or RGB value. Color of bars
#' @param width numeric. Width of bars
#' @param ... all other args are passed to `geom_col`
#'
#' @return ggplot object. The mass spec plot.
#' @export
#' @import ggplot2
plot_isotope <- function(sp, alpha=0.6, fill="blue3", width = 1, ...) {
  requireNamespace("ggplot2")
  p <- ggplot(sp[["isotope_spectrum"]], aes(x=mz, y=rel)) + geom_col(fill=fill, alpha=alpha, width=width, ...) +
    geom_hline(yintercept = 0) + geom_text(aes(label=mz),  hjust = 0, size=4, alpha=alpha) +
    labs(caption=paste0("Precursor m/z ",sp[["precursor"]]))
  print(p)
  return(p)
}


#' Find adducts for molecular mass
#'
#' @param mass numeric. Molecular mass of the precursor
#' @param mode regular expression. Which polarization to use.
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
#' @param data a data.frame object. The feature data of the MS peaks
#' @param var character string. Name of the feature to extract
#' @param mode regular expression. Optional. Which polarization to use. Extracts from feature name by default.
#' @param mz_col character string. Name of the m/z column
#'
#' @return data.frame. The molecular weights derived from different assumed adducts
#' @export
find_mws <- function(data, var, mode=var, mz_col="Average.*Mz") {
  mz_col <- grep(mz_col, colnames(data), ignore.case = T, value = T)[1]
  mz = data[var,mz_col]
  #load("R_files/adducts_pos.Rda")
  #load("R_files/adducts_neg.Rda")
  addcts <- list(adducts.neg,adducts.pos)[[grepl("pos", mode, ignore.case = T)+1]]
  vals <- data.frame(calculated_mass=as.numeric((mz-addcts$Mass)/addcts$Mult), row.names=addcts$Ion.name)
  return(vals)
}

#' Find molecular weights for adduct m/z
#'
#' @param mz numeric. m/z of the precursor
#' @param mode regular expression. Which polarization to use.
#'
#' @return data.frame. The molecular weights derived from different assumed adducts
#' @export
find_mws2 <- function(mz, mode="pos") {
  #load("R_files/adducts_pos.Rda")
  #load("R_files/adducts_neg.Rda")
  addcts <- list(adducts.neg,adducts.pos)[[grepl("pos", mode, ignore.case = T)+1]]
  vals <- data.frame(calculated_mass=as.numeric((mz-addcts$Mass)/addcts$Mult), row.names=addcts$Ion.name)
  return(vals)
}

#' Calculate accurate mass from a molecular formula
#'
#' @param mf character string. Molecular formula.
#' @param monoisotopic logical. Whether to return monoisotopic mass (default) or molar mass.
#'
#' @return numeric. The mass of the molecule. NA if any of the symbols are not recognized or contain non-stable elements.
#' @export
#'
find_mwfromformula <- function(mf, monoisotopic=TRUE) {
  #gsub("\\d", " ", mf)
  multipl <- as.numeric(regmatches(mf, regexpr("(\\d*)", mf)))
  multipl <- ifelse(is.na(multipl), 1, multipl)
  m <- gregexpr("([A-Z][a-z]*)(\\d*)", mf)
  elements <- regmatches(mf, m)[[1]]

  parsed_formula <- list()
  for (element in elements) {
    # Extract the element symbol and count
    symbol <- sub("\\d.*", "", element)
    count <- as.numeric(sub("\\D*", "", element))
    # If count is missing, 1
    if (is.na(count)) {
      count <- 1
    }
    parsed_formula[[symbol]] <- ifelse(is.null(parsed_formula[[symbol]]),yes = 0,no = parsed_formula[[symbol]]) + count
  }
  #multiply
  parsed_formula <- lapply(parsed_formula, FUN = function(x)x*multipl)
  masses <- c()
  #count the total mass, return na if unintelligible
  for (atom in names(parsed_formula)) {
    mass <- atoms[atom,"Mass"]*parsed_formula[[atom]]
    if (!monoisotopic) mass <- atoms[atom,"Atomic.weight"]*parsed_formula[[atom]]
    if (is.na(mass)) {
      print("Unparseable molecular formula")
      masses <- NA
      break
    }
    masses <- c(masses, mass)
  }
  return(sum(masses))
}

#' Find shared molecular masses for two molecular features
#'
#' @description
#' Checks whether the two molecular features have a shared molecular masses within given m/z tolerance.
#' Additionally, can use mass2adduct to find molecular fragments between the two features.
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param var1 character string. Name of the 1st feature to extract
#' @param var2 character string. Name of the 2nd feature to extract
#' @param mode1 regular expression. Optional. Polarization of 1st adduct. Extracts from feature name by default.
#' @param mode2 regular expression. Optional. Polarization of 2nd adduct. Extracts from feature name by default.
#' @param tol numeric. Δ m/z tolerance for adducts
#' @param use.mass2adduct boolean. Whether to also get fragment data from the package mass2adduct
#' @param mz_col regular expression. Name of the column containing m/z values
#'
#' @export
#' @return a data.frame. Suggested adducts for each feature, accuracy and calculated average molecular mass
#'
find_common_mw <- function(data, var1, var2, mode1=var1, mode2=var2, tol=0.005, use.mass2adduct=FALSE, mz_col= ".*mz") {
  mz_col <- grep(mz_col, colnames(data), ignore.case = T, value = T)[1]
  if (abs(data[var1,mz_col] - data[var2,mz_col]) < tol) return(paste0("Same m/z within tolerance of ",tol," Da"))
  mws1 <- find_mws(data, var1, mode1, mz_col = mz_col)
  mws2 <- find_mws(data, var2, mode2, mz_col = mz_col)

  diffm <- mapply(mws2$calculated_mass, FUN=function(x) {abs(mws1$calculated_mass - x)})
  rownames(diffm) <- rownames(mws1)
  colnames(diffm) <- rownames(diffm)
  avrg <- mapply(mws2$calculated_mass, FUN=function(x) {
    mapply(mws1$calculated_mass, FUN=function(y){mean(c(x,y))})
  })
  avrg <- as.data.frame(as.table(avrg))

  diffm <- as.data.frame(as.table(diffm))
  diffm[["A"]] <- avrg[["Freq"]]
  colnames(diffm) <- c(var1,var2,"Accuracy(Da)","Calculated_MW(Da)")

  #hits <- cbind(diffm[diffm<tol],name_m[diffm<tol],rownames(mws1)[diffm<tol], rownames(mws2)[diffm<tol])
  #print( hits[order(hits[,1]),] )

  if (use.mass2adduct) {
    if (!requireNamespace("mass2adduct",quietly = T))
    {
      stop("Package \"mass2adduct\" is required for this function")
    }
    diffs <- mass2adduct::adducts
    diffs[["diff"]] <- abs(diffs$mass - abs(data[var1,mz_col] - data[var2,mz_col]))
    diffs <- subset(diffs, diff <= tol)
    print(diffs)
  }
  return(diffm[diffm[["Accuracy(Da)"]] < tol,])
}

#'Find features by presence of MS2 peaks
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param peaks numeric or vector of numeric values. The m/z values of the MS2 peaks to search
#' @param operator c(all, any) operator to use with multiple given peaks. Defaults to all.
#' @param tol numeric. m/z tolerance for peak matching
#' @param lim numeric. Relative abundance cutoff of peaks
#' @param ms2_col regular expression. Name of the MS2 spectrum column
#'
#' @return a vector of hits (feature IDs)
#'
#' @export
find_by_MS2peaks <- function(data, peaks, operator=all, tol=0.005, lim=0.03, ms2_col="MS.*MS.*spectrum") {
  ms2_col <- grep(ms2_col, colnames(data), ignore.case = T, value = T)[1]
  hits <- c()
  for (r in rownames(data)) {
    spectrum <- data[r,ms2_col]
    if (is.na(spectrum)) next
    spectrum <- strsplit(spectrum,split = " ")[[1]]
    spectrum <- mapply(spectrum, FUN=strsplit, MoreArgs = list(split=":"))
    spectrum <- matrix(apply(as.data.frame(spectrum), MARGIN=1, FUN=as.numeric), ncol=2)
    spectrum <- as.data.frame(cbind(spectrum, spectrum[,2]/max(spectrum[,2])))
    colnames(spectrum) <- c("mz","ab","rel")
    spectrum <- spectrum[spectrum$rel>lim,]
    peaks_hits <- mapply(peaks, FUN=function(p) any(abs(p - spectrum$mz) < tol))
    if (operator(peaks_hits)) hits <- c(hits, r)
  }
  return(hits)
}

#'Find features by molecular weight or molecular formula
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param query numeric/string. Either the molecular mass or the molecular formula for adduct search
#' @param tol numeric. m/z tolerance for peak matching
#' @param mz_col regular expression. Name of the m/z containing column
#'
#' @return matrix. feature_IDs, matched_adducts and accuracy
#'
#' @export
find_by_mw <- function(data, query, tol=0.005, mz_col=".*mz") {
  mz_col <- grep(mz_col, colnames(data), ignore.case = T, value = T)[1]
  if (!is.numeric(query)) {
    query = find_mwfromformula(query)
  }
  qquery = find_adducts(query)

  hits <- c()
  for (q in rownames(qquery)) {
    diff <- abs(data[[mz_col]] - qquery[q,"mz"])
    newhits <- rownames(data)[diff < tol]
    if (length(newhits) > 0) {
      hits <- rbind(hits, cbind(newhits, diff[diff < tol], q))
    }
  }
  colnames(hits) <- c("Feature_ID","accuracy_Da","adduct_match")
  return(hits)
}


#' Find features by m/z value
#'
#' @param data a data.frame object. The feature data of the MS peaks
#' @param mz numeric. The m/z value for peak matching
#' @param tol numeric. m/z tolerance for peak matching
#' @param mz_col regular expression. Name of the m/z containing column
#'
#' @returns a vector of hits (feature IDs)
#' @export
#'
find_by_mz <- function(data, mz, tol=0.005, mz_col=".*mz") {
  mz_col <- grep(mz_col, colnames(data), ignore.case = T, value = T)[1]
  hits <- rownames(data)[abs(data[[mz_col]] - mz) < tol]
  return(hits)
}
