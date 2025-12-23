
#' Estimate proportion of variance of liability explained by SNP in general population
#'
#' This uses equation 10 in Lee et al. A Better Coefficient of Determination for Genetic Profile Analysis. 
#' Genetic Epidemiology 36: 214–224 (2012) <https://doi.org/10.1002/gepi.21614>.
#'
#' @param lor Vector of Log odds ratio.
#' @param af Vector of allele frequencies.
#' @param ncase Vector of Number of cases.
#' @param ncontrol Vector of Number of controls.
#' @param prevalence Vector of Disease prevalence in the population.
#' @param model Is the effect size estimated from the `"logit"` (default) or `"probit"` model.
#' @param correction Scale the estimated r by correction value. The default is `FALSE`.
#'
#' @export
#' @return Vector of signed r values
get_r_from_lor <- function(lor, af, ncase, ncontrol, prevalence, model="logit", correction=FALSE)
{
  
  # print(lor)
  # print(af)
  # print(ncase)
  # print(ncontrol)
  # print(prevalence)
  
  stopifnot(length(lor) == length(af))
  stopifnot(length(ncase) == 1 | length(ncase) == length(lor))
  stopifnot(length(ncontrol) == 1 | length(ncontrol) == length(lor))
  stopifnot(length(prevalence) == 1 | length(prevalence) == length(lor))
  if(length(prevalence) == 1 & length(lor) != 1)
  {
    prevalence <- rep(prevalence, length(lor))
  }
  if(length(ncase) == 1 & length(lor) != 1)
  {
    ncase <- rep(ncase, length(lor))
  }
  if(length(ncontrol) == 1 & length(lor) != 1)
  {
    ncontrol <- rep(ncontrol, length(lor))
  }
  
  nsnp <- length(lor)
  r <- array(NA, nsnp)
  for(i in 1:nsnp)
  {
    if(model == "logit")
    {
      ve <- pi^2/3
    } else if(model == "probit") {
      ve <- 1
    } else {
      stop("Model must be probit or logit")
    }
    popaf <- get_population_allele_frequency(af[i], ncase[i] / (ncase[i] + ncontrol[i]), exp(lor[i]), prevalence[i])
    vg <- lor[i]^2 * popaf * (1-popaf)
    r[i] <- vg / (vg + ve)
    if(correction)
    {
      r[i] <- r[i] / 0.58
    }
    r[i] <- sqrt(r[i]) * sign(lor[i])
  }
  return(r)
}

#' Obtain 2x2 contingency table from marginal parameters and odds ratio
#'
#' Columns are the case and control frequencies.
#' Rows are the frequencies for allele 1 and allele 2.
#'
#' @param af Allele frequency of effect allele.
#' @param prop Proportion of cases.
#' @param odds_ratio Odds ratio.
#' @param eps tolerance, default is `1e-15`.
#'
#' @export
#' @return 2x2 contingency table as matrix
contingency <- function(af, prop, odds_ratio, eps=1e-15)
{
  a <- odds_ratio-1
  b <- (af+prop)*(1-odds_ratio)-1
  c_ <- odds_ratio*af*prop
  
  if (abs(a) < eps)
  {
    z <- -c_ / b
  } else {
    d <- b^2 - 4*a*c_
    if (d < eps*eps) 
    {
      s <- 0
    } else {
      s <- c(-1,1)
    }
    z <- (-b + s*sqrt(max(0, d))) / (2*a)
  }
  y <- vapply(z, function(a) zapsmall(matrix(c(a, prop-a, af-a, 1+a-af-prop), 2, 2)), matrix(0.0, 2, 2))
  i <- apply(y, 3, function(u) all(u >= 0))
  return(y[,,i])
}

#' Estimate allele frequency from SNP
#'
#' @param g Vector of 0/1/2
#'
#' @export
#' @return Allele frequency 
allele_frequency <- function(g)
{
  (sum(g == 1) + 2 * sum(g == 2)) / (2 * sum(!is.na(g)))
}


#' Estimate the allele frequency in population from case/control summary data
#'
#' @param af Effect allele frequency (or MAF)
#' @param prop Proportion of samples that are cases
#' @param odds_ratio Odds ratio
#' @param prevalence Population disease prevalence
#'
#' @export
#' @return Population allele frequency
get_population_allele_frequency <- function(af, prop, odds_ratio, prevalence)
{
  stopifnot(length(af) == length(odds_ratio))
  stopifnot(length(prop) == length(odds_ratio))
  for(i in 1:length(odds_ratio))
  {
    co <- contingency(af[i], prop[i], odds_ratio[i])
    af_controls <- co[1,2] / (co[1,2] + co[2,2])
    af_cases <- co[1,1] / (co[1,1] + co[2,1])
    af[i] <- af_controls * (1 - prevalence) + af_cases * prevalence
  }
  return(af)
}
