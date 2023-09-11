#' Parametric Two Factor Fixed Effects ANOVA for censored data
#'
#' @description Computes tests for each of the two factors and optionally for their interaction using likelihood ratio tests. p-values will not be identical to the usual method of moments ANOVA tests but will be similar.
#' @param y1 The column of data values plus detection limits
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y1` column, and 0 (or `FALSE`) indicates a detected value in `y1`.
#' @param fac1 The first grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param fac2 The second grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param LOG A logical variable indicating whether natural logs are to be taken of the 'y1' column data. Default is TRUE.
#' @param interact A logical variable indicating whether to compute an interaction term between the two variables. Default is TRUE.
#' #' @keywords two-way two-factor factorial ANOVA analysis of variance censored
#' @return Q-Q plots of residuals. Likelihood ratio test statistics ("chisquare"), degrees of freedom ("df") and p-values (pval) for two factors and optionally the interaction. Data on the underlying models, including AIC and R2 are also provided.
#' @export
#' @details Tests are computed using Maximum Likelihood Estimation. When a gaussian	distribution model is	used (LOG=FALSE) modeled values may fall below zero, producing unreal p-values (often lower than they should be).  Because of this, testing in log units is preferable and is the default unless you are transforming the y values prior to running the function (such as taking cube roots to approximate a gamma distribution).  The 'delta.lr0x2' stat output is the -2loglikehood for the test of the model versus an intercept-only model.
#' @seealso [survival::survreg]
#' @importFrom EnvStats qqPlotCensored gofTestCensored qqPlot gofTest
#' @importFrom survival survreg Surv
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#' Millard, S.P., 2013. EnvStats: An R Package for Environmental Statistics. Springer-Verlag, New York.
#'
#'@examples
#' data(Gales_Creek)
#' Gales_Creek$Period <- c(rep("early", 35), rep("middle", 12), rep("late", 16))
#' with(Gales_Creek,cen2way(TCr, CrND, Season, Period))
#'


cen2way <- function(y1, y2, fac1, fac2, LOG = TRUE, interact = TRUE) {
  yname <- deparse(substitute(y1))
  name.fac1 <- deparse(substitute(fac1))
  name.fac2 <- deparse(substitute(fac2))

  # for both log and original units
  Fact1 <- as.factor(fac1)
  Fact2 <- as.factor(fac2)
  facts <- data.frame (Fact1, Fact2)
  nonas <- na.omit(cbind(y1, y2, facts))     # for nonas[,2] 1 = censored, 0 = detected
  n.rows <- nrow(nonas)   # number of observations
  xnona <- nonas[,-(1:2)]
  AIC <- matrix(0, nrow = 4, ncol = 1)
  BIC <- matrix(0, nrow = 4, ncol =1)
  Rescaled_Loglik_R2 <- matrix(0, nrow = 4, ncol =1)

  # create factor design matrices e, d, int
    # factor 1
     levels.1 <- levels(xnona[,1])  # names of levels in factor 1
      nlev.1  <- nlevels(xnona[,1])  # number of levels in factor 1
      ncol.1 <- nlev.1 - 1    # number of columns in matrix for factor 1 -> j
      fac1.names <- matrix("a", nrow = 1, ncol= ncol.1)
      e <- matrix(0, nrow = n.rows, ncol = ncol.1)
       for (j in 1:ncol.1) { fac1.names[j] <- paste(name.fac1, as.character(j), sep = "")
           for (i in 1:n.rows) {
           if (xnona[i,1] == levels.1[nlev.1]) {e[i,j]<- -1}
           if (xnona[i,1] == levels.1[j]) {e[i,j]<- 1}
    } }
  # factor 2
   levels.2 <- levels(xnona[,2])  # names of levels in factor 2
    nlev.2  <- nlevels(xnona[,2])  # number of levels in factor 2
    ncol.2 <- nlev.2 - 1    # number of columns in matrix for factor 2 -> k
    fac2.names <- matrix("a", nrow = 1, ncol= ncol.2)
    d <- matrix(0, nrow = n.rows, ncol = ncol.2)
        for (k in 1:ncol.2) { fac2.names[k] <- paste(name.fac2, as.character(k), sep = "")
    for (i in 1:n.rows) {
      if (xnona[i,2] == levels.2[nlev.2]) {d[i,k]<- -1}
      if (xnona[i,2] == levels.2[k]) {d[i,k]<- 1}
    } }
    colnames(e) <- fac1.names
#  print (e)
    colnames(d) <- fac2.names
#  print (d)

    # interaction terms
    i <- 0
    ncol.int <- ncol.1 * ncol.2
    int <- matrix(0, nrow = n.rows, ncol = ncol.int)
    int.name <- matrix("a", nrow = 1, ncol= ncol.int)
      for (j in 1:ncol.1) {
      for (k in 1:ncol.2) {
        i <- i+1
        int[,i] <- e[,j] * d[,k]
        int.name[i] <- paste("int", as.character(i), sep = "")
      }
    }
    colnames(int) <- int.name
#    print (int)
    x.all <- data.frame (e, d, int)
#    print (x.all)

# REGRESSION
# ln units for LOG = TRUE
 if (LOG == TRUE) {
 lnvar <- log(nonas[,1])    # Y in log units (default)
flip.log <- max(lnvar) +1 - lnvar
surv.log <- Surv(flip.log, as.logical(1-nonas[,2]) )
vtext<- paste("Quantiles of", yname, "residuals (log units)")

#  reg.out is model with both factors and interaction
  reg.out <- survreg(surv.log ~ ., data = x.all, dist = "gaussian")
  cn <- names(reg.out$coefficients[-1])
  xvars.txt <- cn[1]
  for (i in 1:length(cn))  {j <-(i+1)
  if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
  reg.out$call[3] <- xvars.txt
  }
ylog.pred <- max(lnvar) +1 - reg.out$linear.predictors
reg.out$call[2] <- paste("log(", yname, ")", sep = "")
reg.out$coefficients <- reg.out$coefficients * (-1)
reg.out$coefficients[1] <- max(lnvar)+1 + reg.out$coefficients[1]   # coeffs in cenreg lognormal

ylog.resi <- lnvar - ylog.pred
reg.out$linear.predictors <- ylog.pred
reg.out$resids <- ylog.resi

Rescaled_Loglik_R2[4] <- signif((1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))) /(1-exp(2*reg.out$loglik[1]/length(reg.out$y))),4)
AIC[4] <- -2*reg.out$loglik[2] + (2*reg.out$df +1)
BIC[4] <- -2*reg.out$loglik[2] + log(reg.out$df+reg.out$df.residual+1)*reg.out$df

# both factors but no interaction model. Used to test for interaction effect. Also basis for no interaction factor tests.
x.noint <- data.frame (e, d)
# print (x.noint)
reg.noint <- survreg(surv.log ~ ., data = x.noint, dist = "gaussian")
cn <- names(reg.noint$coefficients[-1])
xvars.txt <- cn[1]
for (i in 1:length(cn))  {j <-(i+1)
if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
reg.noint$call[3] <- xvars.txt
}
noint.pred <- max(lnvar) +1 - reg.noint$linear.predictors
reg.noint$call[2] <- paste("log(", yname, ")", sep = "")
reg.noint$coefficients <- reg.noint$coefficients * (-1)
reg.noint$coefficients[1] <- max(lnvar)+1 + reg.noint$coefficients[1]   # coeffs in cenreg lognormal

noint.resi <- lnvar - noint.pred
reg.noint$linear.predictors <- noint.pred
reg.noint$resids <- noint.resi

Rescaled_Loglik_R2[3] <- signif((1-exp(-2*(reg.noint$loglik[2]-reg.noint$loglik[1])/length(reg.noint$y))) /(1-exp(2*reg.noint$loglik[1]/length(reg.noint$y))),4)
AIC[3] <- -2*reg.noint$loglik[2] + (2*reg.noint$df +1)
BIC[3] <- -2*reg.noint$loglik[2] + log(reg.noint$df+reg.noint$df.residual+1)*reg.noint$df


if (interact == TRUE) {
  header <- c("       ANOVA Table with interaction")
  if (sum(nonas[,2]) == 0) {testnorm <- gofTest(ylog.resi, distribution = "norm")}
    else {testnorm <- gofTestCensored(ylog.resi,nonas[,2])}
  ptext <- paste(name.fac1,"+",name.fac2,"+ interaction model","\n", "Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
  if (sum(nonas[,2]) == 0) {qqPlot(ylog.resi, add.line = T, ylab = vtext, main = "Lognormal Q-Q Plot of residuals")}
    else {qqPlotCensored(ylog.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Lognormal Q-Q Plot of residuals")}
  mtext(ptext, cex = 0.7)

# factor 1 but not factor 2. Used to test for factor 2 effect
    x.fac1 <- data.frame (e, int)
    reg.fac1 <- survreg(surv.log ~ ., data = x.fac1, dist = "gaussian")
    cn <- names(reg.fac1$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    reg.fac1$call[3] <- xvars.txt
    }
    fac1.pred <- max(lnvar) +1 - reg.fac1$linear.predictors
    reg.fac1$call[2] <- paste("log(", yname, ")", sep = "")
    reg.fac1$coefficients <- reg.fac1$coefficients * (-1)
    reg.fac1$coefficients[1] <- max(lnvar)+1 + reg.fac1$coefficients[1]   # coeffs in cenreg lognormal
    fac1.resi <- lnvar - fac1.pred
    reg.fac1$linear.predictors <- fac1.pred
    reg.fac1$resids <- fac1.resi

    Rescaled_Loglik_R2[1] <- signif((1-exp(-2*(reg.fac1$loglik[2]-reg.fac1$loglik[1])/length(reg.fac1$y))) /(1-exp(2*reg.fac1$loglik[1]/length(reg.fac1$y))),4)
    AIC[1] <- -2*reg.fac1$loglik[2] + (2*reg.fac1$df +1)
    BIC[1] <- -2*reg.fac1$loglik[2] + log(reg.fac1$df+reg.fac1$df.residual+1)*reg.fac1$df

# factor 2 but not factor 1. Used to test for factor 1 effect
    x.fac2 <- data.frame (d, int)
    reg.fac2 <- survreg(surv.log ~ ., data = x.fac2, dist = "gaussian")
    cn <- names(reg.fac2$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    reg.fac2$call[3] <- xvars.txt
    }
    fac2.pred <- max(lnvar) +1 - reg.fac2$linear.predictors
    reg.fac2$call[2] <- paste("log(", yname, ")", sep = "")
    reg.fac2$coefficients <- reg.fac2$coefficients * (-1)
    reg.fac2$coefficients[1] <- max(lnvar)+1 + reg.fac2$coefficients[1]   # coeffs in cenreg lognormal
    fac2.resi <- lnvar - fac2.pred
    reg.fac2$linear.predictors <- fac2.pred
    reg.fac2$resids <- fac2.resi

    Rescaled_Loglik_R2[2] <- signif((1-exp(-2*(reg.fac2$loglik[2]-reg.fac2$loglik[1])/length(reg.fac2$y))) /(1-exp(2*reg.fac2$loglik[1]/length(reg.fac2$y))),4)
    AIC[2] <- -2*reg.fac2$loglik[2] + (2*reg.fac2$df +1)
    BIC[2] <- -2*reg.fac2$loglik[2] + log(reg.fac2$df+reg.fac2$df.residual+1)*reg.fac2$df

# run the tests with interaction
    chisquare <- c(signif(2*(reg.out$loglik[2]-reg.fac2$loglik[2]), 4), signif(2*(reg.out$loglik[2]-reg.fac1$loglik[2]),4), signif(2*(reg.out$loglik[2]-reg.noint$loglik[2]),4))
 }

if (interact == FALSE) {
  header <- c("       ANOVA Table without interaction")
  if (sum(nonas[,2]) == 0) {testnorm <- gofTest(noint.resi, distribution = "norm")}
    else {testnorm <- gofTestCensored(noint.resi,nonas[,2])}
  ptext <- paste(name.fac1,"+",name.fac2,"factors","\n", "Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
  if (sum(nonas[,2]) == 0) {qqPlot(noint.resi, add.line = T, ylab = vtext, main = "Lognormal Q-Q Plot of residuals")}
    else {qqPlotCensored(noint.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Lognormal Q-Q Plot of residuals")}
  mtext(ptext, cex = 0.7)

  # factor 2 only, no interaction, use to test for factor 1
  x.fac2 <- data.frame (d)
  # print (x.fac2)
  reg.fac2 <- survreg(surv.log ~ ., data = x.fac2, dist = "gaussian")
  cn <- names(reg.fac2$coefficients[-1])
  xvars.txt <- cn[1]
  for (i in 1:length(cn))  {j <-(i+1)
  if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
  reg.fac2$call[3] <- xvars.txt
  }
  ylog.pred <- max(lnvar) +1 - reg.fac2$linear.predictors
  reg.fac2$call[2] <- paste("log(", yname, ")", sep = "")
  reg.fac2$coefficients <- reg.fac2$coefficients * (-1)
  reg.fac2$coefficients[1] <- max(lnvar)+1 + reg.fac2$coefficients[1]   # coeffs in cenreg lognormal
  ylog.resi <- lnvar - ylog.pred
  reg.fac2$linear.predictors <- ylog.pred
  reg.fac2$resids <- ylog.resi

  Rescaled_Loglik_R2[2] <- signif((1-exp(-2*(reg.fac2$loglik[2]-reg.fac2$loglik[1])/length(reg.fac2$y))) /(1-exp(2*reg.fac2$loglik[1]/length(reg.fac2$y))),4)
  AIC[2] <- -2*reg.fac2$loglik[2] + (2*reg.fac2$df +1)
  BIC[2] <- -2*reg.fac2$loglik[2] + log(reg.fac2$df+reg.fac2$df.residual+1)*reg.fac2$df

  # factor 1 only, no interaction, use to test for factor 2
  x.fac1 <- data.frame (e)
  # print (x.fac1)
  reg.fac1 <- survreg(surv.log ~ ., data = x.fac1, dist = "gaussian")
  cn <- names(reg.fac1$coefficients[-1])
  xvars.txt <- cn[1]
  for (i in 1:length(cn))  {j <-(i+1)
  if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
  reg.fac1$call[3] <- xvars.txt
  }
  ylog.pred <- max(lnvar) +1 - reg.fac1$linear.predictors
  reg.fac1$call[2] <- paste("log(", yname, ")", sep = "")
  reg.fac1$coefficients <- reg.fac1$coefficients * (-1)
  reg.fac1$coefficients[1] <- max(lnvar)+1 + reg.fac1$coefficients[1]   # coeffs in cenreg lognormal
  ylog.resi <- lnvar - ylog.pred
  reg.fac1$linear.predictors <- ylog.pred
  reg.fac1$resids <- ylog.resi

  Rescaled_Loglik_R2[1] <- signif((1-exp(-2*(reg.fac1$loglik[2]-reg.fac1$loglik[1])/length(reg.fac1$y))) /(1-exp(2*reg.fac1$loglik[1]/length(reg.fac1$y))),4)
  AIC[1] <- -2*reg.fac1$loglik[2] + (2*reg.fac1$df +1)
  BIC[1] <- -2*reg.fac1$loglik[2] + log(reg.fac1$df+reg.fac1$df.residual+1)*reg.fac1$df
  chisquare <- c(signif(2*(reg.noint$loglik[2]-reg.fac2$loglik[2]), 4), signif(2*(reg.noint$loglik[2]-reg.fac1$loglik[2]),4), signif(2*(reg.out$loglik[2]-reg.noint$loglik[2]),4))
}  # end of no interaction tests

  cat("Two-way Fixed Effects ANOVA by MLE Likelihood Ratio Tests for Censored Data", "\n")
  cat("ln(",yname,")", " modeled by Factor 1: ", name.fac1, "  and Factor 2: ", name.fac2, "\n", sep="")
 }   # end of log units

# y in original units
    else{
      y.low <- nonas[,1]*(1-nonas[,2])                  #  0 for low end of all NDs
      surv.norm <- Surv(y.low, nonas[,1], type="interval2")
      vtext<- paste("Quantiles of", yname, "residuals")

      # reg.out.  2 factors plus interaction model
        reg.out <- survreg(surv.norm ~ ., data = x.all, dist = "gaussian")
        cn <- names(reg.out$coefficients[-1])
        xvars.txt <- cn[1]
        for (i in 1:length(cn))  {j <-(i+1)
        if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
        }
        reg.out$call[3] <- xvars.txt
        reg.out$call[2] <- yname
        ynorm.pred <- reg.out$linear.predictors
        ynorm.resi <- as.numeric(y.low - ynorm.pred)
        reg.out$linear.predictors <- ynorm.pred
        reg.out$resids <- ynorm.resi

        Rescaled_Loglik_R2[4] <- signif((1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))) /(1-exp(2*reg.out$loglik[1]/length(reg.out$y))),4)
        AIC[4] <- -2*reg.out$loglik[2] + (2*reg.out$df +1)
        BIC[4] <- -2*reg.out$loglik[2] + log(reg.out$df+reg.out$df.residual+1)*reg.out$df

        # no interaction model used to test both interaction and noint
        x.noint <- data.frame (e, d)

        # factors 1 and 2, no interaction
        reg.noint <- survreg(surv.norm ~ ., data = x.noint, dist = "gaussian")
        cn <- names(reg.noint$coefficients[-1])
        xvars.txt <- cn[1]
        for (i in 1:length(cn))  {j <-(i+1)
        if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
        }
        reg.noint$call[3] <- xvars.txt
        ynorm.pred <- reg.noint$linear.predictors
        reg.noint$call[2] <- yname
        ynorm.resi <- y.low - ynorm.pred
        reg.noint$linear.predictors <- ynorm.pred
        reg.noint$resids <- ynorm.resi

        Rescaled_Loglik_R2[3] <- signif((1-exp(-2*(reg.noint$loglik[2]-reg.noint$loglik[1])/length(reg.noint$y))) /(1-exp(2*reg.noint$loglik[1]/length(reg.noint$y))),4)
        AIC[3] <- -2*reg.noint$loglik[2] + (2*reg.noint$df +1)
        BIC[3] <- -2*reg.noint$loglik[2] + log(reg.noint$df+reg.noint$df.residual+1)*reg.noint$df

        if (interact == TRUE) {
        header <- c("       ANOVA Table with interaction")
        if (sum(nonas[,2]) == 0) {testnorm <- gofTest(ynorm.resi, distribution = "norm")}
        else {testnorm <- gofTestCensored(ynorm.resi,nonas[,2])}
        ptext <- paste(name.fac1,"+",name.fac2,"+ interaction model","\n", "Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
        if (sum(nonas[,2]) == 0) {qqPlot(ynorm.resi, add.line = T, ylab = vtext, main = "Normal Q-Q Plot of residuals")}
          else {qqPlotCensored(ynorm.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Normal Q-Q Plot of residuals")}
        mtext(ptext, cex = 0.7)

        # factor 1 but not factor 2. Used to test for factor 2 effect
        x.fac1 <- data.frame (e, int)
        reg.fac1 <- survreg(surv.norm ~ ., data = x.fac1, dist = "gaussian")
        cn <- names(reg.fac1$coefficients[-1])
        xvars.txt <- cn[1]
        for (i in 1:length(cn))  {j <-(i+1)
        if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
        reg.fac1$call[3] <- xvars.txt
        }
        reg.fac1$call[3] <- xvars.txt
        reg.fac1$call[2] <- yname
        ynorm.pred <- reg.fac1$linear.predictors
        ynorm.resi <- as.numeric(y.low - ynorm.pred)
        reg.fac1$linear.predictors <- ynorm.pred
        reg.fac1$resids <- ynorm.resi

        Rescaled_Loglik_R2[1] <- signif((1-exp(-2*(reg.fac1$loglik[2]-reg.fac1$loglik[1])/length(reg.fac1$y))) /(1-exp(2*reg.fac1$loglik[1]/length(reg.fac1$y))),4)
        AIC[1] <- -2*reg.fac1$loglik[2] + (2*reg.fac1$df +1)
        BIC[1] <- -2*reg.fac1$loglik[2] + log(reg.fac1$df+reg.fac1$df.residual+1)*reg.fac1$df

        # factor 2 but not factor 1. Used to test for factor 1 effect
        x.fac2 <- data.frame (d, int)
        reg.fac2 <- survreg(surv.norm ~ ., data = x.fac2, dist = "gaussian")
        cn <- names(reg.fac2$coefficients[-1])
        xvars.txt <- cn[1]
        for (i in 1:length(cn))  {j <-(i+1)
        if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
        reg.fac2$call[3] <- xvars.txt
        }
        reg.fac2$call[2] <- yname
        ynorm.pred <- reg.fac2$linear.predictors
        ynorm.resi <- as.numeric(y.low - ynorm.pred)
        reg.fac2$linear.predictors <- ynorm.pred
        reg.fac2$resids <- ynorm.resi

        Rescaled_Loglik_R2[2] <- signif((1-exp(-2*(reg.fac2$loglik[2]-reg.fac2$loglik[1])/length(reg.fac2$y))) /(1-exp(2*reg.fac2$loglik[1]/length(reg.fac2$y))),4)
        AIC[2] <- -2*reg.fac2$loglik[2] + (2*reg.fac2$df +1)
        BIC[2] <- -2*reg.fac2$loglik[2] + log(reg.fac2$df+reg.fac2$df.residual+1)*reg.fac2$df

        # run the tests with interaction
        chisquare <- c(signif(2*(reg.out$loglik[2]-reg.fac2$loglik[2]), 4), signif(2*(reg.out$loglik[2]-reg.fac1$loglik[2]),4), signif(2*(reg.out$loglik[2]-reg.noint$loglik[2]),4))
          }  # end if with interaction

        if (interact == FALSE)  {
        header <- c("       ANOVA Table without interaction")
        if (sum(nonas[,2]) == 0) {testnorm <- gofTest(ynorm.resi, distribution = "norm")}
          else {testnorm <- gofTestCensored(ynorm.resi,nonas[,2])}
        ptext <- paste(name.fac1,"+",name.fac2,"factors","\n", "Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
        if (sum(nonas[,2]) == 0) {qqPlot(ynorm.resi, add.line = T, ylab = vtext, main = "Normal Q-Q Plot of residuals")}
          else {qqPlotCensored(ynorm.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Normal Q-Q Plot of residuals")}
        mtext(ptext, cex = 0.7)

      # factor 2 only, no interaction, use to test for factor 1
      x.fac2 <- data.frame (d)
      reg.fac2 <- survreg(surv.norm ~ ., data = x.fac2, dist = "gaussian")
      cn <- names(reg.fac2$coefficients[-1])
      xvars.txt <- cn[1]
      for (i in 1:length(cn))  {j <-(i+1)
      if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
      }

      reg.fac2$call[3] <- xvars.txt
      reg.fac2$call[2] <- yname
      ynorm.pred <- reg.fac2$linear.predictors
      ynorm.resi <- y.low - ynorm.pred
      reg.fac2$linear.predictors <- ynorm.pred
      reg.fac2$resids <- ynorm.resi

       Rescaled_Loglik_R2[ 2] <- signif((1-exp(-2*(reg.fac2$loglik[2]-reg.fac2$loglik[1])/length(reg.fac2$y))) /(1-exp(2*reg.fac2$loglik[1]/length(reg.fac2$y))),4)
      AIC[2] <- -2*reg.fac2$loglik[2] + (2*reg.fac2$df +1)
      BIC[2] <- -2*reg.fac2$loglik[2] + log(reg.fac2$df+reg.fac2$df.residual+1)*reg.fac2$df

      # factor 1 only, no interaction, use to test for factor 2
      x.fac1 <- data.frame (e)
      # print (x.fac1)
      reg.fac1 <- survreg(surv.norm ~ ., data = x.fac1, dist = "gaussian")
      cn <- names(reg.fac1$coefficients[-1])
      xvars.txt <- cn[1]
      for (i in 1:length(cn))  {j <-(i+1)
      if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
      }
      reg.fac1$call[3] <- xvars.txt
      reg.fac1$call[2] <- yname
      ynorm.pred <- reg.fac1$linear.predictors
      ynorm.resi <- y.low - ynorm.pred
      reg.fac1$linear.predictors <- ynorm.pred
      reg.fac1$resids <- ynorm.resi

    Rescaled_Loglik_R2[1] <- signif((1-exp(-2*(reg.fac1$loglik[2]-reg.fac1$loglik[1])/length(reg.fac1$y))) /(1-exp(2*reg.fac1$loglik[1]/length(reg.fac1$y))),4)
    AIC[1] <- -2*reg.fac1$loglik[2] + (2*reg.fac1$df +1)
    BIC[1] <- -2*reg.fac1$loglik[2] + log(reg.fac1$df+reg.fac1$df.residual+1)*reg.fac1$df
    chisquare <- c(signif(2*(reg.noint$loglik[2]-reg.fac2$loglik[2]), 4), signif(2*(reg.noint$loglik[2]-reg.fac1$loglik[2]),4), signif(2*(reg.out$loglik[2]-reg.noint$loglik[2]),4))
        }  # end of no interaction tests

        cat("Two-way Fixed Effects ANOVA by MLE Likelihood Ratio Tests for Censored Data", "\n")
        cat(yname, " modeled by Factor 1: ", name.fac1, "  and Factor 2: ", name.fac2, "\n", sep="")
    }  # end of original units

    # final printing
        cat(header,"\n")
     FACTORS <- c(name.fac1, name.fac2, "interaction")
    df <- c(ncol.1, ncol.2, ncol.int)
    pv1 <- signif(pchisq(chisquare[1], df[1], lower.tail = FALSE),4)
    pv2 <- signif(pchisq(chisquare[2], df[2], lower.tail = FALSE),4)
    pv3 <- signif(pchisq(chisquare[3], df[3], lower.tail = FALSE),4)
    pval <- c(pv1, pv2, pv3)
   anova.table <- data.frame(FACTORS, chisquare, df, pval)
   if (interact == FALSE) {anova.table <- anova.table[-3,]}
   print(anova.table, print.gap = 3, row.names = FALSE)
   cat("\n")

   loglikelihood <- c(reg.fac1$loglik[2], reg.fac2$loglik[2], reg.noint$loglik[2], reg.out$loglik[2])
   delta.lr0x2 <- c(signif(-2*(reg.fac1$loglik[1]-reg.fac1$loglik[2]), 4), signif(-2*(reg.fac2$loglik[1]-reg.fac2$loglik[2]),4), signif(-2*(reg.noint$loglik[1]-reg.noint$loglik[2]),4), signif(-2*(reg.out$loglik[1]-reg.out$loglik[2]),4))
   MODELS <- c(paste(name.fac1,"only"), paste(name.fac2,"only"), "both factors", "both + interaction")
    the.end <- data.frame (MODELS, loglikelihood, delta.lr0x2, AIC, BIC, Rescaled_Loglik_R2)
    if (interact == FALSE) {the.end <- the.end[-4,]}
    print(the.end, print.gap = 3, row.names = FALSE)
}
