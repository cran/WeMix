require(WeMix)
require(testthat)

# for sleepstudy, calls in here to lmer
library(lme4)
# for grad function
library(numDeriv)

options(width = 500)
options(useFancyQuotes = FALSE)

### Data Used: sleepstudy
tolerance <- 2E-7
data("sleepstudy")
### Unweigted =====================================
sleepstudyU <- sleepstudy
sleepstudyU$weight1L1 <- 1
sleepstudyU$weight1L2 <- 1

context("model runs")
test_that("lnl agrees with lme4 1", {
  
  system.time(wm0 <- mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), nQuad=13, verbose=FALSE, fast=TRUE,run=TRUE))
  wm1 <- mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), nQuad=13, verbose=FALSE, run=FALSE)
  expect_equal(wm0$lnl, -897.039321502613, tolerance=tolerance*897) # value from  lmer(Reaction ~ Days + (1 | Subject), data=sleepstudy, REML=FALSE)
  
  expect_equal(unname(c(wm0$coef, wm0$vars)), unname(wm1$parlme), tolerance=tolerance*wm1$parlme)
  # agrees with lme4
  # check coef
  lme1 <- lmer(Reaction ~ Days + (1|Subject), data=sleepstudy, REML=FALSE)
  expect_equal(coef(wm0),
               expected = getME(lme1, "fixef"),
               tolerance = tolerance * coef(wm0))  
  # check vars
  lmevars1 <- data.frame(summary(lme1)$varcor)$sdcor
  expect_equal(unname(wm0$vars),
               expected = unname(lmevars1)^2,
               tolerance = tolerance * wm0$vars)  
  # agrees with GLLAMM
  #source: logit_command.do
  gllamm_model1 <- c("(Intercept)" = 251.4051,
                     "Days"         = 10.46729,
                     "Residual"     = 1296.8768,
                     "Subject"      = 954.52789,
                     "lnl"          = -897.03932)
  expect_equal(unname(coef(wm0)),
               expected = unname(gllamm_model1[1:2]),
               tolerance = abs(tolerance * gllamm_model1[1:2]))
  expect_equal(unname(wm0$vars),
               expected = unname(gllamm_model1[3:4]),
               tolerance = tolerance * wm0$vars)
  expect_equal(unname(wm0$lnl),
               expected=unname(gllamm_model1[5]),
               tolerance=abs(tolerance*wm0$lnl))
})

context("unweighted model with 1 random effect")
test_that("lnl agrees with lme4 2", {
  wm1 <- mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), nQuad=13, verbose=FALSE, run=FALSE)
  system.time(wm1lnl <- wm1$lnlf(wm1$parlme)) # speed test 1, for a longer test, increase nQuad on prev line
  expect_equal(wm1lnl, -897.039321502613, tolerance=tolerance*897) # value from  lmer(Reaction ~ Days + (1 | Subject), data=sleepstudy, REML=FALSE)
})

context("unweighted model with 2 random effects")
test_that("agrees with lme4 3", {
  lme2 <- lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), data=sleepstudyU, REML=FALSE)
  wm2 <- mix(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), nQuad=13L,fast=TRUE, verbose=FALSE, run=FALSE)
  expect_equal(wm2$lnlf(wm2$parlme), as.numeric(logLik(lme2)), tol=1E-7*abs(as.numeric(logLik(lme2))))
  system.time(grad <- grad(wm2$lnlf,wm2$parlme))
  expect_equal(grad, rep(0, length(wm2$parlme)), tolerance=1E-5)
  
  #test actal values on data subset for speed
  ss_sample <- sleepstudyU[row.names(sleepstudy) %in% seq(1,60),]
  
  lme2Sample <- lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), data=ss_sample, REML=FALSE)
  wm2Sample <- mix(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), data=ss_sample, weights=c("weight1L1", "weight1L2"), fast=TRUE, nQuad=13, verbose=FALSE, run=TRUE)
  # check coef
  expect_equal(coef(wm2Sample),
               expected = getME(lme2Sample, "fixef"),
               tolerance = tolerance * coef(wm2Sample))
  # check vars
  lmewm2vars <- data.frame(summary(lme2Sample)$varcor)$sdcor
  expect_equal(unname(wm2Sample$vars),
               expected = unname(lmewm2vars)^2,
               tolerance = tolerance * wm2Sample$vars)
})

ss <- sleepstudy
ss1 <- ss
ss2 <- ss
doubles <- c(308, 309, 310) # subject with double obs
ss2 <- rbind(ss, subset(ss, Subject %in% doubles))

ss1$W1 <- ifelse(ss1$Subject %in% doubles, 2, 1)
ss1$W2 <- 1

ss2$W2 <- ss2$W1 <- 1

doubles <- c(308, 309, 310) # subject with double obs
ss30 <- subset(ss, Subject %in% doubles)
ss30$Subject <- as.numeric(as.character(ss30$Subject)) + 1000
ss0 <- ss
ss0$Subject <- as.numeric(as.character(ss$Subject))
ss3 <- rbind(ss0, ss30)
ss3$Subject <- as.factor(ss3$Subject)

ss3$W2 <- 1
ss3$W1 <- 1

ss4 <- ss
ss4$W2 <- ifelse(ss4$Subject %in% doubles, 2, 1)
ss4$W1 <- 1

context("repeating is the same as weighting: L1 replicate vs weigting")
test_that("L1 replicate vs weigting", {
  # mix for L1, weighted
  wmeL1W <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss1,
                weights=c("W1", "W2"), nQuad=13, run=FALSE, fast=TRUE, verbose=FALSE)
  
  # mix for L1, duplicated
  system.time(wmeL1D <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss2,
                            weights=c("W1", "W2"), nQuad=13, run=FALSE,fast=TRUE,  verbose=FALSE))
  #statares <- c(251.4619, 10.40726, 1000.7466, 1338.0865) # not used
  
  # check weighted agrees with duplicated lme4 results
  expect_equal(wmeL1W$lnlf(wmeL1D$parlme), -1048.34318418762, tolerance=1050*tolerance)
  grd <- numDeriv::grad(wmeL1W$lnlf, wmeL1D$parlme)
  expect_equal(grd, rep(0,length(wmeL1D$parlme)), tolerance=tolerance)
  
  # check duplicated agrees with duplicated lme4 results
  expect_equal(wmeL1D$lnlf(wmeL1D$parlme), -1048.34318418762, tolerance=1050*tolerance)
  grd <- numDeriv::grad(wmeL1D$lnlf, wmeL1D$parlme)
  expect_equal(grd, rep(0,length(wmeL1D$parlme)), tolerance=tolerance)
})

context("grouping factor not sorted")
test_that("grouping factor not sorted", {
  ss1_mixed <- ss1[c(125:180,1,100,2:99,101:124),]
  row.names(ss1_mixed) <- NULL
  # mix for L1, weighted
  wmeL1W <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss1_mixed,
                weights=c("W1", "W2"), nQuad=13, run=FALSE,fast=TRUE,  verbose=FALSE)
  
  # mix for L1, duplicated
  system.time(wmeL1D <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss2,
                            weights=c("W1", "W2"), nQuad=13, run=FALSE,fast=TRUE,  verbose=FALSE))
  #statares <- c(251.4619, 10.40726, 1000.7466, 1338.0865) # not used
  
  # check weighted agrees with duplicated lme4 results
  expect_equal(wmeL1W$lnlf(wmeL1D$parlme), -1048.34318418762, tolerance=1050*tolerance)
  grd <- numDeriv::grad(wmeL1W$lnlf, wmeL1D$parlme)
  expect_equal(grd, rep(0,length(wmeL1D$parlme)), tolerance=tolerance)
  
  # check final results
  mix1 <-  mix(formula=Reaction ~ Days + (1 | Subject), data=ss1_mixed,
               weights=c("W1", "W2"), nQuad=13, run=TRUE, verbose=FALSE, fast = TRUE)
  mix1REF <-  mix(formula=Reaction ~ Days + (1 | Subject), data=ss1,
                  weights=c("W1", "W2"), nQuad=13, run=TRUE, verbose=FALSE, fast = TRUE)
  expect_equal(mix1$coef, mix1REF$coef)
  expect_equal(mix1$vars, mix1REF$vars)
  expect_equal(mix1$lnl, mix1REF$lnl)
  
})

context("repeating is the same as weighting: L2 replicate vs weigting")
test_that("L2 replicate vs weigting", {
  # mix for L2, duplicated
  system.time(wmeL2D <- mix(formula=Reaction ~ Days + (1 | Subject),
                            data=ss3, weights=c("W1", "W2"),
                            nQuad=13, run=FALSE, verbose=FALSE))
  
  # mix for L2, weighted
  system.time(wmeL2W <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss4,
                            weights=c("W1", "W2"), nQuad=13, run=FALSE, verbose=FALSE))
  
  expect_equal(wmeL2W$lnlf(wmeL2D$parlme), -1055.34690957995, tolerance=1050*2E-7)
  grd <- numDeriv::grad(wmeL2W$lnlf, wmeL2D$parlme)
  expect_equal(grd, rep(0,length(grd)), tolerance=tolerance*100) # note larger tollerance
  expect_equal(wmeL2D$lnlf(wmeL2D$parlme), -1055.34690957995, tolerance=1050*2E-7)
  grd <- numDeriv::grad(wmeL2D$lnlf, wmeL2D$parlme)
  expect_equal(grd, rep(0,length(grd)), tolerance=tolerance)
})

context("repeating is the same as weighting: L1 replicate vs weigting, 2 REs")
test_that("L1 replicate vs weigting, 2 REs", {
  # mix for L1, weighted, 2 REs
  wmeL1WRE2 <- mix(formula=Reaction ~ Days + (1 | Subject) + (0+Days|Subject),
                   data=ss1, weights=c("W1", "W2"), nQuad=13, run=FALSE, verbose=FALSE)
  
  # mix for L1, duplicated, 2 REs
  wmeL1DRE2 <- mix(formula=Reaction ~ Days + (1 | Subject) + (0+Days|Subject),
                   data=ss2, weights=c("W1", "W2"),nQuad=13, run=FALSE, verbose=FALSE)
  
  expect_equal(wmeL1WRE2$lnlf(wmeL1DRE2$parlme), -1018.29298875158, tolerance=1050*2E-7)
  grd <- numDeriv::grad(wmeL1WRE2$lnlf, wmeL1DRE2$parlme)
  expect_equal(grd, rep(0,length(grd)), tolerance=tolerance)
  
  expect_equal(wmeL1DRE2$lnlf(wmeL1DRE2$parlme), -1018.29298875158, tolerance=1050*2E-7)
  grd <- numDeriv::grad(wmeL1DRE2$lnlf, wmeL1DRE2$parlme)
  expect_equal(grd, rep(0,length(grd)), tolerance=tolerance)
})

context("repeating is the same as weighting: L2 replicate vs weigting, 2 REs")
test_that("L2 replicate vs weigting, 2 REs", {
  # mix for L1, duplicated, 2 REs
  system.time(wmeL2DRE2 <- mix(formula=Reaction ~ Days + (1 | Subject) + (0+Days|Subject), data=ss3,
                               weights=c("W1", "W2"),nQuad=13, run=FALSE, verbose=FALSE))
  
  # mix for L1, weighted, 2 REs
  system.time(wmeL2WRE2 <- mix(formula=Reaction ~ Days + (1 | Subject) + (0+Days|Subject), data=ss4,
                               weights=c("W1", "W2"),nQuad=13, run=FALSE, verbose=FALSE))
  
  expect_equal(wmeL2DRE2$lnlf(wmeL2DRE2$parlme), -1027.89702466404, tolerance=1050*2E-7)
  grd <- numDeriv::grad(wmeL2DRE2$lnlf, wmeL2DRE2$parlme)
  expect_equal(grd, rep(0,length(grd)), tolerance=tolerance)
  
  expect_equal(wmeL2WRE2$lnlf(wmeL2DRE2$parlme), -1027.89702466404, tolerance=1050*2E-7)
  grd <- numDeriv::grad(wmeL2WRE2$lnlf, wmeL2DRE2$parlme)
  expect_equal(grd, rep(0,length(grd)), tolerance=tolerance*100) # note larger tolerance
})


ssB <- sleepstudy
set.seed(2)
ssB$Reaction <- ssB$Days * 3.141 + rnorm(nrow(ssB))
ssB$W2 <- 1
ssB$W1 <- 1:nrow(ssB)

context("Zero variance estimate in lmer")
test_that("simple model with zero variance estimate", {
  # this has 0 variance estimate in lmer
  # lmeB <- summary(lmer(Reaction ~ Days + (1|Subject), data=ssB))
  suppressWarnings(mixB <- mix(formula=Reaction ~ Days + (1 | Subject), data=ssB,
                               weights=c("W1", "W2"),  nQuad=13, run=TRUE, fast=TRUE,  verbose=FALSE))
  expect_equal(coef(mixB), structure(c(0.266730029137198, 3.09719970530335), .Names = c("(Intercept)", "Days")))
  expect_equal(mixB$vars, structure(c(0.0794570533955131, 1.18418346832652), .Names = c("Subject:(Intercept)", "Residual")))
  expect_equal(mixB$CMEAN, list(NULL, structure(c(0.277638347944171, 0.134014825632883, 0.18924017723366, -0.380043840177177, -0.0991323869545116, 0.183286617901737, -0.453643346335159, 0, -0.444084506833831, -0.223777920669184, 0.293971931396789, 0.260900976526689, -0.216206172798319, 0.405651983183519, -0.206689174876937, 0.209466503252012, -0.326652547838431, -0.199452174848381), .Dim = c(18L, 1L))))
})
