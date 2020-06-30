# Load appropriate packages
library(rethinking)
library(R.matlab)

# LOCAL SENSITIVITY ANALYSIS

# Read in 
SADlo.orig <- readMat('SensitivityAnalysis_LoDens_Wig1.mat')
SADlo.pars <- data.frame(SADlo.orig$SAD[,,1]) # These are the parameter combos.
SADlo.eq <- data.frame(SADlo.orig$SAD[,,2])   # These are the resulting equilibria
SADhi.orig <- readMat('SensitivityAnalysis_HiDens_Wig1.mat')
SADhi.pars <- data.frame(SADhi.orig$SAD[,,1]) # These are the parameter combos.
SADhi.eq <- data.frame(SADhi.orig$SAD[,,2])   # These are the resulting equilibria

# Check that all parameter combinations led to the correct number of equilibria
which(SADlo.eq$X21 != 1)
which(SADhi.eq$X21 != 2)

# Reorder the equilibria in SADhi.eq so that the equilibrium in which Th1 wins always comes first
for (i in 1:nrow(SADhi.eq)) {
  if (SADhi.eq[i,1] < 1) {
    SADhi.eq[i, ] <- SADhi.eq[i, c(5,6,7,8,1,2,3,4,9:21)]
  }
}

# Delete all but the first 4 columns of SADlo.eq (since there is exactly 1 eq per parameter combo)
# and all but the first 8 columns of SADhi.eq (since there are exactly 2 eq per parameter combo)
SADlo.eq <- SADlo.eq[ , 1:4]
names(SADlo.eq) <- c("coor1", "coor2", "coor3", "coor4")
SADhi.eq <- SADhi.eq[ , 1:8]
names(SADhi.eq) <- c("coor1a", "coor2a", "coor3a", "coor4a", "coor1b", "coor2b", "coor3b", "coor4b")

# Delete all the even numbered columns of SADlo.pars (since parameter symmetry is enforced)
SADlo.pars <- SADlo.pars[ , seq(1,21,2)]
names(SADlo.pars) <- c("B","H","F","E","W","D","G","V","Q","L","Asym")
SADhi.pars <- SADhi.pars[ , seq(1,21,2)]
names(SADhi.pars) <- c("B","H","F","E","W","D","G","V","Q","L","Asym")

# Finally, merge the two data frames into one.
SADlo <- cbind(SADlo.pars, SADlo.eq) 
SADhi <- cbind(SADhi.pars, SADhi.eq) 

# Standardize all predictor variables.
SADlo.s <- SADlo
SADhi.s <- SADhi
for (i in seq(from=1, to=11, by=1)) {
  SADlo.s[ , i] <- (SADlo[ , i] - mean(SADlo[ , i]))/sd(SADlo[ , i])
  SADhi.s[ , i] <- (SADhi[ , i] - mean(SADhi[ , i]))/sd(SADhi[ , i])
}

# Linear regression of coordinate 1 of the equilibrium on the predictor variables.
SADlo.s.mod <- map(alist(coor1 ~ dnorm(mu, sigma),
                         mu <- a + bb*B + bh*H + bf*F + be*E + bw*W + bd*D + bg*G + bv*V + bq*Q + bl*L +ba*Asym,
                         a ~ dnorm(2.3, 1),
                         bb ~ dnorm(0,1),
                         bh ~ dnorm(0,1),
                         bf ~ dnorm(0,1),
                         be ~ dnorm(0,1),
                         bw ~ dnorm(0,1),
                         bd ~ dnorm(0,1),
                         bg ~ dnorm(0,1),
                         bv ~ dnorm(0,1),
                         bq ~ dnorm(0,1),
                         bl ~ dnorm(0,1),
                         ba ~ dnorm(0,1),
                         sigma ~ dunif(0,1)),
                   data=SADlo.s)
precis(SADlo.s.mod, digits=3)

SADhi.s.mod <- map(alist(coor1a ~ dnorm(mu, sigma),
                         mu <- a + bb*B + bh*H + bf*F + be*E + bw*W + bd*D + bg*G + bv*V + bq*Q + bl*L +ba*Asym,
                         a ~ dnorm(2.3, 1),
                         bb ~ dnorm(0,1),
                         bh ~ dnorm(0,1),
                         bf ~ dnorm(0,1),
                         be ~ dnorm(0,1),
                         bw ~ dnorm(0,1),
                         bd ~ dnorm(0,1),
                         bg ~ dnorm(0,1),
                         bv ~ dnorm(0,1),
                         bq ~ dnorm(0,1),
                         bl ~ dnorm(0,1),
                         ba ~ dnorm(0,1),
                         sigma ~ dunif(0,1)),
                   data=SADhi.s)
precis(SADhi.s.mod, digits=3)

# Prepare data for bar plot giving effect sizes for sensitivity analysis
SADlo.eff <- coef(SADlo.s.mod)
SADlo.eff <- SADlo.eff[c(2:12)] #Remove the intercept coefficient.
SADlo.eff <- abs(SADlo.eff)     #Take the absolute value of effect size.
names(SADlo.eff) <- c("B","H","F","E","W","D","G","V","Q","L","Asym")
SADlo.eff <- sort(SADlo.eff, decreasing = TRUE)
SADlo.eff <- 100*SADlo.eff/mean(SADlo.s$coor1) #Put all effect sizes in terms of their percent effect on the equilibrium.

SADhi.eff <- coef(SADhi.s.mod)
SADhi.eff <- SADhi.eff[c(2:12)] #Remove the intercept coefficient.
SADhi.eff <- abs(SADhi.eff)     #Take the absolute value of effect size.
names(SADhi.eff) <- c("B","H","F","E","W","D","G","V","Q","L","Asym")
SADhi.eff <- sort(SADhi.eff, decreasing = TRUE)
SADhi.eff <- 100*SADhi.eff/mean(SADhi.s$coor1a) #Put all effect sizes in terms of their percent effect on the equilibrium.

# Bar plot giving effect sizes for sensitivity analysis
barplot(SADlo.eff, ylim=c(0,6), ylab="Effect Size (% change in equilibrium position)", xlab="Non-Dimensional Parameter", 
        main="Sensitivity of Equilibrium Position to Parameters at Low Cell Density")
barplot(SADhi.eff, ylim=c(0,6), ylab="Effect Size (% change in equilibrium position)", xlab="Non-Dimensional Parameter", 
        main="Sensitivity of Equilibria Positions to Parameters at High Cell Density")

# PSEUDO-GLOBAL EQUILIBRIUM ANALYSIS

# Read in all the necessary data
Lo.1 <- readMat('SensitivityAnalysis_LoDens_Wig1.mat')
Lo.2 <- readMat('SensitivityAnalysis_LoDens_Wig2.mat')
Lo.3 <- readMat('SensitivityAnalysis_LoDens_Wig3.mat')
Lo.4 <- readMat('SensitivityAnalysis_LoDens_Wig4.mat')
Lo.5 <- readMat('SensitivityAnalysis_LoDens_Wig5.mat')
Lo.6 <- readMat('SensitivityAnalysis_LoDens_Wig6.mat')
Lo.7 <- readMat('SensitivityAnalysis_LoDens_Wig7.mat')
Lo.8 <- readMat('SensitivityAnalysis_LoDens_Wig8.mat')
Lo.9 <- readMat('SensitivityAnalysis_LoDens_Wig9.mat')
Hi.1 <- readMat('SensitivityAnalysis_HiDens_Wig1.mat')
Hi.2 <- readMat('SensitivityAnalysis_HiDens_Wig2.mat')
Hi.3 <- readMat('SensitivityAnalysis_HiDens_Wig3.mat')
Hi.4 <- readMat('SensitivityAnalysis_HiDens_Wig4.mat')
Hi.5 <- readMat('SensitivityAnalysis_HiDens_Wig5.mat')
Hi.6 <- readMat('SensitivityAnalysis_HiDens_Wig6.mat')
Hi.7 <- readMat('SensitivityAnalysis_HiDens_Wig7.mat')
Hi.8 <- readMat('SensitivityAnalysis_HiDens_Wig8.mat')
Hi.9 <- readMat('SensitivityAnalysis_HiDens_Wig9.mat')

# Select out the useful data for this analysis.
Lo.eq <- data.frame(Lo.1$SAD[,21,2], Lo.2$SAD[,21,2], Lo.3$SAD[,21,2], Lo.4$SAD[,21,2], Lo.5$SAD[,21,2],
                    Lo.6$SAD[,21,2], Lo.7$SAD[,21,2], Lo.8$SAD[,21,2], Lo.9$SAD[,21,2])
names(Lo.eq) <- c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%")
Hi.eq <- data.frame(Hi.1$SAD[,21,2], Hi.2$SAD[,21,2], Hi.3$SAD[,21,2], Hi.4$SAD[,21,2], Hi.5$SAD[,21,2],
                    Hi.6$SAD[,21,2], Hi.7$SAD[,21,2], Hi.8$SAD[,21,2], Hi.9$SAD[,21,2])
names(Hi.eq) <- c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%")

# Reframe the data into a stacked-barplot-friendly format.
Lo.eq.plot <- data.frame(matrix(ncol=2, nrow=9*1000))
Hi.eq.plot <- data.frame(matrix(ncol=2, nrow=9*1000))
names(Lo.eq.plot) <- c("NumEq", "ParPercVar")
names(Hi.eq.plot) <- c("NumEq", "ParPercVar")
Lo.eq.plot$NumEq <- c(Lo.eq[,1], Lo.eq[,2], Lo.eq[,3], Lo.eq[,4], Lo.eq[,5], 
                       Lo.eq[,6], Lo.eq[,7], Lo.eq[,8], Lo.eq[,9])
Lo.eq.plot$ParPercVar <- c(rep.int(10,1000), rep.int(20,1000), rep.int(30,1000),rep.int(40,1000), 
                           rep.int(50,1000), rep.int(60,1000), rep.int(70,1000), rep.int(80,1000),
                           rep.int(90,1000))
Hi.eq.plot$NumEq <- c(Hi.eq[,1], Hi.eq[,2], Hi.eq[,3], Hi.eq[,4], Hi.eq[,5], 
                       Hi.eq[,6], Hi.eq[,7], Hi.eq[,8], Hi.eq[,9])
Hi.eq.plot$ParPercVar <- c(rep.int(10,1000), rep.int(20,1000), rep.int(30,1000), rep.int(40,1000), 
                           rep.int(50,1000), rep.int(60,1000), rep.int(70,1000), rep.int(80,1000),
                           rep.int(90,1000))
Lo.Table <- table(Lo.eq.plot$NumEq, Lo.eq.plot$ParPercVar)
Hi.Table <- table(Hi.eq.plot$NumEq, Hi.eq.plot$ParPercVar)

# Make stacked bar plots
barplot(Lo.Table, yaxt="n", ylab="Proportion of Samples", xlab="Percent Variation in Parameter Values", 
        main = "Randomly Sampled Parameter Combinations at Low Cell Density", col=gray.colors(5))
axis(side=2, at=c(0,200,400,600,800,1000), labels=c("0%","20%","40%","60%","80%","100%"))
legend("bottomright", inset=0.05, legend=c("1","2", "3", "4", "5"), title="Number of Stable Equilibria",
       horiz=TRUE, fill=gray.colors(5))
barplot(Hi.Table, yaxt="n", ylab="Proportion of Samples", xlab="Percent Variation in Parameter Values", 
        main = "Randomly Sampled Parameter Combinations at High Cell Density", col=gray.colors(5))
axis(side=2, at=c(0,200,400,600,800,1000), labels=c("0%","20%","40%","60%","80%","100%"))
legend("bottomright", inset=0.05, legend=c("1","2", "3", "4", "5"), title="Number of Stable Equilibria",
       horiz=TRUE, fill=gray.colors(5))











