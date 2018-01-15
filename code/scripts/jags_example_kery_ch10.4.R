
setwd("/Users/ali/Documents/REPOS/seasonal_occupancy/code/scripts/jags_test/")


# Choose sample sizes and prepare obs. data array y
set.seed(1)                   # So we all get same data set
M <- 100                      # Number of sites
J <- 3                        # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for occupancy model and compute occupancy
beta0 <- 0                    # Logit-scale intercept
beta1 <- 3                    # Logit-scale slope for vegHt
psi <- plogis(beta0 + beta1 * vegHt) # Occupancy probability
# plot(vegHt, psi, ylim = c(0,1), type = "l", lwd = 3) # Plot psi relationship

# Now visit each site and observe presence/absence perfectly
z <- rbinom(M, 1, psi)        # True presence/absence

# Look at data so far
table(z)

# Plot the true system state
par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
plot(vegHt, z, xlab="Vegetation height", ylab="True presence/absence (z)", frame = F, cex = 1.5)
plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red")

# Create a covariate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2                        # Logit-scale intercept
alpha1 <- -3                        # Logit-scale slope for wind
p <- plogis(alpha0 + alpha1 * wind) # Detection probability
# plot(p ~ wind, ylim = c(0,1))     # Look at relationship

# Take J = 3 presence/absence measurements at each site
for(j in 1:J) {
    y[,j] <- rbinom(M, z, p[,j])
}
sum(apply(y, 1, max))               # Number of sites with observed presences

# Plot observed data and true effect of wind on detection probability
plot(wind, y, xlab="Wind", ylab="Observed det./nondetection data (y)", frame = F, cex = 1.5)
plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")

# Look at the data: occupancy, true presence/absence (z), and measurements (y)
cbind(psi=round(psi,2), z=z, y1=y[,1], y2=y[,2], y3=y[,3])

# Create factors
time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
hab <- c(rep("A", 33), rep("B", 33), rep("C", 34))  # Must have M = 100


# Load unmarked, format data and summarize
library(unmarked)
umf <- unmarkedFrameOccu(
   y = y,                                            # Pres/Abs measurements
   siteCovs = data.frame(vegHt = vegHt, hab = hab),  # site-specific covs.
   obsCovs = list(wind = wind, time = time))         # obs-specific covs.
summary(umf)


# Fit model and extract estimates
# Detection covariates follow first tilde, then occupancy covariates
summary(fm.occ1 <- occu(~wind ~vegHt, data=umf))

# Fit model p(time+wind), psi(hab+vegHt)
summary(fm2.occ <- occu(~time+wind-1 ~hab+vegHt-1, data=umf))

# Bundle and summarize data set
str( win.data <- list(y = y, vegHt = vegHt, wind = wind, M = nrow(y), J = ncol(y), XvegHt = seq(-1, 1, length.out=100), Xwind = seq(-1, 1, length.out=100)) )

# Specify model in BUGS language
sink("model.txt")
cat("
model {

# Priors
mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
alpha0 <- logit(mean.p)      # Detection intercept
alpha1 ~ dunif(-20, 20)      # Detection slope on wind
mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
beta0 <- logit(mean.psi)     # Occupancy intercept
beta1 ~ dunif(-20, 20)       # Occupancy slope on vegHt

# Likelihood
for (i in 1:M) {
   # True state model for the partially observed true state
   z[i] ~ dbern(psi[i])      # True occupancy z at site i
   logit(psi[i]) <- beta0 + beta1 * vegHt[i]
   for (j in 1:J) {
      # Observation model for the actual observations
      y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j
      p.eff[i,j] <- z[i] * p[i,j]   # 'straw man' for WinBUGS
      logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
   }
}

# Derived quantities
N.occ <- sum(z[])       # Number of occupied sites among sample of M
psi.fs <- N.occ/M       # Proportion of occupied sites among sample of M
for(k in 1:100){
   logit(psi.pred[k]) <- beta0 + beta1 * XvegHt[k] # psi predictions
   logit(p.pred[k]) <- alpha0 + alpha1 * Xwind[k]  # p predictions
}
}
",fill = TRUE)
sink()

# Initial values: must give for same quantities as priors given !
zst <- apply(y, 1, max)        # Avoid data/model/inits conflict
inits <- function(){list(z = zst, mean.p = runif(1), alpha1 = runif(1), mean.psi = runif(1), beta1 = runif(1))}

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "N.occ", "psi.fs", "psi.pred", "p.pred", "z") # Also estimate z = "conditional occ. prob."

# MCMC settings
ni <- 25000   ;   nt <- 10   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART 2 min) and summarize posteriors
# out1B <- bugs(win.data, inits, params, "model.txt", n.chains = nc,
# n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
# print(out1B, dig = 3)

out1B <- jags.model("model.txt", data=win.data, inits=inits,
                n.chains = nc, n.adapt=nb, quiet=FALSE)

# Compare truth with MLEs and bayesian posterior inference in table ...
truth <- c(alpha0, alpha1, beta0, beta1, sum(z), sum(z)/M)
tmp <- summary(fm.occ1)
MLEs <- rbind(tmp[[2]][1:2,1:2], tmp[[1]][1:2,1:2], sumZ = c(mean(pb.fs@t.star[,1]), sd(pb.fs@t.star[,1])), psi.fs = c(mean(pb.fs@t.star[,2]), sd(pb.fs@t.star[,2])))
print(cbind(truth, MLEs, out1B$summary[1:6, 1:2]))

# .... and in a graph (Fig. 10-4)
par(mfrow = c(2, 2), mar = c(4, 5, 2, 2), las = 1, cex.lab = 1, cex = 1.2)
plot(vegHt, z, xlab="Vegetation height", ylab="Occupancy prob. (MLE)", ylim = c(0, 1), frame = F)   # True presence/absence
lines(seq(-1,1, 0.01), pred.occ[,1], col = "blue", lwd = 2)
matlines(seq(-1,1, 0.01), pred.occ[,3:4], col = "grey", lty = 1)
lines(vegHt, psi, lwd=3, col="red")   # True psi
plot(wind, y, xlab="Wind", ylab="Detection prob. (MLE)", ylim = c(0,1), frame=F)
lines(seq(-1, 1, 0.1), pred.det[,1], col = "blue", lwd = 2)
matlines(seq(-1, 1, 0.1), pred.det[,3:4], col = "grey", lty = 1)
plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")
plot(vegHt, z, xlab="Vegetation height", ylab="Occupancy prob. (P. mean)", las = 1, frame = F)   # True presence/absence
lines(vegHt, psi, lwd=3, col="red")   # True psi
lines(win.data$XvegHt, out1B$summary[7:106,1], col="blue", lwd = 2)
matlines(win.data$XvegHt, out1B$summary[7:106,c(3,7)], col="grey", lty = 1)
plot(wind, y, xlab="Wind", ylab="Detection prob. (P. mean)", frame = F)
plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")
lines(win.data$Xwind, out1B$summary[107:206,1], col="blue", lwd = 2)
matlines(win.data$Xwind, out1B$summary[107:206,c(3,7)], col="grey", lty = 1)


# Plot posterior distribution of number of occupied sites (see Fig. 10-3, right)
hist(out1B$sims.list$N.occ, col = "grey", breaks = 60, xlim = c(20, 100),
main = "", freq = F)
abline(v = out1B$mean$N.occ, col = "blue", lwd = 3)   # add point estimate
abline(v = out1B$summary[5,c(3,7)], col = "grey", lwd = 3)   # add 95% CRI
abline(v = sum(apply(y, 1, max)), lty = 2, lwd = 3)   # observed #occ sites
abline(v = sum(z), col = "red", lwd = 2)              # true #occ sites

