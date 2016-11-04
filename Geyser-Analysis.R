library(MASS)
library(rgl)
library(sm)
detach(faithful)


############### Data Preparation #####################

# Load data set and show some preliminary summary
g <- geyser
str(g)
summary(g)

l <- length(g$duration)

# waiting time is measured starting at the starting point of the former
# eruption. We get rid of this mixing. We cannot do it for the first
# waiting time, so we leave it as it is.
g$nextwaiting <- c(g$waiting[1], g$waiting[-1] - g$duration[-l])
g$durgrouping <- factor(ifelse(duration <= 3.1, "short", "long"))
g$waitgrouping <- factor(ifelse(waiting <= 67, "short", "long"))

try(detach(g))
attach(g)

# predictive function used to predict the next waiting time.
# takes the duration of the last eruption.
fpred <- function(d) {
  
   30 + 10 *d
}

dwait <- function(x, p, u1, s1, u2, s2) {
  
  e <- p * dnorm(x, mean = u1, sd = s1) + (1-p)* dnorm(x, mean = u2, sd = s2)
}

dwaitoptim <- function(p,x) {
  
  e <- p[1]*dnorm(x, p[2], p[3]) + (1-p[1])*dnorm(x, p[4], p[5])
  if(any( e <= 0)) Inf else -sum(log(e))
}

############### Preliminary Plots ######################

#Plot waiting time
plot(waiting)

# Plot waiting time over duration of geyser
plot(duration[-l], waiting[-1])
abline(v=3.1)

plot(waiting[-l], waiting[-1])
abline(v=67)
plot3d(g$waiting[-l], g$duration[-l], g$waiting[-1])

# Same plots with nextwaiting
plot(nextwaiting)
plot(duration[-l], nextwaiting[-1])
plot(nextwaiting[-l], nextwaiting[-1])
plot3d(nextwaiting[-l], duration[-l], nextwaiting[-1])

# Plot time series itself
# Seems that long and short duratoin seem to alternate
plot(duration[1:100], type="l")
plot(nextwaiting[1:100], type="l")

resid.fpred <- waiting[-1] - fpred(duration[-l])
plot(resid.fpred)

# Density plot of stuff
sm.density(cbind(waiting = g$waiting[-l], 
                 nextwaiting = g$waiting[-1],duration = g$duration[-l]))

############ Linear Models ########################

# 1) unreduced waiting time first
# 1.1) duration
# 1.1a) no grouping
lm.duration <- lm(waiting[-1] ~ duration[-l], data =g)
summary(lm.duration) # -14.6940
plot(duration[-l], waiting[-1])
abline(lm.duration)
plot(resid(lm.duration))
oldpar <- par(mfrow=c(3,2))
plot(lm.duration,which=1:6)
par(oldpar)


# 1.1b) with grouping
lm.durgrouping <- lm(waiting[-1] ~ durgrouping[-l] * duration[-l], data =g)
summary(lm.durgrouping) # -14.8288
plot(resid(lm.durgrouping))
oldpar <- par(mfrow=c(3,2))
plot(lm.durgrouping,which=1:6)
par(oldpar)
sort(resid(lm.durgrouping)) # 186, 210, 258 to -11.153

coloring <- rep("black", l)
coloring[186] <- "red"
coloring[210] <- "green"
coloring[258] <- "pink"
plot(duration[-l], waiting[-1], col = coloring)
plot(waiting[-l], waiting[-1], col = coloring)
plot3d(waiting[-l], duration[-l], waiting[-1], col = coloring)

# 1.2) duration + waitingtime
# 1.2a) nogrouping
lm.durwait <- lm(waiting[-1] ~ duration[-l] + waiting[-l], data = g)
summary(lm.durwait) # -16.5225
oldpar <- par(mfrow=c(3,2))
plot(lm.durwait,which=1:6)
par(oldpar)
plot(resid(lm.durwait))
which.min(resid(lm.durwait)) # 62
sort(resid(lm.durwait)) # 62, 133 (-12.65981947)

coloring <- rep("black", l)
coloring[62] <- "red" 
coloring[133] <- "green"
plot(duration[-l], waiting[-1], col = coloring)
plot(waiting[-l], waiting[-1], col = coloring)
plot3d(waiting[-l], duration[-l], waiting[-1], col = coloring)
g[c(62, 133),]

# 1.2b) duration grouping
lm.durwaitgrouping1 <- lm(waiting[-1] ~ durgrouping[-l] * 
                            (duration[-l] + waiting[-l]), data = g)
summary(lm.durwaitgrouping1) # -12.61968
oldpar <- par(mfrow=c(3,2))
plot(lm.durwaitgrouping1,which=1:6)
par(oldpar)
plot(resid(lm.durwaitgrouping1))
sort(resid(lm.durwaitgrouping1)) # 133, 186, 62 (without down to -8.45)

coloring <- rep("black", l)
coloring[62] <- "red" 
coloring[133] <- "green"
coloring[186] <- "orange"
plot(duration[-l], waiting[-1], col = coloring)
plot(waiting[-l], waiting[-1], col = coloring)
plot3d(waiting[-l], duration[-l], waiting[-1], col = coloring)
g[c(62, 133, 186),]

# 1.2b) duration and waiting grouping
lm.durwaitgrouping2 <- lm(waiting[-1] ~ durgrouping[-l] * waitgrouping[-l]
                          * (duration[-l] + waiting[-l]), data = g)
summary(lm.durwaitgrouping2) # 12.6197
oldpar <- par(mfrow=c(3,2))
plot(lm.durwaitgrouping2,which=1:6)
par(oldpar)
plot(resid(lm.durwaitgrouping2))
sort(resid(lm.durwaitgrouping2)) # 133, 186, 62 down to -8.79
# worse than only duration grouping let that be.

# 1.2c) waiting grouping
lm.durwaitgrouping3 <- lm(waiting[-1] ~ waitgrouping[-l]
                          * (duration[-l] + waiting[-l]), data = g)
summary(lm.durwaitgrouping3) # - 13.1815
oldpar <- par(mfrow=c(3,2))
plot(lm.durwaitgrouping3,which=1:6)
par(oldpar)
plot(resid(lm.durwaitgrouping3))
sort(resid(lm.durwaitgrouping3)) # 133, 111, 186 only to -10.49

# 1.2d) additiv
lm.durwaitgrouping4 <- lm(waiting[-1] ~ waitgrouping[-l]*duration[-l] 
                          + waiting[-l], data =g)
summary(lm.durwaitgrouping4) # -12.9146
oldpar <- par(mfrow=c(3,2))
plot(lm.durwaitgrouping4,which=1:6)
par(oldpar)
plot(resid(lm.durwaitgrouping4))

# 1.2e) two groupings, but differently distributed
lm.durwaitgrouping5 <- lm(waiting[-1] ~ durgrouping[-l] * (duration[-l] 
                          + waitgrouping[-l]* waiting[-l]), data = g)
summary(lm.durwaitgrouping5) # 12.6197
plot(resid(lm.durwaitgrouping5))
sort(resid(lm.durwaitgrouping5)) # 133, 186, 62 then -8.68

# 1.3) second level
# 1.3a) 4 variables no grouping
first <- c(-1,-2)
last <- c(-(l-1), -l)
lm.2full <- lm(waiting[first] ~ duration[last] + waiting[last]
               + duration[c(-1,-l)] + waiting[c(-1,-l)], data = g)
summary(lm.2full) # -15.7157
plot(resid(lm.2full))
sort(resid(lm.2full)) # 61, 132, 202 get to -10.88

# 1.3b) no last waiting
lm.2nowaiting <- lm(waiting[first] ~ duration[last] 
               + duration[c(-1,-l)] + waiting[c(-1,-l)], data = g)
summary(lm.2nowaiting) # -16.814
plot(resid(lm.2nowaiting)) # not good enough

# 1.3c) add duration grouping
lm.2fulldurgrouping <- lm(waiting[first] ~ durgrouping[c(-1,-l)] * 
                            (duration[last] + waiting[last] +
                               duration[c(-1,-l)] + waiting[c(-1,-l)]))
summary(lm.2fulldurgrouping) # -12.401
plot(resid(lm.2fulldurgrouping))
sort(resid(lm.2fulldurgrouping))# 132, 185, 257, 61 to -8.17

coloring <- rep("black", l)
coloring[132] <- "red"
coloring[185] <- "green"
coloring[257] <- "orange"
coloring[61] <- "pink"
plot3d(waiting[-l], duration[-l], waiting[-1], col = coloring)
plot(waiting[-l], waiting[-1], col = coloring)
plot(duration[-l], waiting[-1], col = coloring)

# 1.3d) two duration groupings
lm.2fulldurgrouping2 <- lm(waiting[first] ~ durgrouping[last] * 
                             (duration[last] + waiting[last]) +
                             durgrouping[c(-1,-l)] *
                                (duration[c(-1,-l)] + waiting[c(-1,-l)]))
summary(lm.2fulldurgrouping2) # -12.3733
plot(resid(lm.2fulldurgrouping2))
sort(resid(lm.2fulldurgrouping2)) # 185, 132, 66, 257 to -7.98
oldpar <- par(mfrow=c(3,2))
plot(lm.2fulldurgrouping2,which=1:6)
par(oldpar)

# 2) reduced waiting time

# Look at entry 62
g[62,]
coloring <- rep("black", l)
coloring[62] <- "red"
plot(duration[-l], waiting[-1], col = coloring)
plot(waiting[-l], waiting[-1], col =coloring)
plot3d(waiting[-l], duration[-l], waiting[-1], col = coloring)
# seems to have a very long duration. Might be an outlier.

# only waiting grouping

############### Fit Distribution ############################

truehist(waiting, h=5, ymax = 0.04)
wait.dns <- density(waiting, n = 512, width = "SJ")
lines(wait.dns)

p0 <- c(p = mean(waiting < 70), u1 = 50, s1 = 5, u2 = 80, s2 = 5)

p1 <- optim(p0, dwaitoptim, x = waiting)$par
detach(as.list(p1))
attach(as.list(p1))

wait.fdns <- list(x = wait.dns$x, 
                  y = dwait(wait.dns$x, p, u1, s1, u2, s2))
lines(wait.fdns)

