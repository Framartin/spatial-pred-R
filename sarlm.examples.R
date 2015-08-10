library(spdep)
source("sarlm.R")
source("nb2mat.R")


data(oldcol)
lw <- nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw, type="mixed")

# I - In-sample predictors

# defaut predictor: trend-signal decomposition
p1A = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F)
p1B = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "TS") # same call
# X: only the trend
p2 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "trend")
# TC
p3 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "TC")
# BP
p4 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "BP")


# forecast
p1 = predict.sarlm(COL.lag.eig, listw = lw, zero.policy = F, newdata = COL.OLD+1, type = "trend")
p2 = predict.sarlm(COL.lag.eig, listw = lw, zero.policy = F, newdata = COL.OLD+1, type = "TC")

p1 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, newdata = COL.OLD+1, type = "trend")
p2 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, newdata = COL.OLD+1, type = "TC")



# III - Out-of-sample predictions


# decomposing data in in-sample and out-of-sample
listw = lw
region.id = attr(listw, "region.id")
region.id.data = region.id[!region.id %in% c("1042", "1043", "1044", "1045")]
region.id.newdata = region.id[region.id %in% c("1042", "1043", "1044", "1045")]
listw.d = .listw.decompose(listw, region.id.data, region.id.newdata, type = c("Wss", "Wos", "Wso", "Woo")) # submatrices in the order of region.id.data and region.id.newdata

# fit model
data = COL.OLD[!rownames(COL.OLD) %in% c("1042", "1043", "1044", "1045"), ]
listw.sub <- subset(listw, !region.id %in% c("1042", "1043", "1044", "1045"))
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=data, listw.sub)
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=data, listw.sub, type="mixed")
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=data, listw.sub)
COL.SDEMW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=data, listw.sub, etype = "emixed")



newdata = COL.OLD[c("1042", "1043", "1044", "1045"), c("INC", "HOVAL")]
newdata = COL.OLD[c("1042", "1043", "1044", "1045"), ] # work with newdata having other columns in a different order than formula

# out-sample predictors
COL.OLD[region.id.newdata, "CRIME"] # true values

po1 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T)
po2A = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TS1")
po2B = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "KP4")
po3 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC")
po4 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC", all.data = T) # returns for all spatial units
po5 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "BP")
po6 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "trend")
po7A = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC1")
po7B = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "KP1")
po7C = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC1", power = FALSE)

all(po3 - po4[46:49] < 0.0001)

system.time(po3 <- predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC", all.data = F))
system.time(po4 <- predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC", all.data = T))

# aliased coefs
data2 <- data
data2[,1] <- 1
COL.mix.eig <- lagsarlm(CRIME ~ AREA_PL + INC + HOVAL, data=data2, listw.sub, type="mixed")
COL.mix.eig$aliased
po1 = predict.sarlm(COL.mix.eig, listw = lw, newdata = newdata, zero.policy = T)
po2 = predict.sarlm(COL.mix.eig, listw = lw, newdata = newdata, zero.policy = T, legacy.mixed = T)
po3 = predict.sarlm(COL.mix.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC")
po4 = predict.sarlm(COL.mix.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC", legacy.mixed = T)


