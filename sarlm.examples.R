library(spdep)
source("sarlm.R")

data(oldcol)
lw <- nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw, type="mixed")

# I - In-sample predictors

# defaut predictor: trend-signal decomposition
p1 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F)
p1_bis = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "TS") # same call
# X: only the trend
p2 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "trend")
# TC
p3 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "TC")
# BP
p4 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "BP")


# TODO: Prediction


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

newdata = COL.OLD[c("1042", "1043", "1044", "1045"), c("INC", "HOVAL")]
newdata = COL.OLD[c("1042", "1043", "1044", "1045"), ] # work with newdata having other columns in a different order than formula

# out-sample predictors
COL.OLD[region.id.newdata, "CRIME"] # true values

po1 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T)
po2 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TS1")
po3 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC")
po4 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC", all.data = T) # returns for all spatial units
#TODO: need to be check. We do not recover the same predictions. This is link to the fact that lw is used directly in the wrong order
po5 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "BP")
po6 = predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "trend")


system.time(po3 <- predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC", all.data = F))
system.time(po4 <- predict.sarlm(COL.lag.eig, listw = lw, newdata = newdata, zero.policy = T, type = "TC", all.data = T))
