data(oldcol)
lw <- nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw, type="mixed")
p1 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F)
p2 = predict.sarlm(COL.mix.eig, listw = lw, zero.policy = F, type = "BP")

# loo
newdata = matrix(c(1,2), nrow=1)
rownames(newdata) = "1003"
colnames(newdata) = c('INC', 'HOVAL')


# decomposing data in in-sample and out-of-sample
region.id = attr(listw, "region.id")
region.id.data = region.id[!region.id %in% c("1003", "1005", "1042")]
region.id.newdata = region.id[region.id %in% c("1003", "1005", "1042")]
listw.d = .listw.decompose(listw, region.id.data, region.id.newdata, type = c("Wss", "Wos", "Wso", "Woo"))
