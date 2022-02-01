library(mvabund)
data(reveg)
reveg$abundMV=mvabund(reveg$abund)
abundMV=mvabund(reveg$abund)
ft_reveg=manyglm(abundMV~treatment+offset(log(pitfalls)),
                 family="negative.binomial", data=reveg) # offset included as in Exercise 10.9

tr=reveg$treatment
logPt=log(reveg$pitfalls)
ft_reveg=manyglm(abundMV~tr,offset=logPt, family="negative.binomial") # offset included as in Exercise 10.9
plotenvelope(ft_reveg,which=1,n.sim=99)

abund = mvabund(reveg$abund)
X = data.frame(treatment=reveg$treatment, pitfalls=reveg$pitfalls)
ft_reveg=manyglm(abund~treatment,offset=log(pitfalls),
                 family="negative.binomial", data=X)
#plotenvelope(ft_reveg,which=1)
library(ecoCopula)
cord_reveg=cord(ft_reveg)
simulate(cord_reveg)


#library(ecostats)
#plotenvelope(glm.spid, which=1, n.sim=59)
#simulate(ft_cord)

example(manyglm)
#c=cord(glm.spid)
#plotenvelope(glm.spid,which=2,n.sim=99)

ft_reveg=manyglm(abundMV~treatment+offset(log(pitfalls)),
                 family="negative.binomial", data=reveg) # offset included as in Exercise 10.9
#plotenvelope(ft_reveg, which=1, n.sim=59)

                 
data(iris)
irisY = as.matrix(iris[,1:4])
ft=lm(log(irisY[,1:2])~Species,data=iris)
plotenvelope(ft,which=1:2)
# should return no error

library(psych)
fa_iris <- fa(iris[,1:4], nfactors=2,fm="ml",rotate="varimax")

iVar=1
plotenvelope(lm(iris[,iVar]~fa_iris$scores), which=1:3, col=iris$Species)
# should return no error
