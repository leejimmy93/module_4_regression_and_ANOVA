### Needed packages: multcomp, gee 

## read data set ----------------------------------------------------------------
cholesterol = read.table("http://faculty.washington.edu/rhubb/sisg/SISG-Data-cholesterol.txt", header=T)
attach(cholesterol)

## Exploratory data analysis ----------------------------------------------------
## graphical display: boxplot 
boxplot(chol ~ as.factor(rs174548))
boxplot(chol ~ as.factor(rs4775401))
boxplot(chol ~ as.factor(APOE))

## alternative graphical display: graph of means 
par(mfrow=c(1,3))
plot.design(chol ~ as.factor(rs174548))
plot.design(chol ~ as.factor(rs4775401))
plot.design(chol ~ as.factor(APOE))

## numeric descriptives 
tapply(chol, as.factor(rs174548), mean)
tapply(chol, as.factor(rs174548), sd)

tapply(chol, as.factor(rs4775401), mean)
tapply(chol, as.factor(rs4775401), sd)

tapply(chol, as.factor(APOE), mean)
tapply(chol, as.factor(APOE), sd)

## an alternative way to obtain descriptives by group
by(chol, as.factor(rs174548), mean)
by(chol, as.factor(rs174548), sd)

by(chol, as.factor(rs4775401), mean)
by(chol, as.factor(rs4775401), sd)

by(chol, as.factor(APOE), mean)
by(chol, as.factor(APOE), sd)

## a little fancier way to obtain mean and standard deviation at once
tapply(chol, as.factor(rs174548), function(x){c(mean=mean(x), sd=sd(x))})
by(chol, as.factor(rs174548), function(x){c(mean=mean(x), sd=sd(x))})


## Inferential data analysis ----------------------------------------------------

## One-way ANOVA by hand (creating your own dummy variables)
dummy1 = 1*(rs174548==1)
dummy2 = 1*(rs174548==2)
fitd = lm(chol ~ dummy1 + dummy2)
summary(fitd)
anova(fitd)

## telling R that a variable is categorical
fit1 = lm(chol ~ -1+as.factor(rs174548))
summary(fit1)
fit0 = lm(chol ~1)
anova(fit0,fit1)

## here the predictor is treated as "continuous"
fit2 = lm(chol ~ rs174548)
summary(fit2)
anova(fit2)

## now, let's consider multiple comparisons for the oneway anova model -- need package!
library(multcomp)

## all pairwise comparisons
K = contrMat(table(rs174548), type="Tukey")

fit1 = lm(chol ~ -1+as.factor(rs174548))

mc = glht(fit1, linfct =K)
summary(mc, test=adjusted("none"))
summary(mc, test=adjusted("bonferroni"))


## alternatively, define contrast matrix directly
contr = rbind("mean(C/G+G/G) - mean(C/C)" = c(-2, 1, 1))
mc2 = glht(fit1, linfct =contr)
summary(mc2, test=adjusted("none"))

## more than one contrast (again user-defined)
contr2 = rbind("mean(C/G+G/G) - mean(C/C)" = c(-2, 1, 1),
               "mean(C/C+G/G) - mean(C/G)" = c(1, -2, 1))
mc3 = glht(fit1, linfct =contr2)
summary(mc3, test=adjusted("none"))
summary(mc3, test=adjusted("bonferroni"))

## testing variances
bartlett.test(chol ~ as.factor(rs174548))

## non-parametric ANOVA
kruskal.test(chol ~ as.factor(rs174548))

## One-way (not assuming equal variances)
oneway.test(chol ~ as.factor(rs174548))

## Using robust standard errors
library(gee)
summary(gee(chol ~ as.factor(rs174548), id=seq(1,length(chol))))


## Two-way ANOVA ------------------------------------------------------------
fit1 = lm(chol ~ as.factor(sex) + as.factor(rs174548))
summary(fit1)
anova(fit1)

fit2 = lm(chol ~ as.factor(sex) * as.factor(rs174548))
summary(fit2)
anova(fit2)

anova(fit1, fit2)


## Analysis of Covariance ----------------------------------------------------
fit0 = lm(chol ~ as.factor(rs174548))
summary(fit0)
anova(fit0)

fit1 = lm(chol ~ as.factor(rs174548) + age)
summary(fit1)
anova(fit1)
 
fit2 = lm(chol ~ as.factor(rs174548) * age)
summary(fit2)
anova(fit2)

anova(fit0, fit1, fit2)


## show the fitted lines for the model without an interaction ---------------------------------------
plot(age[rs174548 == 0],chol[rs174548 == 0], xlim = c(min(age),max(age)),
	ylim = c(min(chol),max(chol)), xlab = "Age (years)", ylab = "Total cholesterol (mg/dl)")
# add points for C/G
points(age[rs174548 == 1],chol[rs174548 == 1], col = "red")
# add points for G/G
points(age[rs174548 == 2],chol[rs174548 == 2], col = "blue")

# superimpose regression lines
abline(coef(fit1)[1],coef(fit1)[4])
abline(sum(coef(fit1)[c(1,2)]),coef(fit1)[4], col = "red")
abline(sum(coef(fit1)[c(1,3)]),coef(fit1)[4], col = "blue")
legend(30,240, lty = c(1,1), col = c("black","red","blue"), legend = c("C/C","C/G", "G/G"))


## show the fitted lines for the model with an interaction ------------------------------------------
plot(age[rs174548 == 0],chol[rs174548 == 0], xlim = c(min(age),max(age)),
	ylim = c(min(chol),max(chol)), xlab = "Age (years)", ylab = "Total cholesterol (mg/dl)")
# add points for C/G
points(age[rs174548 == 1],chol[rs174548 == 1], col = "red")
# add points for G/G
points(age[rs174548 == 2],chol[rs174548 == 2], col = "blue")

# superimpose regression lines
abline(coef(fit2)[1],coef(fit2)[4])
abline(sum(coef(fit2)[c(1,2)]),sum(coef(fit2)[c(4,5)]), col = "red")
abline(sum(coef(fit2)[c(1,3)]),sum(coef(fit2)[c(4,6)]), col = "blue")
legend(30,240, lty = c(1,1), col = c("black","red","blue"), legend = c("C/C","C/G", "G/G"))

## mean cholesterol for different genotypes
predict(fit0, new=data.frame(rs174548=0))
predict(fit0, new=data.frame(rs174548=1))
predict(fit0, new=data.frame(rs174548=2))

## mean cholesterol for different genotypes adjusted by age
predict(fit1, new=data.frame(age=mean(age),rs174548=0))
predict(fit1, new=data.frame(age=mean(age),rs174548=1))
predict(fit1, new=data.frame(age=mean(age),rs174548=2))

