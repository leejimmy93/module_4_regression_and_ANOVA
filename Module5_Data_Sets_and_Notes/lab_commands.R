###### LAB 1 ######################################################################
## read data set
cholesterol = read.table("http://faculty.washington.edu/rhubb/sisg/SISG-Data-cholesterol.txt", header=T)
attach(cholesterol) 

# compute univariate summary statistics for triglycerides
mean(TG)
median(TG)
sd(TG)
summary(TG)

# graphical displays for triglycerides
boxplot(TG)
stem(TG)

hist(TG)
 
# create a binary indicator for BMI > 25
ibmi = ifelse(BMI > 25, 1, 0)

# compute univariate summary statistics for triglycerides for BMI > 25 and BMI <= 25
tapply(TG,ibmi,mean)
tapply(TG,ibmi,median)
tapply(TG,ibmi,sd)

# plot boxplots for triglycerides separately by BMI > 25 and BMI <= 25
boxplot(TG ~ ibmi)
t.test(TG ~ ibmi)

# scatterplot of triglycerides vs BMI
plot(BMI, TG)

# fit linear regression models for the association between triglycerides and BMI
fit1 = lm(TG ~ BMI)
summary(fit1)

# get prediction intervals at BMI = 23
predict(fit1, newdata = data.frame(BMI = 23), interval = "confidence")
predict(fit1, newdata = data.frame(BMI = 23), interval = "prediction")

###### LAB 2 ######################################################################
# load the gee() package for robust standard errors
library(gee)
 
# identify outliers in scatterplot of triglycerides vs BMI
plot(BMI,TG)
bmi37 = which(BMI<=37)

# excluding subjects with BMI > 37
fit2 = lm(TG[bmi37] ~ BMI[bmi37])
summary(fit2)

# analyze residuals from the regression analysis of triglycerides and BMI
plot(fit2$fitted, fit2$residuals)
abline(0,0)
qqnorm(fit2$residuals)
qqline(fit2$residuals)

# fit a linear regression model with robust standard errors
fit.gee = gee(TG ~ BMI, id = seq(1,length(TG)))
summary(fit.gee)

# calculate p-values for robust regression
z = abs(fit.gee$coef/sqrt(diag(fit.gee$robust)))
2*(1-pnorm(z))

# fit a regression model for log transformed triglycerides and BMI
fit.log = lm(log(TG) ~ BMI)
summary(fit.log)

# analyze residuals from the regression analysis of log transformed 
# triglycerides and BMI
par(mfrow = c(1,2))
plot(fit.log$fitted, fit.log$residuals)
abline(0,0)
qqnorm(fit.log$residuals)
qqline(fit.log$residuals)

# binary variable indicating presence of ApoE4
APOE4 = ifelse(APOE %in% c(3,5,6), 1, 0)

# scatterplot with subjects stratified by ApoE4
par(mfrow = c(1,1))
plot(BMI[APOE4 == 0], TG[APOE4 == 0], pch = 19, xlab = "BMI", ylab = "triglycerides")
points(BMI[APOE4 == 1], TG[APOE4 == 1], pch = 1)
 
# multiple linear regression of triglycerides on BMI, ApoE4, and interaction
fit3 = lm(TG ~ BMI + APOE4 + BMI*APOE4)
summary(fit3)


###### LAB 3 ######################################################################
library(gee)
library(multcomp)
 
## read data set ----------------------------------------------------------------
cholesterol = read.table("http://faculty.washington.edu/rhubb/sisg/SISG-Data-cholesterol.txt", header=T)
attach(cholesterol)
 
## Exploratory data analysis ----------------------------------------------------
## graphical display: boxplot 
par(mfrow = c(1,2))
boxplot(chol ~ as.factor(rs4775401))
boxplot(chol ~ as.factor(APOE))

## alternative graphical display: graph of means 
par(mfrow = c(2,1))
plot.design(chol ~ as.factor(rs4775401))
plot.design(chol ~ as.factor(APOE))


## numeric descriptives 
tapply(chol, as.factor(rs4775401), mean)
tapply(chol, as.factor(rs4775401), sd)
tapply(chol, as.factor(APOE), mean)
tapply(chol, as.factor(APOE), sd)

## Inferential data analysis ----------------------------------------------------
fit1 = lm(chol ~ as.factor(rs4775401))
summary(fit1)
anova(fit1)

fit2 = lm(chol ~ as.factor(APOE))
summary(fit2)
anova(fit2)

M2 = contrMat(table(APOE), type="Tukey")
fit3 = lm(chol ~ -1 + as.factor(APOE))
mc2 = glht(fit3, linfct =M2)
summary(mc2, test=adjusted("none"))
summary(mc2, test=adjusted("bonferroni"))
summary(mc2, test=adjusted("holm"))
summary(mc2, test=adjusted("hochberg"))
summary(mc2, test=adjusted("hommel"))
summary(mc2, test=adjusted("BH"))
summary(mc2, test=adjusted("BY"))
summary(mc2, test=adjusted("fdr"))


## One-way (not assuming equal variances)
oneway.test(chol ~ as.factor(rs4775401))
oneway.test(chol ~ as.factor(APOE))

## Using robust standard errors
summary(gee(chol ~ as.factor(rs4775401), id=seq(1,length(chol))))
summary(gee(chol ~ as.factor(APOE), id=seq(1,length(chol))))

## non-parametric ANOVA
kruskal.test(chol ~ as.factor(rs4775401))
kruskal.test(chol ~ as.factor(APOE))


###### LAB 4 ######################################################################
## Two-way ANOVA ------------------------------------------------------------
## exploratory data analysis
table(rs174548, APOE)
tapply(chol, list(as.factor(rs174548), as.factor(APOE)), mean)
tapply(chol, list(as.factor(rs174548), as.factor(APOE)), sd)
plot.design(chol ~ as.factor(rs174548) + as.factor(APOE))

## model with interaction
fit1 = lm(chol ~ as.factor(rs174548)*as.factor(APOE))
summary(fit1)
 
## model without interaction
fit2 = lm(chol ~ as.factor(rs174548) + as.factor(APOE))
summary(fit2)

## compare models with and without interaction
anova(fit2,fit1)

###### LAB 5 ######################################################################
by(cbind(chol,age), APOE, cor, method="pearson")
by(cbind(chol,age), APOE, cor, method="spearman")

plot(age, chol, xlab="AGE (yrs)", ylab="CHOLESTEROL (mg/dl)", type="n")
for (i in 1:6){
  lines(lowess(age[APOE==i], chol[APOE==i]), col=i)
  points(age[APOE==i], chol[APOE==i], col=i, pch=16)
 }
legend(min(age), max(chol), legend=paste("APOE", seq(1,6)), col=seq(1,6), pch=16, lty=1)

## ANCOVA Model with an interaction
fit1 = lm(chol ~ as.factor(APOE) * age)
summary(fit1)

## ANCOVA Model without an interaction
fit2 = lm(chol ~ as.factor(APOE) + age)
summary(fit2)

## compare models with and without interaction
anova(fit2, fit1)

## ONE-WAY ANOVA model 
fit3 = lm(chol ~ as.factor(APOE))
summary(fit3)
anova(fit3, fit2)

## mean cholesterol for different genotypes
predict(fit3, new=data.frame(APOE=1))
predict(fit3, new=data.frame(APOE=2))
predict(fit3, new=data.frame(APOE=3))
predict(fit3, new=data.frame(APOE=4))
predict(fit3, new=data.frame(APOE=5))
predict(fit3, new=data.frame(APOE=6))

## mean cholesterol for different genotypes adjusted by age
predict(fit2, new=data.frame(age=mean(age),APOE=1))
predict(fit2, new=data.frame(age=mean(age),APOE=2))
predict(fit2, new=data.frame(age=mean(age),APOE=3))
predict(fit2, new=data.frame(age=mean(age),APOE=4))
predict(fit2, new=data.frame(age=mean(age),APOE=5))
predict(fit2, new=data.frame(age=mean(age),APOE=6))
