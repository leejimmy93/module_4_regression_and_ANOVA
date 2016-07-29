## R commands for linear regression analysis of cholesterol data set
  
## read data set ----------------------------------------------------------------------
cholesterol = read.table("http://faculty.washington.edu/rhubb/sisg/SISG-Data-cholesterol.txt", header=T)
attach(cholesterol)

## visualize a few rows of data--------------------------------------------------------
head(cholesterol)

# scatterplot of total cholesterol vs age
plot(age,chol, xlab = "Age (years)", ylab = "Total cholesterol (mg/dl)")
# with regression line superimposed
abline(lm(chol~age)$coef, col = "red")

## dichotomize age and compare cholesterol in two groups-------------------------------
group = 1*(age > 55)t.test(chol ~ group)
# boxplot of cholesterol by age group
boxplot(chol ~ group, names = c("30 - 55", "56 - 80"), xlab = "Age", ylab = "Cholesterol (mg/dl)")

## fit a linear model to the age and cholesterol data----------------------------------
fit = lm(chol ~ age)summary(fit)confint(fit)

## prediction intervals-----------------------------------------------------------------
predict.lm(fit, newdata=data.frame(age=c(46,47,48)), interval="confidence")predict.lm(fit, newdata=data.frame(age=c(46,47,48)), interval="prediction")


## anova table--------------------------------------------------------------------------
anova(fit)

## residuals plots----------------------------------------------------------------------
# residuals vs fitted values
plot(fit$fitted, fit$residuals)
# horizontal reference line
abline(0,0)

# Qunatile-quantile plot of residuals
qqnorm(fit$residuals, main = "")

## robust regression--------------------------------------------------------------------
# plot triglycerides vs age
plot(age, TG)
abline(lm(TG~age)$coef, col = "red")

# fit simple linear regression for age and TG
fit.tg = lm(TG ~ age)

# plot residuals vs fitted values
plot(fit.tg$fitted, fit.tg$residuals)
# horizontal reference line
abline(0,0)

# fit robus regression
library(gee)
fit.tg.ese = gee(TG ~ age, id = seq(1, length(age)))

## transformations----------------------------------------------------------------------
#plot log transformed TG vs age
plot(age,log(TG))
abline(lm(log(TG)~age)$coef, col = "red")

# fit regression model for log transformed outcomes
fit.tg.ln = lm(log(TG) ~ age)

# plot residuals from log transformed model
plot(fit.tg.ln$fitted, fit.tg.ln$residuals)
# horizontal reference line
abline(0,0)

## investigate delta betas--------------------------------------------------------------
dfb = dfbeta(fit)

## multiple regression for effect of age and BMI on cholesterol-------------------------
par(mfrow = c(1,2))
plot(age,BMI)
plot(BMI,chol)

fit2 = lm(chol ~ age + BMI)
summary(fit2)

# compare age only and age,birthweight models
anova(fit2)
anova(fit,fit2)

## effects modification by sex?---------------------------------------------------------
plot(age[sex == 0],chol[sex==0], xlab = "age", ylab = "chol")
points(age[sex == 1],chol[sex==1], col = "red")
legend(65,140, lty = c(1,1), col = c("black","red"), legend = c("Male","Female"))

# first fit model with age and sex main effects
fit2 = lm(chol ~ age+sex)

# next add sex interaction
fit3 = lm(chol ~ age*sex)
anova(fit3)
anova(fit2,fit3)

## plot effects separately by gender----------------------------------------------------
# scatterplot with points for males only
plot(age[sex == 0],chol[sex == 0], xlim = c(min(age),max(age)),
	ylim = c(min(chol),max(chol)), xlab = "Age (years)", ylab = "Total cholesterol (mg/dl)")
# add points for females
points(age[sex == 1],chol[sex == 1], col = "red")
# superimpose regression line for males
abline(coef(fit3)[1],coef(fit3)[2])
# and regression line for females
abline(sum(coef(fit3)[c(1,3)]),sum(coef(fit3)[c(2,4)]), col = "red")
legend(30,240, lty = c(1,1), col = c("black","red"), legend = c("Male","Female"))
