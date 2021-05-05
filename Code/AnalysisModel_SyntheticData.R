#Mirthe Hendriks 
# 5-5-2021

library(mice) #imputations
library(readr)

library(magrittr) # piping 
library(dplyr) #data manipulation
library(furrr) #parallel mapping
library(corrplot)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(purrr)

library(broom.mixed) #pipe for mira? 

#simulation  
set.seed(123) #seed for reproducibility 

setwd("~/ADS/Thesis")
truedata <- read_csv("diabetes.csv", 
                 col_types = cols(Pregnancies = col_integer(), 
                                  Glucose = col_integer(), 
                                  BloodPressure = col_integer(), 
                                  SkinThickness = col_integer(), 
                                  Insulin = col_integer(), 
                                  Age = col_integer(), 
                                  Outcome = col_factor(levels = c("1", "0"))))

datad <- read_csv("diabetes.csv")

#looking at the data 
summary(truedata)
dim(truedata)
str(truedata)

#complete model; logistic regression, all variables included
model1<-glm(Outcome~Pregnancies+Glucose+BloodPressure+SkinThickness+Insulin+BMI+DiabetesPedigreeFunction+Age, data=truedata, family = "binomial")
summary(model1)

### Formulate an analysis model 

#model: logistic regression, three selected predictors 
model_true<-glm(Outcome~Pregnancies+Glucose+BMI, data=truedata, family="binomial")
summary(model_true)
confint(model_true)

#logistic regression model shortened 
model<-function(x){ 
  glm(Outcome~Pregnancies+Glucose+BMI, data=x, family="binomial")
}
model(truedata)

hist(truedata$Outcome)

#obtain the statistical properties of the observed data 
stat_properties<-function(x){ 
  mu_preg<-mean(x$Pregnancies)
  mu_gluc<-mean(x$Glucose)
  mu_BMI<-mean(x$BMI)
  var_preg <- var(x$Pregnancies)
  var_gluc <- var(x$Glucose)
  var_BMI <- var(x$BMI)
  beta_preg <- x %>% model %>% coefficients %>% .[2]
  beta_gluc <- x %>% model %>% coefficients %>% .[3]
  beta_BMI <- x %>% model %>% coefficients %>% .[4]
  resvar <- x %>% model %>% residuals %>% var
  R2 <- x %>% model %>% summary %>% .$r.squared
  return(unlist(c(mu_preg,mu_gluc,mu_BMI,var_preg,var_gluc,var_BMI,beta_preg,beta_gluc,beta_BMI,resvar,R2)))
}

#correlation variables 
cor(datad)
corrplot(cor(datad))

distributions <- function(x){
  ggplot(gather(x),aes(value))+
    geom_histogram(bins=10)+
    facet_wrap(~key, scales = 'free_x')
}

?ggplot

#plot the distributions of the observed data 
distributions_model<-function(x){ 
  hist_preg <- hist(x$Pregnancies)
  hist_gluc <- hist(x$Glucose)
  hist_BMI <- hist(x$BMI)
  return(hist_preg,hist_gluc,hist_BMI)
}

stat_properties(truedata)
distributions(datad)
distributions_model(truedata)

#parallel processing for an increased speed
plan(multisession)

nsim=1000 #number of iterations, later change to 1000 

#specify method 
meth <- make.method(truedata)
meth <- rep("pmm",ncol(truedata))  
names(meth) <- colnames(truedata)
meth['Outcome'] <- "logreg"
#meth['DiabetesPedigreeFunction'] <- "polr"

#cart method, all variables are imputed by means of cart (classification and regression trees)
cart <- rep("cart", ncol(truedata))
names(cart) <- colnames(truedata)

#specify predictor matrix 
pred <- make.predictorMatrix(truedata)

#default synthetic datasets
syn_ds <- future_map(1:nsim, ~ { 
  truedata %>% mice(m=5,
                    method = meth,
                    predictorMatrix = pred,
                    where=matrix(TRUE, nrow(truedata), ncol(truedata)),
                    print=FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=TRUE, .id = "syn")

syn_ds[1]

# in the output the imputation methods are empty "" for all variables, 
# these empty strings imply that the variable is complete and therefore not imputed
# clearly something is has gone wrong here

#make synthetic datasets using cart 
syn_cart <- future_map(1:nsim, ~ { 
  truedata %>% mice(m=5,
                    method = cart,
                    predictorMatrix = pred,
                    where=matrix(TRUE, nrow(truedata), ncol(truedata)),
                    print=FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=T, .id = "syn")


syn_cart[2]

saveRDS(syn_cart, 'syn_cart_model.rds')
syn_cart <- readRDS('syn_cart_model.rds')

# pool function partially synthetic data
pool3.syn <- function(mira) {
  if(class(mira)[1] == "mira") { # if the input object is of class mira
    fitlist <- mira %$% analyses # extract the analyses from the mira object
  } 
  else {
    fitlist <- mira              # and otherwise, just take the input list
  }
  
  m <- length(fitlist)           # number of imputations
  
  pooled <- fitlist %>% 
    map_dfr(broom::tidy) %>%
    group_by(term) %>%
    summarise(est     = mean(estimate),
              bm      = sum((estimate - est)^2) / (m - 1),
              ubar    = mean(std.error^2),
              var     = ubar + bm/m, # new variance estimate
              df      = (m - 1) * (1 + (ubar * m)/bm), # and new df estimate
              lower   = est - qt(.975, df) * sqrt(var),
              upper   = est + qt(.975, df) * sqrt(var), .groups = 'drop')
  pooled
}

map_dfr(syn_ds, function(x) {
  x %$%
    glm(Outcome ~ Pregnancies + Glucose + BMI, family = binomial) %>% 
    pool3.syn()
})

# this outputs a tibble of 4000 x 8, so not really a pooled output yet 

map_dfr(syn_cart, function(x) {
  x %$%
    glm(Outcome ~ Pregnancies + Glucose + BMI, family = binomial) %>% 
    pool3.syn()
})

syn_pooled <- map_dfr(syn_ds, function(x) {
  x %$%
    glm(Outcome ~ Pregnancies + Glucose + BMI, family = binomial) %>% 
    pool3.syn()
})


pool3.syn(syn_ds)

pool3.syn(syn_cart)


## does not work yet 
# get output
syn_out <- syn_pooled %>%
  map(function(x) summary(x, population.inference = TRUE)) %>%
  map_dfr(function(x) {
    coef(x) %>% 
      as.data.frame %>%
      rownames_to_column(var = "term")})

syn_out




#other stats properties function which includes the model 
stat_properties_M<-function(x){ 
  mu_preg<-mean(x$Pregnancies)
  mu_gluc<-mean(x$Glucose)
  mu_BMI<-mean(x$BMI)
  var_preg <- var(x$Pregnancies)
  var_gluc <- var(x$Glucose)
  var_BMI <- var(x$BMI)
  beta_preg <- glm(Outcome~Pregnancies+Glucose+BMI, data=x, family="binomial") %>% coefficients %>% .[2]
  beta_gluc <- glm(Outcome~Pregnancies+Glucose+BMI, data=x, family="binomial") %>% coefficients %>% .[3]
  beta_BMI <- glm(Outcome~Pregnancies+Glucose+BMI, data=x, family="binomial") %>% coefficients %>% .[4]
  resvar <- beta_preg <- glm(Outcome~Pregnancies+Glucose+BMI, data=x, family="binomial") %>% residuals %>% var
  R2 <- beta_preg <- glm(Outcome~Pregnancies+Glucose+BMI, data=x, family="binomial") %>% summary %>% .$r.squared
  return(unlist(c(mu_preg,mu_gluc,mu_BMI,var_preg,var_gluc,var_BMI,beta_preg,beta_gluc,beta_BMI,resvar,R2)))
}




