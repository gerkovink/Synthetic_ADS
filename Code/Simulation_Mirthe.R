#Mirthe Hendriks 
# 19-5-2021

library(mice) #imputations
library(readr)

library(magrittr) # piping 
library(dplyr) #data manipulation
library(furrr) #parallel mapping
library(corrplot)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(broom)
library(purrr)
library(caret)
library(mboost) # for prediction model with method glmboost for classification
library(DT)

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

truedata <- truedata[!truedata$BMI == 0,] #disregard data with BMI of 0 
truedata <- truedata[!truedata$BloodPressure == 0,] #disregard data with BloodPressure of 0 

#looking at the data 
summary(truedata)
dim(truedata)
str(truedata)

#complete model; logistic regression, all variables included
model1<-glm(Outcome~Pregnancies+Glucose+BloodPressure+SkinThickness+Insulin+BMI+DiabetesPedigreeFunction+Age, data=truedata, family = "binomial")
summary(model1)

### Formulate an analysis model 

#model: logistic regression, three selected predictors 
model_true<-glm(Outcome~BMI+Glucose+Pregnancies, data=truedata, family=binomial(link = "logit"))
summary(model_true)
confint(model_true)

#logistic regression model shortened 
model<-function(x){ 
  glm(Outcome~BMI+Glucose+Pregnancies, data=x, family=binomial(link = "logit"))
}
model(truedata)


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
truedata_num<- data.frame(lapply(truedata, function(x) as.numeric(as.character(x))))
cor(truedata_num)
corrplot(cor(truedata_num))

distributions <- function(x){
  ggplot(gather(x),aes(value))+
    geom_histogram(bins=10)+
    facet_wrap(~key, scales = 'free_x')
}


#plot the distributions of the observed data 
distributions_model<-function(x){ 
  hist_preg <- hist(x$Pregnancies)
  hist_gluc <- hist(x$Glucose)
  hist_BMI <- hist(x$BMI)
}

stat_properties(truedata)
distributions(datad)
distributions_model(truedata)


### SIMULATION 

#parallel processing for an increased speed
plan(multisession)

nsim=1000 #number of iterations, later change to 1000 

#specify method 
meth <- make.method(truedata)
meth <- rep("pmm",ncol(truedata))  
names(meth) <- colnames(truedata)
meth['Outcome'] <- "logreg"

#cart method, all variables are imputed by means of cart (classification and regression trees)
cart <- rep("cart", ncol(truedata))
names(cart) <- colnames(truedata)

#specify predictor matrix 
pred <- make.predictorMatrix(truedata)

#default synthetic datasets; method specified as pmm and logreg for the binary outcome variable
syn_ds <- future_map(1:nsim, ~ { 
  truedata %>% mice(m=5,
                    method = meth,
                    predictorMatrix = pred,
                    where=matrix(TRUE, nrow(truedata), ncol(truedata)),
                    print=FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=TRUE, .id = "syn")
#set seed again in future_options to ensure reproducable results with parallel processing 

syn_ds_10 <- future_map(1:10, ~ { 
  truedata %>% mice(m=5,
                    method = meth,
                    predictorMatrix = pred,
                    where=matrix(TRUE, nrow(truedata), ncol(truedata)),
                    print=FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=TRUE, .id = "syn")


#make synthetic datasets using cart 
syn_cart <- future_map(1:nsim, ~ { 
  truedata %>% mice(m=5,
                    method = cart,
                    predictorMatrix = pred,
                    where=matrix(TRUE, nrow(truedata), ncol(truedata)),
                    print=FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=T, .id = "syn")

syn_cart_maxit <- future_map(1:nsim, ~ { 
  truedata %>% mice(m = 5,
                    maxit = 1,
                    method = cart,
                    predictorMatrix = pred,
                    where = matrix(TRUE, nrow(truedata), ncol(truedata)),
                    print = FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=T, .id = "syn")

syn_cartcp <- future_map(1:nsim, ~ { 
  truedata %>% mice(m = 5,
                    method = cart,
                    predictorMatrix = pred,
                    where = matrix(TRUE, nrow(truedata), ncol(truedata)),
                    cp = 1e-32,
                    minbucket = 3,
                    print = FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=T, .id = "syn")

syn_cartcp_maxit <- future_map(1:nsim, ~ { 
  truedata %>% mice(m = 5,
                    maxit = 1,
                    method = cart,
                    predictorMatrix = pred,
                    where = matrix(TRUE, nrow(truedata), ncol(truedata)),
                    cp = 1e-32,
                    minbucket = 3,
                    print = FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=T, .id = "syn")

syn_ds <- readRDS('syn_ds_model.rds')
syn_cart <- readRDS('syn_cartt_model.rds')
syn_cart_maxit <- readRDS('syn_cart_maxit.rds')
syn_cartcp <- readRDS('syn_cartcp.rds')
syn_cartcp_maxit <- readRDS('syn_cartcp_maxit.rds')

syn_ds[1]
syn_cart[1]
syn_cart_maxit[1]
syn_cartcp[1]
syn_cartcp_maxit[1]

### Pooling 
# pool function partially synthetic data
pool3.syn <- function(mira) {
  if(class(mira)[1] == "mira") { # if the input object is of class mira
    fitlist <- mira %$% analyses # extract the analyses from the mira object
  } 
  else {
    fitlist <- mira              # and otherwise, just take the input list
  }
  
  vars <- fitlist[[1]] %>% coef() %>% names()

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
              upper   = est + qt(.975, df) * sqrt(var), .groups = 'drop') %>% 
    arrange(factor(term, levels=vars)) #order the rows by values of selected columns 
  pooled
}

#pools the m=5 iterations for the 1000 simulations 


# Pooling various synthetic datasets

#outputs a tibble of 4000 x 8 with the pooled parameter estimates for the 1000 simulations;
# an intercept and three regression coefficients (BMI, Glucose, Pregnancies)

pooled_ds <- map_dfr(syn_ds, function(x) {
  x %$%
    glm(Outcome ~ BMI + Glucose + Pregnancies, family = binomial(link = "logit")) %>% 
    pool3.syn()
})

pooled_cart <- map_dfr(syn_cart, function(x) {
  x %$%
    glm(Outcome ~ BMI + Glucose + Pregnancies, family = binomial(link = "logit")) %>% 
    pool3.syn()
})

pooled_cart_maxit <- map_dfr(syn_cart_maxit, function(x) {
  x %$%
    glm(Outcome ~ BMI + Glucose + Pregnancies, family = binomial(link = "logit")) %>% 
    pool3.syn()
})


pooled_cartcp <- map_dfr(syn_cartcp, function(x) {
  x %$%
    glm(Outcome ~ BMI + Glucose + Pregnancies, family = binomial(link = "logit")) %>% 
    pool3.syn()
})

pooled_cartcp_maxit <- map_dfr(syn_cartcp_maxit, function(x) {
  x %$%
    glm(Outcome ~ BMI + Glucose + Pregnancies, family = binomial(link = "logit")) %>% 
    pool3.syn()
})

pooled_ds
pooled_cart
pooled_cart_maxit
pooled_cartcp
pooled_cartcp_maxit

#compute confidence interval coverage 
ci_cov <- function(pooled, true_fit = NULL, coefs = NULL, vars = NULL) {
  
  if (!is.null(true_fit)) {
    coefs <- coef(true_fit)
    vars   <- diag(vcov(true_fit))
  }
 
  nsim <- nrow(pooled) / length(unique(pooled$term))
  
  pooled %>% mutate(true_coef = rep(coefs, times = nsim),
                    true_var  = rep(vars, times = nsim),
                    cover     = lower < true_coef & true_coef < upper) %>%
    group_by(term) %>%
    summarise("True Est" = unique(true_coef),
              "Syn Est"  = mean(est),
              "Bias"     = mean(est - true_coef),
              "True SE"  = unique(sqrt(true_var)),
              "Syn SE"   = mean(sqrt(var)),
              "df"       = mean(df),
              "Lower"    = mean(lower),
              "Upper"    = mean(upper),
              "CIW"      = mean(upper - lower),
              "Coverage" = mean(cover), .groups = "drop")
}

cic_ds <- ci_cov(pooled_ds, true_fit = model_true)
cic_cart <- ci_cov(pooled_cart, true_fit = model_true)
cic_cart_maxit <- ci_cov(pooled_cart_maxit, true_fit = model_true)
cic_cartcp <- ci_cov(pooled_cartcp, true_fit = model_true)
cic_cartcp_maxit <- ci_cov(pooled_cartcp_maxit, true_fit = model_true)

column_names= c('True Est', 'Syn Est', 'Bias', 'True SE','Syn SE', 'df', 'Lower', 'Upper', 'CIW', 'Coverage')

datatable(cic_ds) %>% formatRound(columns = column_names,digits=3)
datatable(cic_cart) %>% formatRound(columns = column_names,digits=3)
datatable(cic_cart_maxit) %>% formatRound(columns = column_names,digits=3)
datatable(cic_cartcp) %>% formatRound(columns = column_names,digits=3)
datatable(cic_cartcp_maxit) %>% formatRound(columns = column_names,digits=3)


ggplot(data = pooled_cart, aes(est)) + geom_histogram(bins=25) + facet_wrap(~term, scales = 'free_x')


#predict: synthetic or true data 

## Prediction function
predict_synth_or_true = function(data, model) {
  set.seed(123) #for reproducibility
  predicted_values = predict(object = model, newdata = data)
  #true_values = data[,ncol(data)]
  
  true_value = c()
  for (i in data$Synth_or_true) {
    true_value = c(true_value, i)
  }
  
  true_value = factor(true_value, levels = c(0, 1))
  
  result = data.frame(true_value, predicted_values = factor(predicted_values, levels = c(0, 1)))
  return(result)
}

#create train and test datasets and predict 

total_result = data.frame() #data frame in which the results of all sampled data will be stored
for (i in seq(1, length(syn_cartcp_maxit), 1)) { #looping through each simulation
  imputation = syn_cartcp_maxit[[i]]
  for (n in seq(1, imputation$m, 1)) { #looping through each imputation of the simulation
    syndata = complete(imputation, action = n)
    
    true_data = sample_n(truedata, (nrow(truedata)*0.5)) #sampling half of the original data
    true_data$Synth_or_true = 0 #indicator: setting the true data to 0
    
    synthetic_data = sample_n(syndata, (nrow(syndata)*0.5)) #sampling half of the synthetic data
    synthetic_data$Synth_or_true = 1 #setting the synthetic data to 1
    
    prediction_data = rbind(true_data, synthetic_data) #combining the two samples into a single data set
    prediction_data$Synth_or_true = factor(prediction_data$Synth_or_true, levels = c(0, 1))
    
    train_index = sample(seq_len(nrow(prediction_data)), size = floor(0.30 * nrow(prediction_data))) #creating a list of indices for training data
    train = prediction_data[train_index,]
    test = prediction_data[-train_index,]
    
    prediction_model = train(Synth_or_true ~ ., data = train, method = "glm") #prediction model with method glm for classification
    result = predict_synth_or_true(data=test, model=prediction_model) #use the prediction function to get a data frame of true and predicted values
    total_result = rbind(total_result, result) #appending the data frame to the overarching data frame
  }
}
  cfmatrix = confusionMatrix(total_result$predicted_value, total_result$true_value) #creating a confusion matrix
print(cfmatrix)
print(prediction_model)



