library(mice) #imputations
library(readr)

library(magrittr) # piping 
library(dplyr) #data manipulation
library(furrr) #parallel mapping

setwd("C:\Users\mirth\OneDrive\Documenten\ADS\Thesis")
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
cor(datad)
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

#plot the distributions of the observed data 
distributions<-function(x){ 
  hist_preg <- hist(x$Pregnancies)
  hist_gluc <- hist(x$Glucose)
  hist_BMI <- hist(x$BMI)
  return(hist_preg,hist_gluc,hist_BMI)
}

stat_properties(truedata)
distributions(truedata)


#simulation  
set.seed(123) #seed for reproducibility 

#parallel processing for an increased speed
plan(multisession)

nsim=1000 #number of iterations, later change to 1000 

#specify method 
meth <- make.method(truedata)
#meth <- rep("pmm",ncol(truedata))  #OR this??
names(meth) <- colnames(truedata)

#??? Do we need to set Outcome to logistic regression
#meth['Outcome'] <- "logreg"

#cart method, all variables are imputed by means of cart (classification and regression trees)
cart <- rep("cart", ncol(truedata))
names(cart) <- colnames(truedata)

#specify predictor matrix 
pred <- make.predictorMatrix(truedata)
pred
# ??? passive imputation needed

#default synthetic datasets
syn_ds <- future_map(1:nsim, ~ { 
  truedata %>% mice(m=5,
                    method = meth,
                    predictorMatrix = pred,
                    where=matrix(TRUE, nrow(truedata), ncol(truedata)),
                    print=FALSE)
}, .options=future_options(seed=as.integer(123)), .progress=TRUE, .id = "syn")

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



# some nodes for further steps
simulation<-function(n){
  
}

result <- replicate(nsim, simulation(n=nsim), simplify = FALSE)









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




