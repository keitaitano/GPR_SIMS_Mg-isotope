# title: Correction of IMF during SIMS Mg isotope analysis
# author: Keita Itano 
# date:2024/08/13

# Library
library(tidyverse)
library(kernlab)
library(plotly)

# data import
df <- read.csv("Fukuda2020_data.csv") #%>% na.omit()
df_O2 <- df %>%  select(Label, session, bias, bias_error, R_Mg_Si, Yield, MgO:MnO,Fo_Num) %>% mutate(X1=R_Mg_Si/MgO, X2=Yield/MgO)

# Feature selection
df_O2 <- df_O2  %>% mutate(Label=factor(Label))
x <- df_O2 %>% mutate(R_FeMg=(FeO/MgO), R_CaMg=(CaO/MgO), R_CrMg=(Cr2O3/MgO), R_MnMg=(MnO/MgO)) %>% select(R_FeMg, R_CaMg, R_CrMg,R_MnMg)
y <- df_O2$bias
measurement_errors <- df_O2$bias_error
measurement_variances <- measurement_errors^2

##### hyper parameter tuning  ####
## Group K-fold 
K <- length(levels(df_O2$Label))
group_folds  <- list()
for (k in 1:K) {
  fold  <- which(df_O2$Label == levels(df_O2$Label)[k])
  group_folds[[length(group_folds) + 1]] <- fold
}
names(group_folds) <- paste0("Fold", formatC(1:K,width=2,flag="0"))

## model training
sigma <- seq(0.01, 1, length.out = 100)
hyperparam_compar <-  data.frame(matrix(ncol = 3, nrow = 0))

for (sigma in sigma) {
  CV_result <- data.frame(matrix(ncol = 2, nrow = 0))
  for (k in 1:K) {
    noise <- measurement_variances[-group_folds[[k]]] %>% mean(.)/sd(y[-group_folds[[k]]])
    
    fit_rbf <- gausspr(x[-group_folds[[k]],],
                          y[-group_folds[[k]]],
                          variance.model=T,
                          var=noise,
                          scaled=T, 
                          kernel="rbfdot",
                          kpar=list(sigma=sigma)
    )
    
    
    summary_pred <-data.frame(true = y[group_folds[[k]]],
                              pred= predict(fit_rbf, x[group_folds[[k]],]),
                              sd =  predict(fit_rbf, x[group_folds[[k]],], type="sdeviation")) %>% 
                  mutate(residual_bias = true-pred) %>% 
                  mutate(residual_bias_squared = residual_bias^2)
    summary_pred_train <- data.frame(true = y[-group_folds[[k]]],
                                     pred= predict(fit_rbf, x[-group_folds[[k]],]),
                                     sd =  predict(fit_rbf, x[-group_folds[[k]],], type="sdeviation")) %>% 
                          mutate(residual_bias = true-pred) %>% 
                          mutate(residual_bias_squared = residual_bias^2)
    
    #rmse_train <- fit_rbf@error
    rmse_train<- summary_pred_train$residual_bias_squared %>% sum()/length(y[-group_folds[[k]]]) %>% .^0.5
    rmse_test <- summary_pred$residual_bias_squared %>% sum()/length(y[group_folds[[k]]]) %>% .^0.5
    CV_result <- rbind(CV_result, data.frame(rmse_train, rmse_test))
  }
  
  hyperparam_compar <- CV_result %>% apply(., 2, mean) %>%  matrix(., ncol = 2) %>% data.frame() %>%
    mutate(sigma=sigma, RMSE_train=X1, RMSE_test=X2) %>% select(-X1, -X2) %>% rbind(hyperparam_compar,.)  
}

  
hyperparam_compar %>% mutate(length_scale=seq(1, 100, length.out = 100)) %>% 
  ggplot()+
  geom_line(aes(x=sigma, y=RMSE_train),col="red")+
  geom_line(aes(x=sigma, y=RMSE_test),col="blue")+
  scale_x_log10()+
  scale_y_continuous(limits = c(0, 1.6))+
  theme_bw()+
  ylab("RMSE")+
  theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),axis.title.x = element_text(size = 9),axis.title.y = element_text(size = 9))

sigma_opt <- hyperparam_compar %>% filter(RMSE_test == min(hyperparam_compar$RMSE_test)) %>% select(sigma) %>% unlist()


##### sigma = sigma_opt, vaiables: FeO/MgO, CaO/MgO, Cr2O3/MgO, MnO/MgO ####

CV_result2 <- data.frame(matrix(ncol = 6, nrow = 0))
for (k in 1:K) {
  noise <- measurement_variances[-group_folds[[k]]] %>% mean(.)/sd(y[-group_folds[[k]]])
  fit_rbf <- gausspr(x[-group_folds[[k]],],
                     y[-group_folds[[k]]],
                     variance.model=T,
                     var=noise,
                     scaled=T, 
                     kernel="rbfdot",
                     kpar=list(sigma=sigma_opt)
  )
  
  summary_pred_test <-data.frame(true = y[group_folds[[k]]],
                                 pred= predict(fit_rbf, x[group_folds[[k]],]),
                                 sd =  predict(fit_rbf, x[group_folds[[k]],], type="sdeviation")) %>% 
    mutate(residual_bias = true-pred) %>% 
    mutate(residual_bias_squared = residual_bias^2)
  
  summary_pred_train <-data.frame(true = y[-group_folds[[k]]],
                                  pred= predict(fit_rbf, x[-group_folds[[k]],]),
                                  sd =  predict(fit_rbf, x[-group_folds[[k]],], type="sdeviation")) %>% 
    mutate(residual_bias = true-pred) %>% 
    mutate(residual_bias_squared = residual_bias^2)
  
  rmse_train <- summary_pred_train$residual_bias %>% .^2 %>% sum(.)/length(y[-group_folds[[k]]]) %>% sqrt()
  rmse_test <- summary_pred_test$residual_bias_squared %>% sum()/length(y[group_folds[[k]]]) %>% .^0.5
  R_squared_train <- summary_pred_train %>% mutate(ave=true-mean(true)) %>% mutate(A=residual_bias^2, B=ave^2) %>% select(A,B) %>% apply(2, sum) %>% matrix(.,1,2) %>% data.frame() %>% mutate(TEST =1-X1/X2) %>% .[1,3]
  R_squared_test <- summary_pred_test %>% mutate(ave=true-mean(true)) %>% mutate(A=residual_bias^2, B=ave^2) %>% select(A,B) %>% apply(2, sum) %>% matrix(.,1,2) %>% data.frame() %>% mutate(TEST =1-X1/X2) %>% .[1,3]
  residuals <- summary_pred_train$residual_bias
  sigma2 <- var(summary_pred_train$residual_bias)
  log_likelihood <- -0.5 * length(residuals) * log(2 * pi * sigma2) - sum(residuals^2) / (2 * sigma2)
  aic_train <- -2 * log_likelihood + 2 * ncol(x)
  residuals <- summary_pred_test$residual_bias
  sigma2 <- var(summary_pred_test$residual_bias)
  log_likelihood <- -0.5 * length(residuals) * log(2 * pi * sigma2) - sum(residuals^2) / (2 * sigma2)
  aic_test <- -2 * log_likelihood + 2 * ncol(x)
  CV_result2 <- rbind(CV_result2, data.frame(rmse_train, rmse_test,R_squared_train,R_squared_test, aic_train, aic_test))
}

CV_result2 %>% apply(., 2, mean)　%>% format(., digits=2)


### using all data ####
noise <- measurement_variances %>% mean(.)/sd(y)
fit_rbf2 <- gausspr(x,
                   y,
                   variance.model=T,
                   var= noise,
                   scaled=T, 
                   kernel="rbfdot",
                   kpar=list(sigma=sigma_opt))

summary_pred_add <- 
  data.frame(true=df_O2$bias,
             pred= predict(fit_rbf2, x),
             sd =  predict(fit_rbf2, x, type="sdeviation")) %>% 
  mutate(residual_bias = true-pred,Fo=df_O2$Fo_Num)

# prediction performance
RMSE <- summary_pred_add$residual_bias %>% .^2 %>% sum(.)/82 %>% sqrt()
R_squared <- summary_pred_add %>% mutate(ave=true-mean(true)) %>% mutate(A=residual_bias^2, B=ave^2) %>% select(A,B) %>% apply(2, sum) %>% matrix(.,1,2) %>% data.frame() %>% mutate(TEST =1-X1/X2) %>% .[1,3]
print(paste("RMSE:", format(RMSE, digits=2), "R^2:", format(R_squared, digits=2)))
# AIC calculation
residuals <- summary_pred_add$residual_bias
sigma2 <- var(summary_pred_add$residual_bias)  # 残差の分散
log_likelihood <- -0.5 * length(residuals) * log(2 * pi * sigma2) - sum(residuals^2) / (2 * sigma2)
aic <- -2 * log_likelihood + 2 * ncol(x)
print(aic)

# residual plot comparison with Fukuda + (2020)
p <- summary_pred_add %>% 
  ggplot(aes(x=Fo, y=residual_bias))+
  geom_point(aes(col=2*sd), size=2, alpha=0.5)+
  theme_bw()+
  scale_colour_gradient(low="blue", high="red")+ 
  geom_hline(aes(yintercept = 0))+
  geom_hline(aes(yintercept = 0.3), linetype="dashed")+
  geom_hline(aes(yintercept = -0.3), linetype="dashed")+
  scale_x_continuous(limits=c(55,105))+
  scale_y_continuous(limits=c(-1.5,1))+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10))+
  theme(legend.position = "top") +
  guides(color = guide_colorbar(direction = "horizontal", title.position = "top"))
ggsave("residual plot.pdf", p, w=80, h=70, units="mm")

p<- summary_pred_add %>% 
  ggplot(aes(x=Fo,y=pred))+
  geom_errorbar(aes(ymin=pred-2*sd, ymax=pred+2*sd)) +
  geom_point(aes(x=Fo, y=true), col="yellow") + 
  geom_point(aes(), size=2, col="red")+
  theme_bw()+
  geom_hline(aes(yintercept = 0))+
  scale_x_continuous(limits=c(55,105))+
  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))
ggplotly(p)

p2<-summary_pred_add %>% 
  data.frame(., df_O2) %>% 
  ggplot(aes(x=R_Mg_Si/Fo*100,y=pred, label=Label))+
  geom_errorbar(aes(ymin=pred-2*sd, ymax=pred+2*sd)) +
  geom_point(aes(x=R_Mg_Si/Fo*100, y=true), col="yellow") + 
  geom_point(aes(), size=2, col="red")+
  theme_bw()+
  geom_hline(aes(yintercept = 0))+
  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))
ggplotly(p2)




