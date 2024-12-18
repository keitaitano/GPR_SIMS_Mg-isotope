# title: Correction of IMF during SIMS Mg isotope analysis
# author: Keita Itano 
# date:2024/08/13

# Library
library(tidyverse)
library(kernlab)
library(plotly)

# data import
df <- read.csv("Fukuda2021_data.csv") %>% mutate(session = as.factor(session))
summary_result <- data.frame()

#####  Session: Wild ####
df_train <- df %>% filter(type == "standard") %>% filter(session == 1)
df_test <- df %>% filter(type == "unknown")%>% filter(session == 1)

# pre-process 
x_train<- df_train %>% mutate(R_FeMg=(FeO/MgO), R_CaMg=(CaO/MgO), R_CrMg=(Cr2O3/MgO),  R_MnMg=(MnO/MgO), signal=Yield) %>% select(R_FeMg, R_CaMg, R_CrMg,R_MnMg)
x_test <- df_test %>% mutate(R_FeMg=(FeO/MgO), R_CaMg=(CaO/MgO), R_CrMg=(Cr2O3/MgO),  R_MnMg=(MnO/MgO), signal=Yield) %>% select(R_FeMg, R_CaMg,R_CrMg,R_MnMg)
measurement_errors <- df_train$bias_error
measurement_variances <- measurement_errors^2
noise <- measurement_variances %>% mean(.)/sd(df_train$bias)

# GPR fitting
fit_rbf <- gausspr(x_train, df_train$bias,
                   variance.model=T,
                   var=noise,
                   scaled=T,
                   kernel="rbfdot",
                   kpar=list(sigma=0.1))

summary_pred_add <- 
  data.frame(true=df_train$bias,
             pred= predict(fit_rbf, x_train),
             sd =  predict(fit_rbf, x_train, type="sdeviation")) %>% 
  mutate(residual_bias = true-pred,Fo=df_train$Fo)


# prediction performance
RMSE <- summary_pred_add$residual_bias %>% .^2 %>% sum(.)/82 %>% sqrt()
R_squared <- summary_pred_add %>% mutate(ave=true-mean(true)) %>% mutate(A=residual_bias^2, B=ave^2) %>% select(A,B) %>% apply(2, sum) %>% matrix(.,1,2) %>% data.frame() %>% mutate(TEST =1-X1/X2) %>% .[1,3]
print(paste("RMSE:", format(RMSE, digits=2), "R^2:", format(R_squared, digits=2)))

# residual plot
summary_pred_add %>% 
  ggplot(aes(x=Fo, y=residual_bias))+
  geom_point(aes(col=sd), size=2, alpha=0.5)+
  theme_bw()+
  scale_colour_gradient(low="blue", high="red")+ 
  geom_hline(aes(yintercept = 0))+
  geom_hline(aes(yintercept = 0.3), linetype="dashed")+
  geom_hline(aes(yintercept = -0.3), linetype="dashed")+
  scale_x_continuous(limits=c(0.55,1.05))+
  scale_y_continuous(limits=c(-1.5,1))

bias_pred <- predict(fit_rbf, x_test)
bias_pred_sd =  predict(fit_rbf, x_test, type="sdeviation")
df_pred_test <- data.frame(x_test, bias=bias_pred, SD_2=2*bias_pred_sd )
summary_result <- rbind(summary_result, data.frame(df_test$ID, df_pred_test))

#### Session: Kaba ####
df_train <- df %>% filter(type == "standard") %>% filter(session == 3)
df_test <- df %>% filter(type == "unknown")%>% filter(session == 2)

# pre-process 
x_train<- df_train %>% mutate(R_FeMg=(FeO/MgO), R_CaMg=(CaO/MgO), R_CrMg=(Cr2O3/MgO),  R_MnMg=(MnO/MgO), signal=Yield) %>% select(R_FeMg, R_CaMg, R_CrMg,R_MnMg)
x_test <- df_test %>% mutate(R_FeMg=(FeO/MgO), R_CaMg=(CaO/MgO), R_CrMg=(Cr2O3/MgO),  R_MnMg=(MnO/MgO), signal=Yield) %>% select(R_FeMg, R_CaMg,R_CrMg,R_MnMg)
measurement_errors <- df_train$bias_error
measurement_variances <- measurement_errors^2
noise <- measurement_variances %>% mean(.)/sd(df_train$bias)

# GPR fitting
fit_rbf <- gausspr(x_train, df_train$bias,
                   variance.model=T,
                   var=noise,
                   scaled=T,
                   kernel="rbfdot",
                   kpar=list(sigma=0.1))

summary_pred_add <- 
  data.frame(true=df_train$bias,
             pred= predict(fit_rbf, x_train),
             sd =  predict(fit_rbf, x_train, type="sdeviation")) %>% 
  mutate(residual_bias = true-pred,Fo=df_train$Fo)


# prediction performance
RMSE <- summary_pred_add$residual_bias %>% .^2 %>% sum(.)/82 %>% sqrt()
R_squared <- summary_pred_add %>% mutate(ave=true-mean(true)) %>% mutate(A=residual_bias^2, B=ave^2) %>% select(A,B) %>% apply(2, sum) %>% matrix(.,1,2) %>% data.frame() %>% mutate(TEST =1-X1/X2) %>% .[1,3]
print(paste("RMSE:", format(RMSE, digits=2), "R^2:", format(R_squared, digits=2)))

# residual plot
summary_pred_add %>% 
  ggplot(aes(x=Fo, y=residual_bias))+
  geom_point(aes(col=sd), size=2, alpha=0.5)+
  theme_bw()+
  scale_colour_gradient(low="blue", high="red")+ 
  geom_hline(aes(yintercept = 0))+
  geom_hline(aes(yintercept = 0.3), linetype="dashed")+
  geom_hline(aes(yintercept = -0.3), linetype="dashed")+
  scale_x_continuous(limits=c(0.55,1.05))+
  scale_y_continuous(limits=c(-1.5,1))

bias_pred <- predict(fit_rbf, x_test)
bias_pred_sd =  predict(fit_rbf, x_test, type="sdeviation")
df_pred_test <- data.frame(x_test, bias=bias_pred, SD_2=2*bias_pred_sd )
summary_result <- rbind(summary_result, data.frame(df_test$ID, df_pred_test))


#### Session: DOM ####
df_train <- df %>% filter(type == "standard") %>% filter(session == 3)
df_test <- df %>% filter(type == "unknown")%>% filter(session == 3)

# pre-process 
x_train<- df_train %>% mutate(R_FeMg=(FeO/MgO), R_CaMg=(CaO/MgO), R_CrMg=(Cr2O3/MgO),  R_MnMg=(MnO/MgO), signal=Yield) %>% select(R_FeMg, R_CaMg, R_CrMg,R_MnMg)
x_test <- df_test %>% mutate(R_FeMg=(FeO/MgO), R_CaMg=(CaO/MgO), R_CrMg=(Cr2O3/MgO),  R_MnMg=(MnO/MgO), signal=Yield) %>% select(R_FeMg, R_CaMg,R_CrMg,R_MnMg)
measurement_errors <- df_train$bias_error
measurement_variances <- measurement_errors^2
noise <- measurement_variances %>% mean(.)/sd(df_train$bias)

# GPR fitting
fit_rbf <- gausspr(x_train, df_train$bias,
                   variance.model=T,
                   var=noise,
                   scaled=T,
                   kernel="rbfdot",
                   kpar=list(sigma=0.1))

summary_pred_add <- 
  data.frame(true=df_train$bias,
             pred= predict(fit_rbf, x_train),
             sd =  predict(fit_rbf, x_train, type="sdeviation")) %>% 
  mutate(residual_bias = true-pred,Fo=df_train$Fo)
print(summary_pred_add)

# prediction performance
RMSE <- summary_pred_add$residual_bias %>% .^2 %>% sum(.)/82 %>% sqrt()
R_squared <- summary_pred_add %>% mutate(ave=true-mean(true)) %>% mutate(A=residual_bias^2, B=ave^2) %>% select(A,B) %>% apply(2, sum) %>% matrix(.,1,2) %>% data.frame() %>% mutate(TEST =1-X1/X2) %>% .[1,3]
print(paste("RMSE:", format(RMSE, digits=2), "R^2:", format(R_squared, digits=2)))

# residual plot
summary_pred_add %>% 
  ggplot(aes(x=Fo, y=residual_bias))+
  geom_point(aes(col=sd), size=2, alpha=0.5)+
  theme_bw()+
  scale_colour_gradient(low="blue", high="red")+ 
  geom_hline(aes(yintercept = 0))+
  geom_hline(aes(yintercept = 0.3), linetype="dashed")+
  geom_hline(aes(yintercept = -0.3), linetype="dashed")+
  scale_x_continuous(limits=c(0.55,1.05))+
  scale_y_continuous(limits=c(-1.5,1))

bias_pred <- predict(fit_rbf, x_test)
bias_pred_sd =  predict(fit_rbf, x_test, type="sdeviation")
df_pred_test <- data.frame(x_test, bias=bias_pred, SD_2=2*bias_pred_sd )
summary_result <- rbind(summary_result, data.frame(df_test$ID, df_pred_test))

view(summary_result)


df_summary <- read.csv("comparison_Fukuda2021.csv")


p_tmp <- 
  df_summary %>% 
  filter(label !="T149F1" ) %>% filter(label!="T77F4") %>% filter(label!= "T77F7") %>% 
  ggplot(.)+
  geom_abline(slope=0, intercept =-0.15, col="black", linetype="dashed")+
  geom_errorbar(aes(x=CaO, ymin=delta25Mg_fukuda-error_2, ymax=delta25Mg_fukuda+error_2, col=label2), width=0)+
  geom_point(aes(x=CaO, y=delta25Mg_fukuda,  col=label2), size=3.5,shape=1) +
  geom_errorbar(aes(x=CaO, ymin=delta25Mg_itano-error_1, ymax=delta25Mg_itano+error_1), width=0.005)+
  geom_point(aes(x=CaO, y=delta25Mg_itano,  col=label2), size=3.5,shape=16) +
  geom_abline(slope = 0, intercept = 0)+
  coord_cartesian(xlim = c(0, 0.3), ylim = c(-5, 1.5), )+
  scale_color_manual(values = c("red", "blue", "green", "pink", "#FFC20E"))+
  scale_x_continuous(breaks = seq(0, 0.3, 0.1), ) +
  scale_y_continuous(breaks = seq(-5, 1, 1)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))+
  theme(legend.position = "top")
p_tmp

ggsave("Fukuda2021.pdf", p_tmp, w=80, h=80, units="mm")
