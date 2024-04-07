# Plotting grouped boxplot with multiple panels 
# Order it like this
#     x - the data 
#     y - ldaps, ann, rf, xgb
#     z - tmax, tmin 
#     t - R2, RMSE, BIAS 

# R2
ldaps_r2_tmax <- scan('./results/ann/ldaps_r2_score_tmax.txt')
ann_r2_tmax   <- scan('./results/ann/ann_r2_score_tmax.txt')
rf_r2_tmax    <- scan('./results/randomForest/rf_r2_score_tmax.txt')
xgb_r2_tmax   <- scan('./results/boost/xgb_r2_score_tmax.txt')

ldaps_r2_tmin <- scan('./results/ann/ldaps_r2_score_tmin.txt')
ann_r2_tmin   <- scan('./results/ann/ann_r2_score_tmin.txt')
rf_r2_tmin    <- scan('./results/randomForest/rf_r2_score_tmin.txt')
xgb_r2_tmin   <- scan('./results/boost/xgb_r2_score_tmin.txt')

# RMSE 
ldaps_rmse_tmax <- scan('./results/ann/ldaps_rmse_score_tmax.txt')
ann_rmse_tmax   <- scan('./results/ann/ann_rmse_score_tmax.txt')
rf_rmse_tmax    <- scan('./results/randomForest/rf_rmse_score_tmax.txt')
xgb_rmse_tmax   <- scan('./results/boost/xgb_rmse_score_tmax.txt')

ldaps_rmse_tmin <- scan('./results/ann/ldaps_rmse_score_tmin.txt')
ann_rmse_tmin   <- scan('./results/ann/ann_rmse_score_tmin.txt')
rf_rmse_tmin    <- scan('./results/randomForest/rf_rmse_score_tmin.txt')
xgb_rmse_tmin   <- scan('./results/boost/xgb_rmse_score_tmin.txt')

# BIAS
ldaps_bs_tmax <- scan('./results/ann/ldaps_bs_score_tmax.txt')
ann_bs_tmax   <- scan('./results/ann/ann_bs_score_tmax.txt')
rf_bs_tmax    <- scan('./results/randomForest/rf_bs_score_tmax.txt')
xgb_bs_tmax   <- scan('./results/boost/xgb_bs_score_tmax.txt')

ldaps_bs_tmin <- scan('./results/ann/ldaps_bs_score_tmin.txt')
ann_bs_tmin   <- scan('./results/ann/ann_bs_score_tmin.txt')
rf_bs_tmin    <- scan('./results/randomForest/rf_bs_score_tmin.txt')
xgb_bs_tmin   <- scan('./results/boost/xgb_bs_score_tmin.txt')

# Setting up the dataframe 
x <- c( c(ldaps_r2_tmax,ann_r2_tmax,rf_r2_tmax,xgb_r2_tmax), 
        c(ldaps_r2_tmin,ann_r2_tmin,rf_r2_tmin,xgb_r2_tmin),
        c(ldaps_rmse_tmax,ann_rmse_tmax,rf_rmse_tmax,xgb_rmse_tmin),
        c(ldaps_rmse_tmin,ann_rmse_tmin,rf_rmse_tmin,xgb_rmse_tmin),
        c(ldaps_bs_tmax,ann_bs_tmax,rf_bs_tmax,xgb_bs_tmax),
        c(ldaps_bs_tmin,ann_bs_tmin,rf_bs_tmin,xgb_bs_tmin) ) 

y <- rep(rep(c('ldaps','ann','rf','xgb'),each=10),6)

z <- rep(rep(c('tmax','tmin'),each=40),3) 

t <- rep(c('R2','RMSE','BIAS'),each=80)

DF <- data.frame(x,y,z,t)

# Plotting the  grouped boxplot with multiple panels 
# Order it like this
#     x - the data 
#     y - ldaps, ann, rf, xgb
#     z - tmax, tmin 
#     t - R2, RMSE, BIAS 

library(ggplot2)

png('BoxPlots.png')

ggplot(DF,aes(x=factor(y),
              y=x,
              color=factor(z))) + geom_boxplot() + facet_wrap(~t,ncol=3) + 
  labs(color='Target')

dev.off()

# 
# Calculating skill scores 
#   R2   - coefficient of determination
#   Bs   - bias 
#   RMSE - root mean square error 
# 
R2 <- function( yhat, y ){
  ybar <- mean(y,na.rm=TRUE)
  SSres <- sum((yhat-y)**2)
  SStot <- sum((y-ybar)**2)
  return(1-SSres/SStot)
}

# Bs
#
BS <- function( yhat, y ){
  return( mean(yhat-y,na.rm=TRUE) )
}

# RMSE
#
RMSE <- function( yhat, y ){
  return( sqrt( mean( (yhat-y)**2, na.rm = TRUE ) ) )
}

# Next, plotting Scatterplots for tmax 
#   1. ldaps forecast vs. observed
#   2. ann forecast vs. observed
#   3. rf forecast vs. observed
#   4. xgb forecast vs. observed
#
# LDAPS
forecast <- scan('./results/ann/uncorrected_forecast_tmax.txt')
observed <- scan('./results/ann/observed_tmax.txt')
r2 <- format(R2(forecast,observed),digits=2)
bs <- format(BS(forecast,observed),digits=2)
rmse <- format(RMSE(forecast,observed),digits=2)
txt <- paste0('R2 = ',r2,'\nBias = ',bs,'\nRMSE = ',rmse)
png('ScatterPlot_LDAPS-Tmax.png')
smoothScatter(observed,forecast,main='LDAPS Forecast vs. Observed Tmax')
abline(coef=c(0,1))
mtext(txt,side=1,line=-2,cex=1.5)
dev.off()

forecast <- scan('./results/ann/uncorrected_forecast_tmin.txt')
observed <- scan('./results/ann/observed_tmin.txt')
r2 <- format(R2(forecast,observed),digits=2)
bs <- format(BS(forecast,observed),digits=2)
rmse <- format(RMSE(forecast,observed),digits=2)
txt <- paste0('R2 = ',r2,'\nBias = ',bs,'\nRMSE = ',rmse)
png('ScatterPlot_LDAPS-Tmin.png')
smoothScatter(observed,forecast,main='LDAPS Forecast vs. Observed Tmin')
abline(coef=c(0,1))
mtext(txt,side=1,line=-2,cex=1.5)
dev.off()


# ANN
forecast <- scan('./results/ann/ann_prediction_tmax.txt')
observed <- scan('./results/ann/observed_tmax.txt')
r2 <- format(R2(forecast,observed),digits=2)
bs <- format(BS(forecast,observed),digits=2)
rmse <- format(RMSE(forecast,observed),digits=2)
txt <- paste0('R2 = ',r2,'\nBias = ',bs,'\nRMSE = ',rmse)
png('ScatterPlot_ANN-Tmax.png')
smoothScatter(observed,forecast,main='ANN Forecast vs. Observed Tmax')
abline(coef=c(0,1))
mtext(txt,side=1,line=-2,cex=1.5)
dev.off()

forecast <- scan('./results/ann/ann_prediction_tmin.txt')
observed <- scan('./results/ann/observed_tmin.txt')
r2 <- format(R2(forecast,observed),digits=2)
bs <- format(BS(forecast,observed),digits=2)
rmse <- format(RMSE(forecast,observed),digits=2)
txt <- paste0('R2 = ',r2,'\nBias = ',bs,'\nRMSE = ',rmse)
png('ScatterPlot_ANN-Tmin.png')
smoothScatter(observed,forecast,main='ANN Forecast vs. Observed Tmin')
abline(coef=c(0,1))
mtext(txt,side=1,line=-2,cex=1.5)
dev.off()


# RF
forecast <- scan('./results/randomForest/rf_prediction_tmax.txt')
observed <- scan('./results/randomForest/observed_tmax.txt')
r2 <- format(R2(forecast,observed),digits=2)
bs <- format(BS(forecast,observed),digits=2)
rmse <- format(RMSE(forecast,observed),digits=2)
txt <- paste0('R2 = ',r2,'\nBias = ',bs,'\nRMSE = ',rmse)
png('ScatterPlot_RF-Tmax.png')
smoothScatter(observed,forecast,main='RF Forecast vs. Observed Tmax')
abline(coef=c(0,1))
mtext(txt,side=1,line=-2,cex=1.5)
dev.off()

forecast <- scan('./results/randomForest/rf_prediction_tmin.txt')
observed <- scan('./results/randomForest/observed_tmin.txt')
r2 <- format(R2(forecast,observed),digits=2)
bs <- format(BS(forecast,observed),digits=2)
rmse <- format(RMSE(forecast,observed),digits=2)
txt <- paste0('R2 = ',r2,'\nBias = ',bs,'\nRMSE = ',rmse)
png('ScatterPlot_RF-Tmin.png')
smoothScatter(observed,forecast,main='RF Forecast vs. Observed Tmin')
abline(coef=c(0,1))
mtext(txt,side=1,line=-2,cex=1.5)
dev.off()




# XGB
forecast <- scan('./results/boost/xgb_prediction_tmax.txt')
observed <- scan('./results/boost/observed_tmax.txt')
r2 <- format(R2(forecast,observed),digits=2)
bs <- format(BS(forecast,observed),digits=2)
rmse <- format(RMSE(forecast,observed),digits=2)
txt <- paste0('R2 = ',r2,'\nBias = ',bs,'\nRMSE = ',rmse)
png('ScatterPlot_XGB-Tmax.png')
smoothScatter(observed,forecast,main='XGBoost Forecast vs. Observed Tmax')
abline(coef=c(0,1))
mtext(txt,side=1,line=-2,cex=1.5)
dev.off()

forecast <- scan('./results/boost/xgb_prediction_tmin.txt')
observed <- scan('./results/boost/observed_tmin.txt')
r2 <- format(R2(forecast,observed),digits=2)
bs <- format(BS(forecast,observed),digits=2)
rmse <- format(RMSE(forecast,observed),digits=2)
txt <- paste0('R2 = ',r2,'\nBias = ',bs,'\nRMSE = ',rmse)
png('ScatterPlot_XGB-Tmin.png')
smoothScatter(observed,forecast,main='XGBoost Forecast vs. Observed Tmin')
abline(coef=c(0,1))
mtext(txt,side=1,line=-2,cex=1.5)
dev.off()


