# Loading data

library(tidyverse)

f <- 'Bias_correction_ucl.csv'

df <- read_csv(f)

# Remove rows where station and Date columns are missing 

df <- df[ -which(is.na(df[,'station']) & is.na(df[,'Date'])), ]

# Plotting box plots to detect unusual values

output_dir <- 'BoxPlots'

if (!dir.exists(output_dir)){
  
  dir.create(output_dir)

}else{

  print(paste(output_dir,'already exists'))

}


Draw_BoxPlot <- function(df,var,title,ylab){
  png(paste0('./BoxPlots/',var,'.png'))
  boxplot(df[var],main=title,ylab=ylab)
  dev.off()
}

Draw_BoxPlot(df,'Present_Tmax','Maximum Daily Temperature', '(degrees Celcius)')

Draw_BoxPlot(df,'Present_Tmin','Minimum Daily Temperature', '(degrees Celcius)')

Draw_BoxPlot(df,'LDAPS_RHmin','Minimum Daily Relative Humidity', '(%)')

Draw_BoxPlot(df,'LDAPS_RHmax','Maximum Daily Relative Humidity', '(%)')

Draw_BoxPlot(df,'LDAPS_Tmax_lapse','Next Day Maximum Daily Temperature', '(%)')

Draw_BoxPlot(df,'LDAPS_Tmin_lapse','Next Day Minimum Daily Temperature', '(%)')

Draw_BoxPlot(df,'LDAPS_WS','Next Day Average Wind Speed', '(m/s)')

Draw_BoxPlot(df,'LDAPS_LH','Next Day Average Latent Heat Flux', '(W/m2)')

Draw_BoxPlot(df,'LDAPS_CC1','Next Day Average Cloud Cover from 0-5h',
             '(Proportion)')

Draw_BoxPlot(df,'LDAPS_CC2','Next Day Average Cloud Cover from 6-11h',
             '(Proportion)')

Draw_BoxPlot(df,'LDAPS_CC3','Next Day Average Cloud Cover from 12-17h',
             '(Proportion)')

Draw_BoxPlot(df,'LDAPS_CC4','Next Day Average Cloud Cover from 18-23h',
             '(Proportion)')

Draw_BoxPlot(df,'LDAPS_PPT1','Next Day Average Precipitation from 0-5h',
             '(mm)')

Draw_BoxPlot(df,'LDAPS_PPT2','Next Day Average Precipitation from 6-11h',
             '(mm)')

Draw_BoxPlot(df,'LDAPS_PPT3','Next Day Average Precipitation from 12-17h',
             '(mm)')

Draw_BoxPlot(df,'LDAPS_PPT4','Next Day Average Precipitation from 18-23h',
             '(mm)')

Draw_BoxPlot(df,'lat','Latitude','(degrees North)')

Draw_BoxPlot(df,'lon','Longitude','(degrees East)')

Draw_BoxPlot(df,'DEM','Elevation','(m)')

Draw_BoxPlot(df,'Slope','Slope','(Degrees)')

Draw_BoxPlot(df,'Solar radiation','Solar Radiation','(Wh/m2)')

Draw_BoxPlot(df,'Next_Tmax','Next-day maximum air temperature',
             '(degrees Celcius)')

Draw_BoxPlot(df,'Next_Tmin','Next-day minimum air temperature',
             '(degrees Celcius)')

# No unusual values detected 

# Counting number of missing values of each column 

percentage_missing <- function(df,column_number){
  return( sum( is.na(df[,column_number]) ) / nrow( df[,column_number] ) * 100 )
}

for (i in 1: ncol(df)){
  pmiss <- percentage_missing(df,i)
  feat <- colnames(df[i])
  print(paste0(format(feat,width=10),' missing = ',format(pmiss,digits=2),'%'))
}

# At most 0.97% missing values

# Imputation done using R MICE package (cart) and R MissForest package
#   then choose one with lowest KS D-statistic
#
library(mice)

library(missForest)

# changing Solar radiation to Solar_radiation column name
#   mice imputation does not work on column names with space 

colnames(df)[which(names(df) == 'Solar radiation')] <- 'Solar_radiation'

# select only numeric variables for imputation
#   which does not include station because it is just a label 
#   and does not include date because it is 'yyyy-mm-dd' string
#     also, the time information is already contained in solar radiation 
#       variable 
#
df_numeric <- select(df,-c('station','Date'))

# Imputing missing values by mice (cart)
# This may take some time, better save result to a file
#
if (file.exists('imputed_mice_cart.data')){
  load(file='imputed_mice_cart.data')
}else{
  df_imputed_mice_cart <- complete(mice(df_numeric,method='cart'))
  save(df_imputed_mice_cart,file='imputed_mice_cart.data')
}

# Imputing missing values by missing forest
# this may take 30 minutes, better to save the results to a file
#
if (file.exists('imputed_miss_forest.data')){
  load(file='imputed_miss_forest.data')
}else{
  df_imputed_miss_forest <- missForest(data.frame(df_numeric),verbose=TRUE)$ximp 
  save(df_imputed_miss_forest,file='imputed_miss_forest.data')
}

# Plotting histograms and calculating KS-Statistic
# 

if (!dir.exists('./histogram_original_vs_imputed_distribution/')){
  
  dir.create('./histogram_original_vs_imputed_distribution/')

}

get_boxplot_ks <- function( original, imputed_mice_cart, 
                            imputed_miss_forest, varname){

  ks_mice_cart <- format(ks.test(original,
                                 imputed_mice_cart)$statistic,
                         digits=3)

  ks_miss_forest <- format(ks.test(original,
                                 imputed_miss_forest)$statistic,digits=3)

  png(paste0('./histogram_original_vs_imputed_distribution/',varname,'.png'))

  par(mfrow=(c(1,3)))

  hist(original,main=paste(varname,'\nOriginal Distribution'),cex.main=1.2,
     cex.lab=1.5)

  hist(imputed_mice_cart,main=paste(varname,'\nMICE-CART Imputed',
                                  '\nKS statistic = ',ks_mice_cart),
     cex.main=1.2,
     cex.lab=1.5)

  hist(imputed_miss_forest,main=paste(varname,
                                    '\nMISS-FOREST Imputed',
                                    '\nKS statistic = ',ks_miss_forest),
    cex.main=1.2,
    cex.lab=1.5)

  dev.off()

}

get_boxplot_ks( df$Present_Tmax, 
                df_imputed_mice_cart$Present_Tmax, 
                df_imputed_miss_forest$Present_Tmax, 
                'Present_Tmax' )

get_boxplot_ks( df$Present_Tmin, 
                df_imputed_mice_cart$Present_Tmin, 
                df_imputed_miss_forest$Present_Tmin, 
                'Present_Tmin' )

get_boxplot_ks( df$LDAPS_RHmin, 
                df_imputed_mice_cart$LDAPS_RHmin, 
                df_imputed_miss_forest$LDAPS_RHmin, 
                'LDAPS_RHmin' )

get_boxplot_ks( df$LDAPS_RHmax, 
                df_imputed_mice_cart$LDAPS_RHmax, 
                df_imputed_miss_forest$LDAPS_RHmax, 
                'LDAPS_RHmax' )

get_boxplot_ks( df$LDAPS_Tmax_lapse, 
                df_imputed_mice_cart$LDAPS_Tmax_lapse, 
                df_imputed_miss_forest$LDAPS_Tmax_lapse, 
                'LDAPS_Tmax_lapse' )

get_boxplot_ks( df$LDAPS_Tmin_lapse, 
                df_imputed_mice_cart$LDAPS_Tmin_lapse, 
                df_imputed_miss_forest$LDAPS_Tmin_lapse, 
                'LDAPS_Tmin_lapse' )

get_boxplot_ks( df$LDAPS_WS, 
                df_imputed_mice_cart$LDAPS_WS, 
                df_imputed_miss_forest$LDAPS_WS, 
                'LDAPS_WS' )

get_boxplot_ks( df$LDAPS_LH, 
                df_imputed_mice_cart$LDAPS_LH, 
                df_imputed_miss_forest$LDAPS_LH, 
                'LDAPS_LH' )

get_boxplot_ks( df$LDAPS_CC1, 
                df_imputed_mice_cart$LDAPS_CC1, 
                df_imputed_miss_forest$LDAPS_CC1, 
                'LDAPS_CC1' )

get_boxplot_ks( df$LDAPS_CC2, 
                df_imputed_mice_cart$LDAPS_CC2, 
                df_imputed_miss_forest$LDAPS_CC2, 
                'LDAPS_CC2' )

get_boxplot_ks( df$LDAPS_CC3, 
                df_imputed_mice_cart$LDAPS_CC3, 
                df_imputed_miss_forest$LDAPS_CC3, 
                'LDAPS_CC3' )

get_boxplot_ks( df$LDAPS_CC4, 
                df_imputed_mice_cart$LDAPS_CC4, 
                df_imputed_miss_forest$LDAPS_CC4, 
                'LDAPS_CC4' )

get_boxplot_ks( df$LDAPS_PPT1, 
                df_imputed_mice_cart$LDAPS_PPT1, 
                df_imputed_miss_forest$LDAPS_PPT1, 
                'LDAPS_PPT1' )

get_boxplot_ks( df$LDAPS_PPT2, 
                df_imputed_mice_cart$LDAPS_PPT2, 
                df_imputed_miss_forest$LDAPS_PPT2, 
                'LDAPS_PPT2' )

get_boxplot_ks( df$LDAPS_PPT3, 
                df_imputed_mice_cart$LDAPS_PPT3, 
                df_imputed_miss_forest$LDAPS_PPT3, 
                'LDAPS_PPT3' )

get_boxplot_ks( df$LDAPS_PPT4, 
                df_imputed_mice_cart$LDAPS_PPT4, 
                df_imputed_miss_forest$LDAPS_PPT4, 
                'LDAPS_PPT4' )

# MICE cart is chosen because it has lower KS-statistic 
#   most of the time 
# 
for (i in 3:ncol(df)){
  df[,i] <- df_imputed_mice_cart[,i-2]
}

# Shuffle 
#
df <- df[sample(1:nrow(df)),]

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

# 10-fold cross-validation with MinMaxScaler 
#   for artificial neural network (ANN)
#

library(neuralnet)

fold <- 10

N <- nrow(df) / fold

index <- seq(from = 1, to = nrow(df)+N, by = N )

prediction_Tmax <- c()

prediction_Tmin <- c()

observed_Tmax <- c()

observed_Tmin <- c()

uncorrected_Tmax <- c()
  
uncorrected_Tmin <- c()

r2_score_tmax <- c()

r2_score_tmin <- c()

bs_score_tmax <- c()

bs_score_tmin <- c()

rmse_score_tmax <- c()

rmse_score_tmin <- c()

ldaps_r2_score_tmax <- c()

ldaps_r2_score_tmin <- c()

ldaps_bs_score_tmax <- c()

ldaps_bs_score_tmin <- c()

ldaps_rmse_score_tmax <- c()

ldaps_rmse_score_tmin <- c()

for (i in 1:(length(index)-1)){
  
  begin <- index[i]
  
  end <- index[i+1]-1

    print(paste('Evaluating test dataset index = ',
              i,
              begin,
              end))

  test_data  <- df[begin:end,]  
  
  train_data <- df[-c(begin:end),]
  
  # scaling 
  
  train_data_scaled <- train_data
  
  test_data_scaled <- test_data
  
  train_data_scaled$Present_Tmax <- 
    ( train_data$Present_Tmax - min(train_data$Present_Tmax) ) / 
    ( max(train_data$Present_Tmax) - min(train_data$Present_Tmax) )
  
  train_data_scaled$Present_Tmin <- 
    ( train_data$Present_Tmin - min(train_data$Present_Tmin) ) / 
    ( max(train_data$Present_Tmin) - min(train_data$Present_Tmin) )
  
  train_data_scaled$LDAPS_RHmin <- 
    ( train_data$LDAPS_RHmin - min(train_data$LDAPS_RHmin) ) / 
    ( max(train_data$LDAPS_RHmin) - min(train_data$LDAPS_RHmin) )

  train_data_scaled$LDAPS_RHmax <- 
    ( train_data$LDAPS_RHmax - min(train_data$LDAPS_RHmax) ) / 
    ( max(train_data$LDAPS_RHmax) - min(train_data$LDAPS_RHmax) )
  
  train_data_scaled$LDAPS_Tmax_lapse <- 
    ( train_data$LDAPS_Tmax_lapse - min(train_data$LDAPS_Tmax_lapse) ) / 
    ( max(train_data$LDAPS_Tmax_lapse) - min(train_data$LDAPS_Tmax_lapse) )

  train_data_scaled$LDAPS_Tmin_lapse <- 
    ( train_data$LDAPS_Tmin_lapse - min(train_data$LDAPS_Tmin_lapse) ) / 
    ( max(train_data$LDAPS_Tmin_lapse) - min(train_data$LDAPS_Tmin_lapse) )
  
  train_data_scaled$LDAPS_WS <- 
    ( train_data$LDAPS_WS - min(train_data$LDAPS_WS) ) / 
    ( max(train_data$LDAPS_WS) - min(train_data$LDAPS_WS) )
  
  train_data_scaled$LDAPS_LH <- 
    ( train_data$LDAPS_LH - min(train_data$LDAPS_LH) ) / 
    ( max(train_data$LDAPS_LH) - min(train_data$LDAPS_LH) )
  
  train_data_scaled$LDAPS_CC1 <- 
    ( train_data$LDAPS_CC1 - min(train_data$LDAPS_CC1) ) / 
    ( max(train_data$LDAPS_CC1) - min(train_data$LDAPS_CC1) )
  
  train_data_scaled$LDAPS_CC2 <- 
    ( train_data$LDAPS_CC2 - min(train_data$LDAPS_CC2) ) / 
    ( max(train_data$LDAPS_CC2) - min(train_data$LDAPS_CC2) )
  
  train_data_scaled$LDAPS_CC3 <- 
    ( train_data$LDAPS_CC3 - min(train_data$LDAPS_CC3) ) / 
    ( max(train_data$LDAPS_CC3) - min(train_data$LDAPS_CC3) )

  train_data_scaled$LDAPS_CC4 <- 
    ( train_data$LDAPS_CC4 - min(train_data$LDAPS_CC4) ) / 
    ( max(train_data$LDAPS_CC4) - min(train_data$LDAPS_CC4) )

  train_data_scaled$LDAPS_PPT1 <- 
    ( train_data$LDAPS_PPT1 - min(train_data$LDAPS_PPT1) ) / 
    ( max(train_data$LDAPS_PPT1) - min(train_data$LDAPS_PPT1) )

  train_data_scaled$LDAPS_PPT2 <- 
    ( train_data$LDAPS_PPT2 - min(train_data$LDAPS_PPT2) ) / 
    ( max(train_data$LDAPS_PPT2) - min(train_data$LDAPS_PPT2) )

  train_data_scaled$LDAPS_PPT3 <- 
    ( train_data$LDAPS_PPT3 - min(train_data$LDAPS_PPT3) ) / 
    ( max(train_data$LDAPS_PPT3) - min(train_data$LDAPS_PPT3) )
  
  train_data_scaled$LDAPS_PPT4 <- 
    ( train_data$LDAPS_PPT4 - min(train_data$LDAPS_PPT4) ) / 
    ( max(train_data$LDAPS_PPT4) - min(train_data$LDAPS_PPT4) )
  
  train_data_scaled$lat <- 
    ( train_data$lat - min(train_data$lat) ) / 
    ( max(train_data$lat) - min(train_data$lat) )
  
  train_data_scaled$lon <- 
    ( train_data$lon - min(train_data$lon) ) / 
    ( max(train_data$lon) - min(train_data$lon) )
  
  train_data_scaled$DEM <- 
    ( train_data$DEM - min(train_data$DEM) ) / 
    ( max(train_data$DEM) - min(train_data$DEM) )
  
  train_data_scaled$Slope <- 
    ( train_data$Slope - min(train_data$Slope) ) / 
    ( max(train_data$Slope) - min(train_data$Slope) )
  
  train_data_scaled$Solar_radiation <- 
    ( train_data$Solar_radiation - min(train_data$Solar_radiation) ) / 
    ( max(train_data$Solar_radiation) - min(train_data$Solar_radiation) )

  train_data_scaled$Next_Tmax <- 
    ( train_data$Next_Tmax - min(train_data$Next_Tmax) ) / 
    ( max(train_data$Next_Tmax) - min(train_data$Next_Tmax) )

  train_data_scaled$Next_Tmin <- 
    ( train_data$Next_Tmin - min(train_data$Next_Tmin) ) / 
    ( max(train_data$Next_Tmin) - min(train_data$Next_Tmin) )
  
  
  test_data_scaled$Present_Tmax <- 
    ( test_data$Present_Tmax - min(train_data$Present_Tmax) ) / 
    ( max(train_data$Present_Tmax) - min(train_data$Present_Tmax) )
  
  test_data_scaled$Present_Tmin <- 
    ( test_data$Present_Tmin - min(train_data$Present_Tmin) ) / 
    ( max(train_data$Present_Tmin) - min(train_data$Present_Tmin) )
  
  test_data_scaled$LDAPS_RHmin <- 
    ( test_data$LDAPS_RHmin - min(train_data$LDAPS_RHmin) ) / 
    ( max(train_data$LDAPS_RHmin) - min(train_data$LDAPS_RHmin) )
  
  test_data_scaled$LDAPS_RHmax <- 
    ( test_data$LDAPS_RHmax - min(train_data$LDAPS_RHmax) ) / 
    ( max(train_data$LDAPS_RHmax) - min(train_data$LDAPS_RHmax) )
  
  test_data_scaled$LDAPS_Tmax_lapse <- 
    ( test_data$LDAPS_Tmax_lapse - min(train_data$LDAPS_Tmax_lapse) ) / 
    ( max(train_data$LDAPS_Tmax_lapse) - min(train_data$LDAPS_Tmax_lapse) )
  
  test_data_scaled$LDAPS_Tmin_lapse <- 
    ( test_data$LDAPS_Tmin_lapse - min(train_data$LDAPS_Tmin_lapse) ) / 
    ( max(train_data$LDAPS_Tmin_lapse) - min(train_data$LDAPS_Tmin_lapse) )
  
  test_data_scaled$LDAPS_WS <- 
    ( test_data$LDAPS_WS - min(train_data$LDAPS_WS) ) / 
    ( max(train_data$LDAPS_WS) - min(train_data$LDAPS_WS) )
  
  test_data_scaled$LDAPS_LH <- 
    ( test_data$LDAPS_LH - min(train_data$LDAPS_LH) ) / 
    ( max(train_data$LDAPS_LH) - min(train_data$LDAPS_LH) )
  
  test_data_scaled$LDAPS_CC1 <- 
    ( test_data$LDAPS_CC1 - min(train_data$LDAPS_CC1) ) / 
    ( max(train_data$LDAPS_CC1) - min(train_data$LDAPS_CC1) )
  
  test_data_scaled$LDAPS_CC2 <- 
    ( test_data$LDAPS_CC2 - min(train_data$LDAPS_CC2) ) / 
    ( max(train_data$LDAPS_CC2) - min(train_data$LDAPS_CC2) )
  
  test_data_scaled$LDAPS_CC3 <- 
    ( test_data$LDAPS_CC3 - min(train_data$LDAPS_CC3) ) / 
    ( max(train_data$LDAPS_CC3) - min(train_data$LDAPS_CC3) )
  
  test_data_scaled$LDAPS_CC4 <- 
    ( test_data$LDAPS_CC4 - min(train_data$LDAPS_CC4) ) / 
    ( max(train_data$LDAPS_CC4) - min(train_data$LDAPS_CC4) )
  
  test_data_scaled$LDAPS_PPT1 <- 
    ( test_data$LDAPS_PPT1 - min(train_data$LDAPS_PPT1) ) / 
    ( max(train_data$LDAPS_PPT1) - min(train_data$LDAPS_PPT1) )
  
  test_data_scaled$LDAPS_PPT2 <- 
    ( test_data$LDAPS_PPT2 - min(train_data$LDAPS_PPT2) ) / 
    ( max(train_data$LDAPS_PPT2) - min(train_data$LDAPS_PPT2) )
  
  test_data_scaled$LDAPS_PPT3 <- 
    ( test_data$LDAPS_PPT3 - min(train_data$LDAPS_PPT3) ) / 
    ( max(train_data$LDAPS_PPT3) - min(train_data$LDAPS_PPT3) )
  
  test_data_scaled$LDAPS_PPT4 <- 
    ( test_data$LDAPS_PPT4 - min(train_data$LDAPS_PPT4) ) / 
    ( max(train_data$LDAPS_PPT4) - min(train_data$LDAPS_PPT4) )
  
  test_data_scaled$lat <- 
    ( test_data$lat - min(train_data$lat) ) / 
    ( max(train_data$lat) - min(train_data$lat) )
  
  test_data_scaled$lon <- 
    ( test_data$lon - min(train_data$lon) ) / 
    ( max(train_data$lon) - min(train_data$lon) )
  
  test_data_scaled$DEM <- 
    ( test_data$DEM - min(train_data$DEM) ) / 
    ( max(train_data$DEM) - min(train_data$DEM) )
  
  test_data_scaled$Slope <- 
    ( test_data$Slope - min(train_data$Slope) ) / 
    ( max(train_data$Slope) - min(train_data$Slope) )
  
  test_data_scaled$Solar_radiation <- 
    ( test_data$Solar_radiation - min(train_data$Solar_radiation) ) / 
    ( max(train_data$Solar_radiation) - min(train_data$Solar_radiation) )
  
  test_data_scaled$Next_Tmax <- 
    ( test_data$Next_Tmax - min(train_data$Next_Tmax) ) / 
    ( max(train_data$Next_Tmax) - min(train_data$Next_Tmax) )

  test_data_scaled$Next_Tmin <- 
    ( test_data$Next_Tmin - min(train_data$Next_Tmin) ) / 
    ( max(train_data$Next_Tmin) - min(train_data$Next_Tmin) )
  
  # creating train and test datasets
  
  train_Tmax <- select(train_data_scaled,-c('Next_Tmin',
                                            'station',
                                            'Date'))
  
  train_Tmin <- select(train_data_scaled,-c('Next_Tmax',
                                            'station',
                                            'Date'))
  
  test_Tmax <- select(test_data_scaled,-c('Next_Tmin',
                                          'station',
                                          'Date',
                                          'Next_Tmax'))
  
  test_Tmin <- select(test_data_scaled,-c('Next_Tmax',
                                          'station',
                                          'Date',
                                          'Next_Tmin'))

  # Neural network 

  model_Tmax <- neuralnet(
    Next_Tmax ~ Present_Tmax + Present_Tmin + 
      LDAPS_RHmin + LDAPS_RHmax +
      LDAPS_Tmax_lapse + LDAPS_Tmin_lapse + 
      LDAPS_WS + LDAPS_LH + 
      LDAPS_CC1 + LDAPS_CC2 + LDAPS_CC3 + LDAPS_CC4 + 
      LDAPS_PPT1 + LDAPS_PPT2 + LDAPS_PPT3 + LDAPS_PPT4 + 
      lat + lon + DEM + Slope + Solar_radiation, 
    data = train_Tmax, 
    hidden = c(3,2),
    linear.output=FALSE,
    lifesign = 'full',
    threshold = 0.05
  )

  model_Tmin <- neuralnet(
    Next_Tmin ~ Present_Tmax + Present_Tmin + 
      LDAPS_RHmin + LDAPS_RHmax +
      LDAPS_Tmax_lapse + LDAPS_Tmin_lapse + 
      LDAPS_WS + LDAPS_LH + 
      LDAPS_CC1 + LDAPS_CC2 + LDAPS_CC3 + LDAPS_CC4 + 
      LDAPS_PPT1 + LDAPS_PPT2 + LDAPS_PPT3 + LDAPS_PPT4 + 
      lat + lon + DEM + Slope + Solar_radiation, 
    data = train_Tmin, 
    hidden = c(3,2),
    linear.output=FALSE,
    lifesign = 'full',
    threshold = 0.05
  )
  
  # Getting the results 
  
  pred_Tmax <- compute(model_Tmax,test_Tmax)$net.result
  
  pred_Tmin <- compute(model_Tmin,test_Tmin)$net.result
  
  # Unscaling the results for Tmax
  
  rng <- max(train_data$Next_Tmax) - min(train_data$Next_Tmax)
  
  mxm <- max(train_data$Next_Tmax)
  
  mnm <- min(train_data$Next_Tmax)
  
  pred_Tmax <- pred_Tmax * rng + mnm 
  
  # Unscaling the results for Tmin
  
  rng <- max(train_data$Next_Tmin) - min(train_data$Next_Tmin)
  
  mxm <- max(train_data$Next_Tmin)
  
  mnm <- min(train_data$Next_Tmin)
  
  pred_Tmin <- pred_Tmin * rng + mnm 
  
  # Remove simpleton dimension
  
  pred_Tmax <- drop(pred_Tmax)
  
  pred_Tmin <- drop(pred_Tmin)
  
  # LDAPS original (uncorrected forecast)
  
  ldaps_Tmax <- test_data$LDAPS_Tmax_lapse
  
  ldaps_Tmin <- test_data$LDAPS_Tmin_lapse
  
  # Observation 
  
  obs_Tmax <- test_data$Next_Tmax
  
  obs_Tmin <- test_data$Next_Tmin
  
  # Calculating R2, BS, and RMSE
  
  r2_tmax <- R2(pred_Tmax,obs_Tmax)
  
  r2_tmin <- R2(pred_Tmin,obs_Tmin)
  
  ldaps_r2_tmax <- R2(ldaps_Tmax,obs_Tmax)
  
  ldaps_r2_tmin <- R2(ldaps_Tmin,obs_Tmin)
  
  bs_tmax <- BS(pred_Tmax,obs_Tmax)
  
  bs_tmin <- BS(pred_Tmin,obs_Tmin)
  
  ldaps_bs_tmax <- BS(ldaps_Tmax,obs_Tmax)
  
  ldaps_bs_tmin <- BS(ldaps_Tmin,obs_Tmin)
  
  rmse_tmax <- RMSE(pred_Tmax,obs_Tmax)
  
  rmse_tmin <- RMSE(pred_Tmin,obs_Tmin)
  
  ldaps_rmse_tmax <- RMSE(ldaps_Tmax,obs_Tmax)
  
  ldaps_rmse_tmin <- RMSE(ldaps_Tmin,obs_Tmin)
  
  # Appending R2, BS, and RMSE
  
  r2_score_tmax <- append(r2_score_tmax,r2_tmax)
  
  r2_score_tmin <- append(r2_score_tmin,r2_tmin)
  
  ldaps_r2_score_tmax <- append(ldaps_r2_score_tmax,ldaps_r2_tmax)
  
  ldaps_r2_score_tmin <- append(ldaps_r2_score_tmin,ldaps_r2_tmin)
  
  bs_score_tmax <- append(bs_score_tmax,bs_tmax)
  
  bs_score_tmin <- append(bs_score_tmin,bs_tmin)
  
  ldaps_bs_score_tmax <- append(ldaps_bs_score_tmax,ldaps_bs_tmax)
  
  ldaps_bs_score_tmin <- append(ldaps_bs_score_tmin,ldaps_bs_tmin)
  
  print(paste(r2_tmax,r2_tmin,ldaps_r2_tmax,ldaps_r2_tmin))
  
  rmse_score_tmax <- append(rmse_score_tmax,rmse_tmax)
  
  rmse_score_tmin <- append(rmse_score_tmin,rmse_tmin)
  
  ldaps_rmse_score_tmax <- append(ldaps_rmse_score_tmax,ldaps_rmse_tmax)
  
  ldaps_rmse_score_tmin <- append(ldaps_rmse_score_tmin,ldaps_rmse_tmin)
  
  # Appending prediction vs. observed for scatterplot later
  
  prediction_Tmax <- append(prediction_Tmax,pred_Tmax)
  
  prediction_Tmin <- append(prediction_Tmin,pred_Tmin)
  
  observed_Tmax <- append(observed_Tmax,obs_Tmax)
  
  observed_Tmin <- append(observed_Tmin,obs_Tmin)
  
  # Appending original LDAPS forecast for LDAPS performance later 
  
  uncorrected_Tmax <- append(uncorrected_Tmax,ldaps_Tmax)
  
  uncorrected_Tmin <- append(uncorrected_Tmin,ldaps_Tmin)
  
}

# Saving the R2, BIAS, and RMSE score for Next_Tmax and Next_Tmin 
# Saving the predictions and observations for Next_Tmax and Next_Tmin
# Inside a text file 
# 
if (!dir.exists('./results/ann/')){
  dir.create('./results/ann/',recursive=TRUE)
}

savetxt <- function(filename,variable){
  print(filename)
  fid <- file(filename,'w')
    for (x in variable){
      print(format(x))
      writeLines(format(x),fid)
    }
  close(fid)
}

savetxt('./results/ann/ann_r2_score_tmax.txt',r2_score_tmax)

savetxt('./results/ann/ann_bs_score_tmax.txt',bs_score_tmax)

savetxt('./results/ann/ann_rmse_score_tmax.txt',rmse_score_tmax)

savetxt('./results/ann/ann_r2_score_tmin.txt',r2_score_tmin)

savetxt('./results/ann/ann_bs_score_tmin.txt',bs_score_tmin)

savetxt('./results/ann/ann_rmse_score_tmin.txt',rmse_score_tmin)

savetxt('./results/ann/ann_prediction_tmax.txt',prediction_Tmax)

savetxt('./results/ann/ann_prediction_tmin.txt',prediction_Tmin)

savetxt('./results/ann/observed_tmax.txt',observed_Tmax)

savetxt('./results/ann/observed_tmin.txt',observed_Tmin)

savetxt('./results/ann/uncorrected_forecast_tmax.txt',uncorrected_Tmax)

savetxt('./results/ann/uncorrected_forecast_tmin.txt',uncorrected_Tmin)

savetxt('./results/ann/ldaps_r2_score_tmax.txt',ldaps_r2_score_tmax)

savetxt('./results/ann/ldaps_r2_score_tmin.txt',ldaps_r2_score_tmin)

savetxt('./results/ann/ldaps_bs_score_tmax.txt',ldaps_bs_score_tmax)

savetxt('./results/ann/ldaps_bs_score_tmin.txt',ldaps_bs_score_tmin)

savetxt('./results/ann/ldaps_rmse_score_tmax.txt',ldaps_rmse_score_tmax)

savetxt('./results/ann/ldaps_rmse_score_tmin.txt',ldaps_rmse_score_tmin)

# 10-fold cross-validation with XGBoost 
#

library(xgboost)

fold <- 10

N <- nrow(df) / fold

index <- seq(from = 1, to = nrow(df)+N, by = N )

prediction_Tmax <- c()

prediction_Tmin <- c()

observed_Tmax <- c()

observed_Tmin <- c()

r2_score_tmax <- c()

r2_score_tmin <- c()

bs_score_tmax <- c()

bs_score_tmin <- c()

rmse_score_tmax <- c()

rmse_score_tmin <- c()

for (i in 1:(length(index)-1)){
  
  begin <- index[i]
  
  end <- index[i+1]-1

    print(paste('Evaluating test dataset index = ',
              i,
              begin,
              end))

  test_data  <- df[begin:end,]  
  
  train_data <- df[-c(begin:end),]
  
  # XGBoost is decision tree ensemble based model
  
  # No need scaling 
  
  # Setting up X - predictor and y - target 
  
	y_train_tmax <- train_data$Next_Tmax 

	y_test_tmax  <- test_data$Next_Tmax 
	
	X_train_tmax <- select(train_data,-c('Next_Tmax','Next_Tmin','station','Date'))
	
	X_test_tmax  <- select(test_data, -c('Next_Tmax','Next_Tmin','station','Date'))
	
	xgb_train_tmax <- xgb.DMatrix(data=as.matrix(X_train_tmax),label=y_train_tmax)
	
	xgb_test_tmax  <- xgb.DMatrix(data=as.matrix(X_test_tmax), label= y_test_tmax)
	
	y_train_tmin <- train_data$Next_Tmin

	y_test_tmin  <- test_data$Next_Tmin
	
	X_train_tmin <- select(train_data,-c('Next_Tmax','Next_Tmin','station','Date'))
	
	X_test_tmin  <- select(test_data, -c('Next_Tmax','Next_Tmin','station','Date'))
	
	xgb_train_tmin <- xgb.DMatrix(data=as.matrix(X_train_tmin),label=y_train_tmin)
	
	xgb_test_tmin  <- xgb.DMatrix(data=as.matrix(X_test_tmin), label= y_test_tmin)
	
	# Configuring model 
	
	xgb_params <- list( 
		booster = 'gbtree',             # gbtree for tree based models 
		eta = 0.01,                     # learning rate
		gamma = 4,                      # minimum loss reduction to further partition 
		subsample = 0.5,                # subsample ratio prior to growing trees 
		sampling_method = 'uniform',    # each training instance same prob selected 
		objective = 'reg:squarederror', # predicted probability of each class  
		eval_metric = 'rmse'            # root mean square error
	)	
	
	# Training model 
	
	xgb_model_tmax <- xgb.train( 
		params = xgb_params, 
		data = xgb_train_tmax,
		verbose = 2,
		nrounds = 1000
	)
	
	xgb_model_tmin <- xgb.train(
		params = xgb_params,
		data = xgb_train_tmin,
		verbose = 2,
		nrounds = 1000
	)
	
	# Predicting
	
	xgb_preds_tmax <- predict(xgb_model_tmax,as.matrix(X_test_tmax),reshape=TRUE)
	
	xgb_preds_tmin <- predict(xgb_model_tmin,as.matrix(X_test_tmin),reshape=TRUE)
	
	# Calculating R2, BS, and RMSE
	
	r2_tmax <- R2(xgb_preds_tmax,y_test_tmax)
	
	r2_tmin <- R2(xgb_preds_tmin,y_test_tmin)
	
	bs_tmax <- BS(xgb_preds_tmax,y_test_tmax)
	
	bs_tmin <- BS(xgb_preds_tmin,y_test_tmin)
	
	rmse_tmax <- RMSE(xgb_preds_tmax,y_test_tmax)
	
	rmse_tmin <- RMSE(xgb_preds_tmin,y_test_tmin)
	
	# Appending R2, BS, and RMSE
	
	r2_score_tmax <- append(r2_score_tmax,r2_tmax)
	
	r2_score_tmin <- append(r2_score_tmin,r2_tmin)
	
	bs_score_tmax <- append(bs_score_tmax,bs_tmax)
	
	bs_score_tmin <- append(bs_score_tmin,bs_tmin)
	
	print(paste(r2_tmax,r2_tmin))
	
	rmse_score_tmax <- append(rmse_score_tmax,rmse_tmax)
	
	rmse_score_tmin <- append(rmse_score_tmin,rmse_tmin)
	
	# Appending prediction vs. observed for scatterplot later 
	
	prediction_Tmax <- append(prediction_Tmax,xgb_preds_tmax)
	
	prediction_Tmin <- append(prediction_Tmin,xgb_preds_tmin)

	observed_Tmax <- append(observed_Tmax,y_test_tmax)
	
	observed_Tmin <- append(observed_Tmin,y_test_tmin)
	
}

# Saving the R2, BIAS, and RMSE score for Next_Tmax and Next_Tmin 
# Saving the predictions and observations for Next_Tmax and Next_Tmin
# Inside a text file 
# 
if (!dir.exists('./results/boost/')){
  dir.create('./results/boost/',recursive=TRUE)
}

savetxt('./results/boost/xgb_r2_score_tmax.txt',r2_score_tmax)

savetxt('./results/boost/xgb_bs_score_tmax.txt',bs_score_tmax)

savetxt('./results/boost/xgb_rmse_score_tmax.txt',rmse_score_tmax)

savetxt('./results/boost/xgb_r2_score_tmin.txt',r2_score_tmin)

savetxt('./results/boost/xgb_bs_score_tmin.txt',bs_score_tmin)

savetxt('./results/boost/xgb_rmse_score_tmin.txt',rmse_score_tmin)

savetxt('./results/boost/xgb_prediction_tmax.txt',prediction_Tmax)

savetxt('./results/boost/xgb_prediction_tmin.txt',prediction_Tmin)

savetxt('./results/boost/observed_tmax.txt',observed_Tmax)

savetxt('./results/boost/observed_tmin.txt',observed_Tmin)

# 10-fold cross-validation with RandomForests
#

library(randomForest)

fold <- 10

N <- nrow(df) / fold

index <- seq(from = 1, to = nrow(df)+N, by = N )

prediction_Tmax <- c()

prediction_Tmin <- c()

observed_Tmax <- c()

observed_Tmin <- c()

r2_score_tmax <- c()

r2_score_tmin <- c()

bs_score_tmax <- c()

bs_score_tmin <- c()

rmse_score_tmax <- c()

rmse_score_tmin <- c()

for (i in 1:(length(index)-1)){
  
  begin <- index[i]
  
  end <- index[i+1]-1

    print(paste('Evaluating test dataset index = ',
              i,
              begin,
              end))

  test_data  <- df[begin:end,]  
  
  train_data <- df[-c(begin:end),]
  
  # RandomForest (RF) is decision tree ensemble based model
  
  # No need scaling 
  
  # Setting up X - predictor and y - target 
  
  y_train_tmax <- train_data$Next_Tmax 
  
  y_test_tmax  <- test_data$Next_Tmax 
  
  X_train_tmax <- select(train_data,-c('Next_Tmax','Next_Tmin','station','Date'))
  
  X_test_tmax  <- select(test_data, -c('Next_Tmax','Next_Tmin','station','Date'))
  
  y_train_tmin <- train_data$Next_Tmin
  
  y_test_tmin  <- test_data$Next_Tmin
  
  X_train_tmin <- select(train_data,-c('Next_Tmax','Next_Tmin','station','Date'))
  
  X_test_tmin  <- select(test_data, -c('Next_Tmax','Next_Tmin','station','Date'))
  
  
  # Setting up randomForest model and training model 
  
	rf_tmax <- randomForest(X_train_tmax,y_train_tmax,do.trace = TRUE,njobs=12)
					
	rf_tmin <- randomForest(X_train_tmin,y_train_tmin,do.trace = TRUE,njobs=12)
												
	# Predicting
	
	rf_preds_tmax <- predict(rf_tmax,test_data)
	
	rf_preds_tmin <- predict(rf_tmin,test_data)
	
	# Calculating R2, BS, and RMSE
	
	r2_tmax <- R2(rf_preds_tmax,y_test_tmax)
	
	r2_tmin <- R2(rf_preds_tmin,y_test_tmin)
	
	bs_tmax <- BS(rf_preds_tmax,y_test_tmax)
	
	bs_tmin <- BS(rf_preds_tmin,y_test_tmin)

	rmse_tmax <- RMSE(rf_preds_tmax,y_test_tmax)
	
	rmse_tmin <- RMSE(rf_preds_tmin,y_test_tmin)

	# Appending R2, BS, and RMSE
	
	r2_score_tmax <- append(r2_score_tmax,r2_tmax)
	
	r2_score_tmin <- append(r2_score_tmin,r2_tmin)

	bs_score_tmax <- append(bs_score_tmax,bs_tmax)
	
	bs_score_tmin <- append(bs_score_tmin,bs_tmin)
	
	print(paste(r2_tmax,r2_tmin))
	
	rmse_score_tmax <- append(rmse_score_tmax,rmse_tmax)
	
	rmse_score_tmin <- append(rmse_score_tmin,rmse_tmin)
	
	# Appending prediction vs. observed for scatterplot analysis later
	
	prediction_Tmax <- append(prediction_Tmax,rf_preds_tmax)
	
	prediction_Tmin <- append(prediction_Tmin,rf_preds_tmin)
	
	observed_Tmax <- append(observed_Tmax,y_test_tmax)
	
	observed_Tmin <- append(observed_Tmin,y_test_tmin)
	
}

# Saving the R2, BIAS, and RMSE score for Next_Tmax and Next_Tmin 
# Saving the predictions and observations for Next_Tmax and Next_Tmin
# Inside a text file 
# 
if (!dir.exists('./results/randomForest/')){
  dir.create('./results/randomForest/',recursive=TRUE)
}

savetxt('./results/randomForest/rf_r2_score_tmax.txt',r2_score_tmax)

savetxt('./results/randomForest/rf_bs_score_tmax.txt',bs_score_tmax)

savetxt('./results/randomForest/rf_rmse_score_tmax.txt',rmse_score_tmax)

savetxt('./results/randomForest/rf_r2_score_tmin.txt',r2_score_tmin)

savetxt('./results/randomForest/rf_bs_score_tmin.txt',bs_score_tmin)

savetxt('./results/randomForest/rf_rmse_score_tmin.txt',rmse_score_tmin)

savetxt('./results/randomForest/rf_prediction_tmax.txt',prediction_Tmax)

savetxt('./results/randomForest/rf_prediction_tmin.txt',prediction_Tmin)

savetxt('./results/randomForest/observed_tmax.txt',observed_Tmax)

savetxt('./results/randomForest/observed_tmin.txt',observed_Tmin)

	