
install.packages('devtools')
devtools::install_github('daniel1noble/metaDigitise', force=TRUE, dependencies=TRUE)
library(metaDigitise)

df<-metaDigitise('C:/Users/z5288536/OneDrive - UNSW/[ PhD UNSW ]/Side project - Thermal fertility experimental evolution/Evol_reproduction_temp/data_extraction/figures') 
