# Copyright Statement for the R code associated with the manuscript titled
# "An Effective Medium Theory-Based Unified Model for Estimating Thermal Conductivity 
# of Unfrozen and Frozen Soils," published in CATENA.
#
# Author: Hailong Ji
# Affiliation: Nanjing Normal University, China
# Email: jihailongnnu@gmail.com
#
# This R code is provided as supplementary material for the above-mentioned manuscript.
# It is intended for academic and non-commercial use only. Any reproduction, distribution, 
# or modification of this code should be done in compliance with the terms specified 
# by the copyright owner. For any inquiries, please contact the author via the provided email.

# End of Copyright Statement

# calc. conductivity of low-conductivity component, refer to Eq. (5)
f.lambda.lcc<-function(porosity){
  # natural mineral soils
  chi<-0.75
  eta<-1.2
  lambda.lcc<-chi*10^(-eta*porosity)
}


# calc. conductivity of solid particles, refer to Eq. (7)
f.lambda.solid<-function(fr.quartz){
  ifelse(fr.quartz>0.2,lambda.others<-2.0,lambda.others<-3.0)
  lambda.quartz<-7.7
  lambda.solid<-(lambda.quartz^fr.quartz)*(lambda.others^(1-fr.quartz))
}


# calc. conductivity of high-conductivity component, refer to Eq. (6)
f.lambda.hcc<-function(alpha, theta.liq.sat, theta.ice, theta.solid, fr.quartz){
  lambda.liq<-0.56
  lambda.ice<-2.22
  lambda.solid<-f.lambda.solid(fr.quartz)
  
  lambda.hcc<-alpha*(lambda.liq^theta.liq.sat)*
    (lambda.ice^theta.ice)*(lambda.solid^theta.solid)
  return(lambda.hcc)
}


# calc. maximum unfrozen water content, refer to Eqs. (9)-(13)
f.uwc.by.matric<-function(soil.texture,temp){
  # sand, clay, sand, 0~1 
  clay<-soil.texture$clay
  silt<-soil.texture$silt
  sand<-soil.texture$sand
  temp<-ifelse(temp>0,-0.001,temp)  #[degree centigrade]
  
  # calc. psi by Clausius-Clapeyron equation
  latent_fusion<-3.34e+05 #[J/kg]
  rho_liq<-1e+03  #[kg/m3]
  con_fac<-1e+03  #[-] convert [Pa] to [kPa]
  t0<-273.15
  #[kPa]
  psi<-latent_fusion*rho_liq*temp/(con_fac*t0)
  psi<-abs(psi)
  
  A<-exp(-4.396-7.15*clay-4.880*sand^2-4.285*sand^2*clay)*100
  B<-(-3.140)-22.2*clay^2-3.484*sand^2*clay
  # the default base of log(...) in R is natural base
  theta<-exp( log(psi/A)/B )
  return(theta)
}


# calc. effective conductivity based on the explicit form of the GEM equation,
# refer to Eqs. (3)-(8) of Sadeghi et al., 2018, DOI: 10.1002/2017WR021714
f.gem.explicit<-function(soil.meas){
  # simplified variable names
  # three parameters
  fr.cr<-soil.meas$fr.cr
  ts<-soil.meas$ts
  alpha<-soil.meas$alpha
  
  fr.all<-soil.meas$porosity
  # the sum of liquid and ice contents
  fr.hcc<-soil.meas$fr.hcc
  lambda.hcc<-soil.meas$lambda.hcc
  lambda.lcc<-soil.meas$lambda.lcc
  
  ts1<-1/ts
  lambda.hcc1<-lambda.hcc^ts1
  lambda.lcc1<-lambda.lcc^ts1
  
  a1<-fr.cr*lambda.hcc1-(fr.all-fr.cr)*lambda.lcc1
  a1<-a1/(lambda.hcc1-lambda.lcc1)
  
  a2<-fr.all-fr.cr
  a2<-a2/(lambda.hcc1-lambda.lcc1)
  
  a3<-(-1)*fr.cr*lambda.hcc1*lambda.lcc1
  a3<-a3/(lambda.hcc1-lambda.lcc1)
  
  lambda.model<-(-1)*(a1-fr.hcc)+sqrt((a1-fr.hcc)^2-4*a2*a3)
  lambda.model<-lambda.model/(2*a2)
  lambda.model<-lambda.model^ts
  return(lambda.model)
}


f.unified.unfrozen<-function(soil.meas){
  ### step 1: process the soil properties
  # only for frozen soil?
  #soil.meas<-soil.meas[soil.meas$temperature<0,]
  # refer to Eq. (8)
  soil.meas$quartz<-ifelse(is.null(soil.meas$quartz),
                           0.5*soil.meas$sand,soil.meas$quartz)  
  # simplified variable names
  sand<-soil.meas$sand
  fr.quartz<-soil.meas$quartz
  porosity<-soil.meas$porosity

  
  ### step 2: calculate the conductivity of high-conductivity component (hcc)
  # conductivity of high-conductivity component
  # alpha=1.028 is determined from the training dataset, which can be tuned for specific sample
  soil.meas$lambda.hcc<-f.lambda.hcc(alpha=1.028, theta.liq.sat=porosity, 
                                     theta.ice=0, theta.solid=1-porosity, 
                                     fr.quartz=fr.quartz)
  
  
  ### step 3: calculate the fraction of hcc
  soil.meas$fr.hcc<-soil.meas$volume.water.content
  soil.meas$fr.hcc<-ifelse(soil.meas$fr.hcc<=soil.meas$porosity,
                           soil.meas$fr.hcc,soil.meas$porosity)
  
  ### step 4: calculate the conductivity of low-conductivity component (lcc)
  soil.meas$lambda.lcc<-f.lambda.lcc(soil.meas$porosity)
  
  ### step 5: calculate the effective conductivity
  # estimate the parameters of the unified model through the built PTFs
  soil.meas$fr.cr<-0.46*porosity-0.16 # refer to Eq. (26)
  soil.meas$ts<-(-0.18)*sand+0.44 # refer to Eq. (27)
  soil.meas$alpha<-1.028 # refer to Eq. (28)
  soil.meas$unified<-f.gem.explicit(soil.meas)
  return(soil.meas)
}


# 
f.unified.frozen<-function(soil.meas){
  ### step 1: process the soil properties
  # only for frozen soil?
  #soil.meas<-soil.meas[soil.meas$temperature<0,]
  # refer to Eq. (8)
  soil.meas$quartz<-ifelse(is.null(soil.meas$quartz),
                           0.5*soil.meas$sand,soil.meas$quartz)
  # simplified variable names
  fr.quartz<-soil.meas$quartz
  porosity<-soil.meas$porosity
  temp<-soil.meas$temperature
  theta.ini<-soil.meas$volume.water.content

  
  ### step 2: calculate the conductivity of high-conductivity component (hcc)
  # the maximum unfrozen water content (UWC) at -40°„C (assumed to be completely frozen)
  uwc.max.40<-f.uwc.by.matric(soil.texture=soil.meas[,c("sand","silt","clay")],
                              temp = -40)
  # the ice content at saturated soil at -40°„C
  ice.con.40<-1.09*(porosity-uwc.max.40)
  
  # conductivity of high-conductivity component
  # alpha=1.001 is determined from the training dataset, which can be tuned for specific sample
  soil.meas$lambda.hcc<-f.lambda.hcc(alpha=1.001, theta.liq.sat=uwc.max.40, 
                             theta.ice=ice.con.40, theta.solid=1-porosity, 
                             fr.quartz=fr.quartz)
  
  
  ### step 3: calculate the fraction of hcc
  # the maximum UWC at a given temperature
  uwc.max<-f.uwc.by.matric(soil.texture = soil.meas[,c("sand","silt","clay")],
                                      temp = temp)
  # refer to Eq. (14)
  soil.meas$liquid.content<-ifelse(theta.ini>=uwc.max, uwc.max, theta.ini)
  soil.meas$ice.content<-1.09*(theta.ini-soil.meas$liquid.content)
  # fraction of hcc
  soil.meas$fr.hcc<-soil.meas$ice.content+soil.meas$liquid.content
  soil.meas$fr.hcc<-ifelse(soil.meas$fr.hcc<=soil.meas$porosity,
                            soil.meas$fr.hcc,soil.meas$porosity)
  
  ### step 4: calculate the conductivity of low-conductivity component (lcc)
  soil.meas$lambda.lcc<-f.lambda.lcc(soil.meas$porosity)
  
  ### step 5: calculate the effective conductivity
  # estimate the parameters of the unified model through the built PTFs
  soil.meas$fr.cr<-1.01*theta.ini-0.01 # refer to Eq. (26)
  soil.meas$ts<-0.23 # refer to Eq. (27)
  soil.meas$alpha<-1.01 # refer to Eq. (28)
  soil.meas$unified<-f.gem.explicit(soil.meas)
  return(soil.meas)
}

#########################################################################################
# We do not provide the complete soil thermal conductivity (STC) measurement dataset due to the copyright restrictions,
# while the interested readers can digitize them from the references listed in Table 1.
# Here, only two test cases for unfrozen and frozen soils are provided.
# Note, the three parameters are estimated by the built PTFs, 
# which might not perform well for the specific case as discussed in the manuscript.

# Description of the measurement dataset
# where <M> denotes that the field must be prepared before running the script
#
# [code]: specific to a sample, with the naming format of "Author-Year-Soil name-Dry bulk density"
# [clay],[silt],[sand]: <M>, represented by 0-1, where their sum is zero, (m3/m3)
# [quartz]: represented by 0-1, which can be estimated by Eq. (8), (m3/m3)
#
# [particle.density]: <M>, solid particle density, (kg/m3)
# [dry.bulk.density]: <M>, dry bulk density, (kg/m3)
# [porosity]: <M>, porosity, (m3/m3)
# note 1: [porosity]=1-[dry.bulk.density]/[particle.density]
#
# [temperature]: <M>, measurement temperature, which is not used in unfrozen soils, (°„C)
# [volume.water.content]: <M>, akin to initial water content for frozen soils, which is measured under unfrozen state, (m3/m3)
# [liquid.content]: equals to [volume.water.content] for unfrozen soils,
#                   missing for most frozen samples, which can be estimated by Eqs. (9)-(14), (m3/m3)
# [thermal.conductivity]: <M>, measured STC, (W/(m K))
# [instrument]: the measurement method used
#
#
#
# After run the following script, 7 variables regarding the model will be appended in the dataset
# [lambada.hcc], the estimated conductivity of high-conductivity component (hcc), (W/(m K))
# [fr.hcc], the fraction of hcc, (m3/m3)
# [lambda.lcc], the estimated conductivity of high-conductivity component (hcc), (W/(m K))
# [fr.lcc], the fraction of lcc, (m3/m3)
# [ts], scaling exponent
# [alpha], compensating factor
# [unified], modeled STC

# test case 1: unfrozen soil sample
# the measurement comes from Lu et al., 2007, DOI:10.2136/sssaj2006.0041
load("samples.unfrozen.test.RData")
result.unified.unfrozen.test<-f.unified.unfrozen(samples.unfrozen.test)

# test case 2: frozen soil sample
# the measurement comes from Xu et al., 2020, DOI: 10.1016/j.rinp.2019.102830
load("samples.frozen.test.RData")
result.unified.frozen.test<-f.unified.frozen(samples.frozen.test)

# plot the measured and modeled STC 
# require(ggplot2)
# ggplot(result.unified.unfrozen.test)+
#   geom_point(aes(fr.hcc,thermal.conductivity,color="meas."))+
#   geom_point(aes(fr.hcc,unified,color="model"))+
#   ggtitle(result.unified.unfrozen.test$code)
# 
# ggplot(result.unified.frozen.test)+
#   geom_point(aes(fr.hcc,thermal.conductivity,color="meas."))+
#   geom_point(aes(fr.hcc,unified,color="model"))+
#   ggtitle(result.unified.frozen.test$code)