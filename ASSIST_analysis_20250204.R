# ASSSIT March 2024 FINAL analysis
# 23/1/25
# R 4.4

library(lme4)
library(pbkrtest)
library(glmmTMB) # for neg bin models
library(ggplot2)
library(DHARMa)
library(emmeans)
library(MASS)
library(ape)
library(lmerTest)
library(fmsb)
library(nlme)
library(tidyverse)
library(bestNormalize)
library(dplyr)
library(car) 
#----------------------------------------

# code for model simpification in pbkrtest
#https://www.rdocumentation.org/packages/lme4/versions/1.1-34/topics/drop1.merMod

KRSumFun <- function(object, objectDrop, ...) {
  krnames <- c("ndf","ddf","Fstat","p.value","F.scaling")
  r <- if (missing(objectDrop)) {
    setNames(rep(NA,length(krnames)),krnames)
  } else {
    krtest <- KRmodcomp(object,objectDrop)
    unlist(krtest$stats[krnames])
  }
  attr(r,"method") <- c("Kenward-Roger via pbkrtest package")
  r
}


# drop1(gm1,test="F",sumFun=KRSumFun)
#------------------------------------------------------------------------------------
setwd("P:\\NEC05829 LTS-M-ASSIST\\WP3_Woodcock\\2. Data\\1. Main arable field experiment\\Paper1_data_and_analysis\\Woodcock_ASSIST_agroecological_enhancement\\")
data_field_mean_T1to3 <- read.csv("Final_data_csv.csv")

# Description of variables in Final_data_csv.csv
#---------------------------
# Year_since_est	 =  Categorical years since the establishment year -  note most sites established 2018, two sites 2019
# Crop	=  Crop type (WW=winter wheat, OSR=oilseed, SB=spring barley, SW=spring wheat, SB=spring barley)
# Crop_Type	 = General crop category (Winter-winter sown cereal; OSR= winter oilseed;Spring=spring sown cereal)
# Winter_Spring_crop = 	Identifies spring or winter sown crops
# Cereal	= identifies crop as cereal or as oilseed
# Site	= site classifier
# Year_callender	= calender year
# Treat	= treatment number, where 1=BAU control, 2= Enhancing-ES, 3= Maximising-ES
# Field_id	= unique field identifier
# Crop_Failed	= identifies three field for particular years where the crop failed 
# Precision_Ag	= identifies sites where precision agriculture yield data is avalliable -  this may only be on certain field in a site
# Seed_weight_0.25	= quadrat based crop seed weight in a 0.25 m quadrat
# Precision_Ag_whole_site	= identifies sites where all three fields for a year have precision yield data
# PrecAg_yield_t_ha	=  Combine reported precision agriculture yield data tonnes ha-1
# CoverCrop_inthepreceedingyear	 =  identifies field and years where a cover crop was grown prior to the spring crop
# FYM	Margin_area_ha	= the total area of wild flower sown field margins in a field
# in_field_strip_area_ha	= the total area of sown wild flower in-field strips in a field.
# Crop_value_GBP_tonne_FOR_CALCULATION =  Average value of crop £ tonnes -  speific to each crop	
# FM_yield_reduction_at_edges	= The proportional reduction of yield at field edges 
# Treat_2_3	= identifies if field is in treatment 2 and 3 (otherwise 0)
# Treat_3	= identifies if field is in treatment  3 (otherwise 0)
# Field_total_Area_ha	= field total area ha
# cropped_area_ha	= field total cropped area in ha (i.e. some field s have field margins and in-field strips)
# Worms_Anecic_ave	= average abundance of worms
# blackgrass_ave	= average count of blackgrass stems
# Weeds_ave	= average count of weed plants
# FM_sown_forb_SR_ave	= average species richness of the sown component of the wild flower field margins
# Drop_disk_ave	= sward structure in the field margins measrued using a dropdisk
# Aphid_cards_plant_aveN_eaten	= predation assessments using aphid cards, average numbers eaten
# Slug_N_LargeBeetle_bite =  average numer fo artifical slugs btiten per field
# Slug_biomass	= average slug biomass
# Snail_biomass	= average snail bioimass
# Parasitized_aphids_hand_search	= average number of parasitised aphids found by hand searching
# Aphids_hand_search	= average number of aphids found by hand searching
# Predaors_hand_search_sum	= average number of crop canopy predators found by hand searching
# Carabidae_Pitfall_N_mean	= average abundance of ground beetles colleted from pitfall traps
# spider_Pitfall_N_mean	= average abundance of spiders colleted from pitfall traps
# staph_Pitfall_N_mean	= average abundance of rove beetles colleted from pitfall traps
# Parasitoid_OSR_N	= average abundance of parasitoids from suction samples identifed as feeding on oilseed pest species
# Parasitoid_Cereal_N	= average abundance of parasitoids from suction samples identifed as feeding on ceral crop pest species
# Bee_N_margin_m2	= average density of bees in field margin areas
# Parasitoid_N_margin_m2		= average density of parasitoids in field margin areas
# Hoverfly_N_margin_m2		= average density of hoverflies in field margin areas
# Shannon_landuse	= shannon land use diversity in the 2km surrounding fields
# Pollination_assessment_OSR	=  Identifies those crops of oilseed rape where a pollination effectiveness assessment was made
# Yield_from_poll	= yield of oilseed rape attributable to pollinators in gramms
# SMD_Ave_value_crop_ha_noAESpayment	= Standard mean difference in profit (within site and for each year) assuming no subsidies  
# SMD_Ave_value_crop_ha_withAESpayment	= Standard mean difference in profit (within site and for each year) assuming AES subsidies  
# SMD_Yield_Tonnes_ha_corrected = Standard mean difference in yield (within site and for each year)   


data_field_mean_T1to3$Year_callender<-as.factor(data_field_mean_T1to3$Year_callender) 
data_field_mean_T1to3$Year_since_est<-as.factor(data_field_mean_T1to3$Year_since_est) 
data_field_mean_T1to3$Treat<-as.factor(data_field_mean_T1to3$Treat) 

# Main overall data set excluding sites with missing yield data (fields that failed in a year)
# Note -  these fields could arguably have been given a 'zero' yield (they failed) however theya re typically resown with something (e.g. a spring bean)
# and so are not technically no yield.  In this context it was considered sensible to just exclude these failed field data points.

data_field_mean_T1to3 <- data_field_mean_T1to3 %>% filter(Crop_Failed == 0)    #  subset data excluding sites where a whole fields crop failed in a year.




####################################################################################
####################  Deriving key metrics for subsequent analysis   ###############
####################################################################################

#####################################################################################
#  Measures of field area, cropped area and area fo field margins and in-field strips

data_field_mean_T1to3$Field_total_Area_ha	# total area of the field in in hectares (including field margins and in-field strips)
data_field_mean_T1to3$cropped_area_ha   # total area where crops are grown (for T2 and T3 this is less than data_field_mean_T1to3$Field_total_Area_ha	)
data_field_mean_T1to3$Margin_area_ha  # area of sown wild flower field margins in a field	
data_field_mean_T1to3$in_field_strip_area_ha	# area of infield strips in a field
data_field_mean_T1to3$gree_infastrucutre_area_ha  <- data_field_mean_T1to3$in_field_strip_area_ha +
  data_field_mean_T1to3$Margin_area_ha # total area of green infrastructure (infield strips and field margins)	
data_field_mean_T1to3$ratio_FM_to_crop	<-data_field_mean_T1to3$Margin_area_ha / data_field_mean_T1to3$Field_total_Area_ha#  ratio of field margin area to total field area
data_field_mean_T1to3$Ratio_IFS_to_crop <-data_field_mean_T1to3$in_field_strip_area_ha/ data_field_mean_T1to3$Field_total_Area_ha #  ratio of in-field strip area to total field area
data_field_mean_T1to3$ratio_cGI_to_crop <-(data_field_mean_T1to3$in_field_strip_area_ha+data_field_mean_T1to3$Margin_area_ha)/ data_field_mean_T1to3$Field_total_Area_ha

#############################################################################################
# Converting average quadrat based yields (0.25 x 0.25 m) in terms of seed weight in grams (Seed_weight_0.25)
# into quadrat predicted yield in tonnes ha based on the cropped area of the field (i.e. the area with crop on it excluding the fields)

data_field_mean_T1to3$Yeild_Tonnes_ha_quadrat_extrapolated<-(data_field_mean_T1to3$Seed_weight_0.25*4)/100

######################################################################################
# Q1:  What is the relationship between predicted yield for fields based on quadrats
# and the actual precision agriculture yield.  This tests whether the directly measured quadrat based
# crop yields are a robust measure of what is seen by farmers.

#subset data for those sites with precision agriculture data that can be related to the 
# quadrat based yield estimates
Precision_ag_data <- data_field_mean_T1to3 %>% filter(Precision_Ag == 1)

# plot showing the relationship between the precision ag yield data and the quadtrat based data
ggplot(Precision_ag_data, aes(x = Yeild_Tonnes_ha_quadrat_extrapolated, y = PrecAg_yield_t_ha)) +
  geom_point(size=4) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +  # Regression line with ribbon
  labs( x = "Quadrat based yield (tonnes / ha)", y = "Combine harvester yield (tonnes / ha.)", color = "System") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 18, face = "bold"),  # Change axis titles font size
    axis.text = element_text(size = 14),   # Change axis labels font size
  )

# relationship between precision agric. yield data (tonnes ha from the cropped area)
# and quadrat based yield estimates (extrapolated from quadrats to tonne ha)
Prec_yield_treat_response<-lm(PrecAg_yield_t_ha ~  Yeild_Tonnes_ha_quadrat_extrapolated  
                              , data=Precision_ag_data)
drop1(Prec_yield_treat_response,test="F",sumFun=KRSumFun)
sim_res <- simulateResiduals(fittedModel = Prec_yield_treat_response)  #  DHARMa tests of goodness of fit
plot(sim_res)   #  meets assumptions of model
# significant relationship
summary(Prec_yield_treat_response)
# extract model parameter estimates
estimates <- coef(Prec_yield_treat_response)
intercept <- estimates["(Intercept)"]
slope <- estimates["Yeild_Tonnes_ha_quadrat_extrapolated"]

#corrected yield (tonnes ha-1 based on cropped area) based on this correlation.
# This adjusts the underestimated quadrat based yield to be in line with that of the precision agriculture measures of 
# average field yield per ha.   Quadrat based yield underestimates precision yield a bit.
# This yield is the basis for the econoimc assessments.

data_field_mean_T1to3$Yield<-  (data_field_mean_T1to3$Yeild_Tonnes_ha_quadrat_extrapolated*slope)+intercept


#  However, another yield corrected for the total area of the field can be derived (Yield_absolute tonnes ha, based
# on the total field area), i.e. not just for the cropped area of which 'Yield' is based on.  This measure can be used to 
# infer whether the use of management practices in T2 and T3 (field margins, infield strips.
# cover crops and manure) increase the total productivity of the land per unit area of
# whole fields.  By definition this measure of yield would be identical to that of 'Yield'
#  for T1 (where no land is taken out of production) but lower than 'Yield' in T2 and T3 
# where land is allocated to the non-cropped areas.

data_field_mean_T1to3$Yield_absolute<- (data_field_mean_T1to3$Yield * data_field_mean_T1to3$cropped_area_ha)  /
  data_field_mean_T1to3$Field_total_Area_ha


################################################################################################################
# Precision Agriculture harvested yield (subset of fields and sites with this data)
#filter overall data set for sites that have precision ag yields for all three fields at a site for a given year

Precision_ag_data_whole_site_only <- data_field_mean_T1to3 %>% filter(Precision_Ag_whole_site == 1)

# In general what is the trend of precision yield with the treatments
Prec_yield_treat_response<-lmer(PrecAg_yield_t_ha ~
                                  Treat +
                                  Year_since_est+
                                  Crop +
                                  Treat*Crop+
                                  Treat*Year_since_est +
                                  Crop*Year_since_est +
                                  Year_since_est*Treat*Crop +
                                  (1|Site/Year_since_est), data=Precision_ag_data_whole_site_only,  na.action = "na.fail", REML = FALSE)


#....................................................
# deletion of least sig effects
Prec_yield_treat_response_m2<-lmer(PrecAg_yield_t_ha ~
                                     Treat +
                                     Crop +
                                     (1|Site/Year_since_est), data=Precision_ag_data_whole_site_only,  na.action = "na.fail")
drop1(Prec_yield_treat_response_m2,test="F",sumFun=KRSumFun)


#DHARMAa
simulationOutput<-simulateResiduals(fittedModel=Prec_yield_treat_response_m2, plot=F)
plot(simulationOutput)   # meets model assumptions

summary(Prec_yield_treat_response_m2)






###########################################################################################
# Yield (tonnes ha)  -  average for cropped area
#  This is the yield is derived for the cropped area.
# Note this is a derived measure using a combination of quadrat based yield estimates, and the 
# relationship between precision agric yield measure and quadrat yield acting as a correction
#  This is shown to be highly correlated (R2=0.58 with precision ag reported field yields)
# see above for derivation

yield<-lmer(Yield ~
              Treat +
              Year_since_est+
              Crop +
              Treat*Crop +
              Treat*Year_since_est +
              Crop*Year_since_est +
              Year_since_est*Treat*Crop +
              (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 



#min adeq model
yield<-lmer(Yield ~
              Treat +
              Year_since_est+
              Crop +
              (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(yield,test="F",sumFun=KRSumFun)


simulationOutput<-simulateResiduals(fittedModel=yield, plot=F)
plot(simulationOutput) # model meets assumptions
summary(yield)

ggplot(data=data_field_mean_T1to3, aes(y=SMD_Yield_Tonnes_ha_corrected, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
    geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = "Management system",                                           # X-axis legend
    y = "Yield (SMD Tonnes / ha.)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text = element_text(size = 12),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )    




###########################################################################################
# Yield_absolute (tonnes ha)  -  average for total area of field (including cropped and non-cropped areas)
#  This is the yield is derived for the cropped area.
# Note this is a derived measure using a combination of quadrat based yield estimates, and the 
# relationship between precision agric yield measure and quadrat yield acting as a correction
#  This is shown to be highly correlated (R2=0.58 with precision ag reported field yields)
# see above for derivation

yield_ab<-lmer(Yield_absolute  ~
                 Treat +
                 Year_since_est+
                 Crop +
                 Treat*Crop +
                 Treat*Year_since_est +
                 Crop*Year_since_est +
                 Year_since_est*Treat*Crop +
                 (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 



#min adeq model
yield_ab2<-lmer(Yield_absolute  ~ 
                  Year_since_est+
                  Crop +
                  (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(yield_ab2,test="F",sumFun=KRSumFun)
simulationOutput<-simulateResiduals(fittedModel=yield_ab2, plot=F)
plot(simulationOutput) # model meets assumptions

summary(yield_ab2)



#######################################################################################################
# Value of crop

# crop value in £ tonne
# Average value from 2018-2021 based on published AHDB prices used as guide for farmers
# https://ahdb.org.uk/cereals-oilseeds/uk-delivered-prices - 
# https://ahdb.org.uk/cereals-oilseeds/ex-farm-prices-summary	  
data_field_mean_T1to3$Crop_value_GBP_tonne_FOR_CALCULATION    # crop value in £ tonnes  

#  Costs to establish field margins and in-field strips
Plough <-	81.16	# NAAC 2023 contractor costs £ ha
Press<-	50.52 #  GBP ha seedbed (pressing)	NAAC 2023 contractor costs
spread<-	33.06 #  GBP £ ha Broadcast seed	NAAC 2023 contractor costs
Ringroll<-	25.56 # GBP £ ha ring rolling of seeds once sown NAAC 2023 contractor costs
cutting<-50.85  #GBP £ ha  NAAC 2023 contractor costs for a single cut of vegetation
# note cut 3x yr 1, 2x y2, then 1x annually
GreenInfastr_AES_payment_GBP_ha<- 798.00	  #  Payments to farmers in £ by government to establish field margins and in-field strips
GI_seed_mix_GBP_ha <- 770.35#  costs to farmers in £ per ha for the seed mix to establish field margins and in-field strips

# green infrastructure costs for establishment only per year (over a 4 year period)
GreenInfaStr_Estab_cost_GBP_ha_per_year<- -(GI_seed_mix_GBP_ha+Plough+Press+spread+Ringroll+(cutting*7))/4


#Field margin yield reduction at edge 
# different values for each cereal and oilseed crops
# See Fincham, W.N.W., Redhead, J.W., Woodcock, B.A., Pywell, R.F., 2023. Exploring drivers of within‐field crop yield variation using a national precision yield network. J. Appl. Ecol., 319–329.
data_field_mean_T1to3$FM_yield_reduction_at_edges

#Cost of field margins in terms of lost crop  (FIELD TOTAL  GBP)
data_field_mean_T1to3$FM_cost_lost_crop_GBP_field_total<- -( data_field_mean_T1to3$Yield  *
                                                               data_field_mean_T1to3$Margin_area_ha    *
                                                               data_field_mean_T1to3$FM_yield_reduction_at_edges *
                                                               data_field_mean_T1to3$Crop_value_GBP_tonne_FOR_CALCULATION)


#Cost of infield strips in terms of lost crop  (FIELD TOTAL  GBP)
data_field_mean_T1to3$IFS_cost_lost_crop_GBP_field_total<- -( data_field_mean_T1to3$Yield  *
                                                                data_field_mean_T1to3$in_field_strip_area_ha	    *
                                                                data_field_mean_T1to3$Crop_value_GBP_tonne_FOR_CALCULATION)

# Cost of cover crops 
Cover_crops_GBP_ha<-57.20   #  cost of seed mif 3 per ha
# based on 2023 NAAC rates  GBP £ ha
cultivation_seebox<- 40 # Single pass light cultivation with seed box
spraying_costs<- 12.54  #  Spraying off cover crops prior to main spring crop cultivation.  Cost for single pass spraying machinery only.
Glyphosate<- 16.8 #  Glyphosate costs for spraying off crop
Cover_crop_est_cost_GBP_ha_ave_4yrs   <-  -(Cover_crops_GBP_ha+cultivation_seebox+spraying_costs+Glyphosate)/4   # cost GBP ha to establish cover crops (excluding offset AES payment) averaged over 4 years


# Agri-environemtn payment for cover crops
covercrop_aes<- 129.00   # GBP  ha payment to farmers only for year with cover crops
Cover_Crop_AES_payment_GBP_ha_ave_4yrs  <- covercrop_aes/4  # assuming one cover crop in 4 years annual


# Farm yard manure costs
FYM_GBP_tonne<- 3.00   # £ cost per tonne
FYM_GBP_ha<- 30 *  FYM_GBP_tonne  # cost GBP per ha assuming 30 tonnes per ha
Delivery<-  15 #  based on 30t load at £1.50/mile and a delivery distance of 10 miles.
Spreading_FYM<- 97.20  # 2023 NAAC rates for (30 tonnes at £3.24 tonne)
Offset_fertilsier_value<- 64.62   #  Based on AHDB 4 figures cattle FYM contains 6 kg N (£1.10 kg), 1.9 kg P (£0.95 kg) and 8.5 kg K (£0.70 kg).  The amount of N P and K available to the crop in the first year is typically 10% of its content in FYM for the first crop and 5% for the following crop 4 .
#  +£64.62 (@ £1.44 per tonne FYM for the first crop and £0.72 for the second crop)
FYM_price_ha_productANDdelivery <- (FYM_GBP_ha+Delivery)/4 #  cost of purchasing and delivering FYM per ha 
FYM_costs_ha_averaged_over4yrs<- (Offset_fertilsier_value-(FYM_GBP_ha+Delivery+Spreading_FYM))/4   # cost of applying FYM , including offset relating to its fertilsier value in GBP ha


# Economic value of crop  (GBP ha-1)
# Note these values are based on the total field area (cropped area + non -cropped field margin /  infield strip area)

data_field_mean_T1to3$Treat_2_3	# binomial descriptor of whether the field is in treatment 2 and 3 and so will have field margins and cover crops

data_field_mean_T1to3$Treat_3	# binomial descriptor for Treatment 3 which also has in-field strips and farm year manure

data_field_mean_T1to3$Tonnes_per_field <-	data_field_mean_T1to3$cropped_area_ha * data_field_mean_T1to3$Yield #  total crop yield tonnes per field

data_field_mean_T1to3$Crop_value_GBP_fieldtotal	<-data_field_mean_T1to3$Tonnes_per_field *  data_field_mean_T1to3$Crop_value_GBP_tonne_FOR_CALCULATION # total value of that crop £ for that field

data_field_mean_T1to3$Management_costs_GBP_field_total <-	# total management costs £ for that field
  (data_field_mean_T1to3$FM_cost_lost_crop_GBP_field_total + data_field_mean_T1to3$IFS_cost_lost_crop_GBP_field_total)+   # cost of forgone yield where field margin aor infield strips are
  (data_field_mean_T1to3$gree_infastrucutre_area_ha*GreenInfaStr_Estab_cost_GBP_ha_per_year)+   # cost per field of field margin and infield strip establishment
  (data_field_mean_T1to3$Treat_2_3  * Cover_crop_est_cost_GBP_ha_ave_4yrs * data_field_mean_T1to3$cropped_area_ha)+  # cost of cover crops
  (data_field_mean_T1to3$cropped_area_ha * data_field_mean_T1to3$Treat_3 *FYM_costs_ha_averaged_over4yrs)   # farm yard manure costs whole field

data_field_mean_T1to3$AES_payment_total_field	<-  #  total AES payments for that field (for field margins, infield strips and cover crops where present)
  (data_field_mean_T1to3$gree_infastrucutre_area_ha *GreenInfastr_AES_payment_GBP_ha)+
  (Cover_Crop_AES_payment_GBP_ha_ave_4yrs *data_field_mean_T1to3$cropped_area_ha *data_field_mean_T1to3$Treat_2_3)	

data_field_mean_T1to3$Ave_value_crop_ha_noAESpayment <-  # Value of crop £ ha (no AES payment) in a year BASED ON TOTAL FIELD AREA (cropped area + non-cropped AES area)
  (data_field_mean_T1to3$Crop_value_GBP_fieldtotal+data_field_mean_T1to3$Management_costs_GBP_field_total )/
  data_field_mean_T1to3$Field_total_Area_ha

data_field_mean_T1to3$Ave_value_crop_ha_withAESpayment <-	# Value of crop £ ha (with AES payment) in a year BASED ON TOTAL FIELD AREA (cropped area + non-cropped AES area)
  (data_field_mean_T1to3$Crop_value_GBP_fieldtotal+data_field_mean_T1to3$Management_costs_GBP_field_total + data_field_mean_T1to3$AES_payment_total_field	)/
  data_field_mean_T1to3$Field_total_Area_ha

data_field_mean_T1to3$Ave_value_crop_ha_noAESpayment_FYMFREE <-	# Value of crop £ ha (no AES payment but free manure) in a year BASED ON TOTAL FIELD AREA (cropped area + non-cropped AES area)
  data_field_mean_T1to3$Ave_value_crop_ha_noAESpayment +
  ((data_field_mean_T1to3$Treat_3 * FYM_price_ha_productANDdelivery *  data_field_mean_T1to3$cropped_area_ha  )
   / data_field_mean_T1to3$Field_total_Area_ha)


data_field_mean_T1to3$Ave_value_crop_ha_withAESpayment_FYMFREE <- # Value of crop £ ha (with AES payment and free manure) in a year BASED ON TOTAL FIELD AREA (cropped area + non-cropped AES area)
  data_field_mean_T1to3$Ave_value_crop_ha_withAESpayment +
  ((data_field_mean_T1to3$Treat_3 * FYM_price_ha_productANDdelivery *  data_field_mean_T1to3$cropped_area_ha  )
   / data_field_mean_T1to3$Field_total_Area_ha)




#########################################################################################
#  Profit (excluding AES payments)

Profit_noAES<-lmer(Ave_value_crop_ha_noAESpayment ~
                     Treat +
                     Year_since_est+
                     Crop +
                     Treat*Crop +
                     Treat*Year_since_est +
                     Crop*Year_since_est +
                     Year_since_est*Treat*Crop +
                     (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 


# model checks


#................................
#min adeq model
Profit_noAES<-lmer(Ave_value_crop_ha_noAESpayment ~
                     Treat +
                     Year_since_est+
                     Crop +
                     (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail") 
drop1(Profit_noAES,test="F",sumFun=KRSumFun) 
summary(Profit_noAES)
emmeans(Profit_noAES, ~ Treat)



simulationOutput<-simulateResiduals(fittedModel=Profit_noAES, plot=F)
plot(simulationOutput)  # Model meets assumptions


ggplot(data=data_field_mean_T1to3, aes(y=Ave_value_crop_ha_noAESpayment, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = "Management system",                                           # X-axis legend
    y = "Profit (£ / ha. / year)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text = element_text(size = 12),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )    


#######################################################################################
# Profit including the AES payments


Profit_withAES<-lmer(Ave_value_crop_ha_withAESpayment ~
                       Treat +
                       Year_since_est+
                       Crop +
                       Treat*Crop +
                       Treat*Year_since_est +
                       Crop*Year_since_est +
                       Year_since_est*Treat*Crop +
                       (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 





#.........................
#min adeq model

Profit_withAES<-lmer(Ave_value_crop_ha_withAESpayment ~
                       Treat +
                       Year_since_est+
                       Crop +
                       (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail") 


drop1(Profit_withAES,test="F",sumFun=KRSumFun)                       

summary(Profit_withAES)
emmeans(Profit_withAES, ~ Treat)

simulationOutput<-simulateResiduals(fittedModel=Profit_withAES, plot=F)
plot(simulationOutput)   # Model meets assumptions


ggplot(data=data_field_mean_T1to3, aes(y=Ave_value_crop_ha_withAESpayment, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = "Management system",                                           # X-axis legend
    y = "Profit (£ / ha. / year)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text = element_text(size = 12),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )  

###################################################################################
#  ES service and biodiversity  measures

# basic model in terms of fixed and random effects

m1<-lmer( XXXXXXXX   ~
            Treat +
            Year_since_est+
            Treat*Year_since_est +
            (1|Site/Year_since_est), data=data_field_mean_T1to3) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)

####################################################################################

#Aphid_cards_plant_aveN_eaten
# average number iof aphids eaten from cards

aphid_card_data<-data_field_mean_T1to3 %>%   filter(!is.na(Aphid_cards_plant_aveN_eaten))
is.numeric(aphid_card_data$Aphid_cards_plant_aveN_eaten)  #  TRUE
any(is.na(aphid_card_data$Aphid_cards_plant_aveN_eaten))
max(aphid_card_data$Aphid_cards_plant_aveN_eaten)
min(aphid_card_data$Aphid_cards_plant_aveN_eaten)

# average numbers eaten so continuous -  use Gaussian model first

m1<-lmer(Aphid_cards_plant_aveN_eaten~
           Treat +
           Year_since_est+
           (1|Site/Year_since_est), data=aphid_card_data) 
drop1(m1,test="user",sumFun=KRSumFun)
summary(m1)


simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # some deviation detected 

# try log transformation

m1<-lmer(log(Aphid_cards_plant_aveN_eaten+1)~
           Treat +
           Year_since_est+
           (1|Site/Year_since_est), data=aphid_card_data) 
drop1(m1,test="user",sumFun=KRSumFun)
summary(m1)


simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # Model meets assumptions

ggplot(data=aphid_card_data, aes(y=Aphid_cards_plant_aveN_eaten, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  coord_cartesian(ylim = c(0, 2.5))  +# Adjusts view without removing data

  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Aphids eaten (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text = element_text(size = 16),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )  



#---------------------------------------------------------------------------------------
#Slug_aveN_LargeBeetle_bite
# average number of slugs attacked by ground beetles, other large beetles
# per sampling point.  Continuous.

m1<-lmer( Slug_aveN_LargeBeetle_bite~
            Treat +
            Year_since_est+
            Treat*Year_since_est +
            (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # some deviations  -  try transformation to deal with


m1<-lmer( log(Slug_aveN_LargeBeetle_bite+1)~
            1 +
            (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # OK


summary(m1)

ggplot(data=data_field_mean_T1to3, aes(y=Slug_aveN_LargeBeetle_bite, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Slugs bitten (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  


#---------------------------------------------------------------------------------------
#Parasitized_aphids_hand_search
m1<-lmer(Parasitized_aphids_hand_search ~
           Year_since_est+
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # poor residuals

# try transformation -  better but still not great
m1<-lmer(log(1+Parasitized_aphids_hand_search) ~
           Year_since_est+
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # better, but still issues

#  data continuous but overdisperesed  - apply tweedie distribution

# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(
  Parasitized_aphids_hand_search ~  Treat +
    Year_since_est+ 
    Treat*Year_since_est +
    (1|Site/Year_since_est),                # Mean structure
  tweedie(link = "log"), # Tweedie family with log link
  data = data_field_mean_T1to3)
# Failed to converge with interactions  - simplify fixed effect structure

model_tweedie2 <- glmmTMB(
  Parasitized_aphids_hand_search ~  Treat +
    Year_since_est+ 
    (1|Site/Year_since_est),                # Mean structure
  tweedie(link = "log"), # Tweedie family with log link
  data = data_field_mean_T1to3)
# converged

model_tweedie3 <- glmmTMB(
  Parasitized_aphids_hand_search ~  Treat +
    (1|Site/Year_since_est),                # Mean structure
  tweedie(link = "log"), # Tweedie family with log link
  data = data_field_mean_T1to3)
# Likelihood ratio test
anova(model_tweedie2, model_tweedie3, test = "Chisq")
plot(simulationOutput) # looks good
summary(model_tweedie3)



model_tweedie4 <- glmmTMB(
  Parasitized_aphids_hand_search ~  
    Year_since_est+ 
    (1|Site/Year_since_est),                # Mean structure
  tweedie(link = "log"), # Tweedie family with log link
  data = data_field_mean_T1to3)
# Likelihood ratio test
anova(model_tweedie2, model_tweedie4, test = "Chisq")
plot(simulationOutput) # looks good
summary(model_tweedie4)


ggplot(data=data_field_mean_T1to3, aes(y=log(Parasitized_aphids_hand_search+1) , x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Aphids (ln N+1)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  



#-----------------------------------------------------------------------------------
# OSR pollination



# subletting for oilseed pollinated crops where the polliantion assessments were made
OSR <- data_field_mean_T1to3 %>% filter(Pollination_assessment_OSR == 1)




m_OSR<-lmer(Yield_from_poll ~                               
              Treat +
              
              (1|Site/Year_since_est), data=OSR,  na.action = "na.fail", REML = FALSE) 

drop1(m_OSR,test="F",sumFun=KRSumFun)
summary(m_OSR)
simulationOutput<-simulateResiduals(fittedModel=m_OSR, plot=F)
plot(simulationOutput)   # QQ plots are fine -  look like convergence for quarantines failed (low replication0), model ok


ggplot(data=OSR, aes(y=Yield_from_poll , x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Seed set (g)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  




#-----------------------------------------------------------------------------
# soil carbon 

soil<-read.csv("Soil_chemistry.csv")
# description of covariates
# Site	= site
# Treatment	= treatments 1 to 3 (1=BU, 2=Enhancing-ES, 3=Maximising-ES)
# C.organic_perc	=  percentage organic soil carbon - excluding inorganic carbon 
# Bulk_Density_g_cm3	= averge field soil bulk density g cm3
# Soil_C_orgnaic_gcm3 = average soil organic carbin content in g cm3


# basic model/ note no year effects as only assessed once in winter of 2021
m1<-lmer( Soil_C_orgnaic_gcm3~
            Treat +
            (1|Site), data=soil,  na.action = "na.fail", REML = FALSE) 



# final model
m1<-lmer( Soil_C_orgnaic_gcm3~
            Treat +
            (1|Site), data=soil,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # QQ plots are fine, model ok

ggplot(data=soil, aes(y=Soil_C_orgnaic_gcm3 , x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 0.13))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Soil Carbon (g.cm^3)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )   



###################################################################################
#  ES supporting taxa

# basic model
m1<-lmer( ~
            Treat +
            Year_since_est+
            Treat*Year_since_est +
            (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)

#  Responses are all continuous metrics as they are based on average sample point 
# catches (or similar reported values.).  They are as such not count data and so poisson /  neg bin
# not suitable. however, at points KS has proved issues so a log transformation of raw data applied
# to improve dHARMAS QQ plot checks of model assumptions


# staph_Pitfall_N_mean

m1<-lmer( log(staph_Pitfall_N_mean+1)~
            1+
            (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # QQ plots are fine, model ok

ggplot(data=data_field_mean_T1to3, aes(y=staph_Pitfall_N_mean , x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  #coord_cartesian(ylim = c(0, 0.13))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Abundance (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  


#Worms_Anecic_ave

m1<-lmer( log(Worms_Anecic_ave+1)~
            Treat +
            Year_since_est+
            Treat*Year_since_est +
            (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # QQ plots are fine, model ok


ggplot(data=data_field_mean_T1to3, aes(y=Worms_Anecic_ave , x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  #coord_cartesian(ylim = c(0, 0.13))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Abundance (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  



# Carabidae_Pitfall_N_mean_ln

m1<-lmer(log(1+Carabidae_Pitfall_N_mean) ~
           1+
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # QQ plots are fine, model ok

ggplot(data=data_field_mean_T1to3, aes(y=Carabidae_Pitfall_N_mean , x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 35
  ))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Abundance (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  


#spider_Pitfall_N_mean_ln

m1<-lmer(log(1+spider_Pitfall_N_mean) ~
           Treat +
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # QQ plots are fine, model ok


ggplot(data=data_field_mean_T1to3, aes(y=spider_Pitfall_N_mean , x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 35
                    ))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Abundance (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  




# Parasitoid_OSR_N
# Limit data set' to OSR crops only
# All parasitoids 
# subletting for oilseed  crops 
OSR_parasitoids <- data_field_mean_T1to3 %>% filter(Cereal == "Flower")


m1<-lmer( log(Parasitoid_OSR_N+1)~
            Year_since_est+
            (1|Site/Year_since_est), data=OSR_parasitoids ,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # QQ plots are fine, model ok


ggplot(data=OSR_parasitoids , aes(y=Parasitoid_OSR_N, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 75))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Abundance (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  

# Parasitoid_Cereal_N
# limit data set to cereal crops only
# subletting for oilseed  crops 
Cereal_parasitoids <- data_field_mean_T1to3 %>% filter(Cereal == "Cereal")

# Crop pest parasitoids of aphids, midges, beetles etc
#  Excludes hyperparasitoids.  Excludes parasitoids of umbellifer aphids adn non-crop aphids

m1<-lmer( log(Parasitoid_Cereal_N+1)~
            Year_since_est+
            (1|Site/Year_since_est), data=Cereal_parasitoids ,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # QQ plots are fine, model ok

ggplot(data=Cereal_parasitoids , aes(y=Parasitoid_Cereal_N, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  #coord_cartesian(ylim = c(0, 75))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Abundance (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  

#--------------------------------------------------------------------------------------------
# Predaors_hand_search_sum_ln

m1<-lmer( log(1+Predaors_hand_search_sum)~
            Treat +
            Year_since_est+
            (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)   # QQ plots suggest some evidence of problem

# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(Predaors_hand_search_sum
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #

model_tweedie2 <- glmmTMB(Predaors_hand_search_sum
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) # 

# Likelihood ratio test
anova(model_tweedie1, model_tweedie2, test = "Chisq")
# no sig interaction

model_tweedie3 <- glmmTMB(Predaors_hand_search_sum
                          ~  Treat +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) #  looks good

# Likelihood ratio test
anova(model_tweedie3, model_tweedie2, test = "Chisq")
# sig


model_tweedie4 <- glmmTMB(Predaors_hand_search_sum
                          ~  
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #

# Likelihood ratio test
anova(model_tweedie2, model_tweedie4, test = "Chisq")

ggplot(data=data_field_mean_T1to3 , aes(y=Predaors_hand_search_sum, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 20))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Abundance (N)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  

#---------------------------------------------------------------------------------------------
# Bee_N_margin_m2_ln


m1<-lmer(log(1+ Bee_N_margin_m2)
         ~
           Treat +
           Year_since_est+
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
emmeans(m1, ~ Treat*Year_since_est)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)  #   some issues

# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(Bee_N_margin_m2
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #

model_tweedie2 <- glmmTMB(Bee_N_margin_m2
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #   Looks much better -some issues possible with quantile deviations, but acceptable

# Likelihood ratio test
anova(model_tweedie1, model_tweedie2, test = "Chisq")
# significant interaction between year and management

ggplot(data=data_field_mean_T1to3 , aes(y=Bee_N_margin_m2, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 1.5))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Density"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  

#----------------------------------------------------------------------------

# Hoverfly_N_margin_m2_ln


m1<-lmer(log(1+Hoverfly_N_margin_m2)
         ~
           Treat +
           Year_since_est+
           Treat*Year_since_est +
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)  #   some issues
summary(m1)
emmeans(m1, ~ Treat*Year_since_est)

# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(Hoverfly_N_margin_m2
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #

model_tweedie2 <- glmmTMB(Hoverfly_N_margin_m2
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #   these look fine now

# Likelihood ratio test
anova(model_tweedie1, model_tweedie2, test = "Chisq")

ggplot(data=data_field_mean_T1to3 , aes(y=Hoverfly_N_margin_m2, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 1.5))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Density"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  

#----------------------------------------------------------------------------------------
# Parasitoid_N_margin_m2_ln
# oulier in EH_T1_2019 removed

MArgin_parasitoids_no_outlier_EH_T1_2019 <- data_field_mean_T1to3[!is.na(data_field_mean_T1to3$Parasitoid_N_margin_m2), ]   #  exclude missing value linked with outlier


m1<-lmer(log(1+Parasitoid_N_margin_m2)
         ~
           Treat +
           Year_since_est+
           
           (1|Site/Year_since_est), data=MArgin_parasitoids_no_outlier_EH_T1_2019,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)# some problems


# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(Parasitoid_N_margin_m2
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = MArgin_parasitoids_no_outlier_EH_T1_2019)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #

model_tweedie2 <- glmmTMB(Parasitoid_N_margin_m2
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = MArgin_parasitoids_no_outlier_EH_T1_2019)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) # 

# Likelihood ratio test
anova(model_tweedie1, model_tweedie2, test = "Chisq")
# ns

model_tweedie3 <- glmmTMB(Parasitoid_N_margin_m2
                          ~  Treat +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = MArgin_parasitoids_no_outlier_EH_T1_2019)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) # looks fine
summary(model_tweedie3 )

# Likelihood ratio test
anova(model_tweedie2, model_tweedie3, test = "Chisq")
#ns


model_tweedie4 <- glmmTMB(Parasitoid_N_margin_m2
                          ~  
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = MArgin_parasitoids_no_outlier_EH_T1_2019)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #  

# Likelihood ratio test
anova(model_tweedie2, model_tweedie4, test = "Chisq")
# sig

ggplot(data=data_field_mean_T1to3 , aes(y=Parasitoid_N_margin_m2, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 0.5))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Density"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  

########################################################################################
# pests

#blackgrass_ave	

data_field_mean_T1to3$Crop_type

Winter_Crop <- data_field_mean_T1to3 %>% filter(Crop_type == "Winter")

m1<-lmer(log(blackgrass_ave+1)	 ~
           1+ 
           (1|Site/Year_since_est), data=Winter_Crop,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) # some problems 

# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(blackgrass_ave
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = Winter_Crop)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #  much bimproved, slight evidence of oulier.  model OK

model_tweedie2 <- glmmTMB(blackgrass_ave
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = Winter_Crop)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) # model looks ok

# Likelihood ratio test
anova(model_tweedie1, model_tweedie2, test = "Chisq")


model_tweedie3 <- glmmTMB(blackgrass_ave
                          ~  Treat +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = Winter_Crop)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) #

# Likelihood ratio test
anova(model_tweedie2, model_tweedie3, test = "Chisq")



model_tweedie4 <- glmmTMB(blackgrass_ave
                          ~  
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = Winter_Crop)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #

# Likelihood ratio test
anova(model_tweedie2, model_tweedie4, test = "Chisq")


ggplot(data=Winter_Crop , aes(y=blackgrass_ave, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 5))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Tiller counts"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  


#-------------------------------------------------------------------------------------

#Weeds_ave
m1<-lmer(log(1+Weeds_ave) ~
           Year_since_est+ 
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) # some problems 

# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(Weeds_ave
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)

# model modle not converge  - needs fixed effects to be simplified


model_tweedie2 <- glmmTMB(Weeds_ave
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) # looks fine


model_tweedie3 <- glmmTMB(Weeds_ave
                          ~  Treat +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) #

# Likelihood ratio test
anova(model_tweedie2, model_tweedie3, test = "Chisq")



model_tweedie4 <- glmmTMB(Weeds_ave
                          ~  
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #

# Likelihood ratio test
anova(model_tweedie2, model_tweedie4, test = "Chisq")

ggplot(data=data_field_mean_T1to3 , aes(y=Weeds_ave, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
 # coord_cartesian(ylim = c(0, 5))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Plant count"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  

 
#------------------------------------------------------------------------------------
#Aphids_hand_search_ln

m1<-lmer( log(1+Aphids_hand_search)~
            Year_since_est+Treat +
            (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)  #  some problems


# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(Aphids_hand_search
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #  model looks fine

model_tweedie2 <- glmmTMB(Aphids_hand_search
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) # 

# Likelihood ratio test
anova(model_tweedie1, model_tweedie2, test = "Chisq")


model_tweedie3 <- glmmTMB(Aphids_hand_search
                          ~  Treat +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) # some over dispersion issues

# Likelihood ratio test
anova(model_tweedie2, model_tweedie3, test = "Chisq")
summary(model_tweedie3)


model_tweedie4 <- glmmTMB(Aphids_hand_search
                          ~  
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #  model is better but has problems

# Likelihood ratio test
anova(model_tweedie2, model_tweedie4, test = "Chisq")


ggplot(data=data_field_mean_T1to3 , aes(y=log(Aphids_hand_search), x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  # coord_cartesian(ylim = c(0, 5))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Abundance ln(N+1)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  


#----------------------------------------------------------------------------------------
#Slug_biomass
m1<-lmer( Slug_biomass~  1+
            (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)   
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)  #  problems

# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(Slug_biomass
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #  Much better -  but still some issues in particular with quantile deviations.  A big improvement though.

model_tweedie2 <- glmmTMB(Slug_biomass
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) # model looks ok

# Likelihood ratio test
anova(model_tweedie1, model_tweedie2, test = "Chisq")


model_tweedie3 <- glmmTMB(Slug_biomass
                          ~  Treat +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) #  model ok

# Likelihood ratio test
anova(model_tweedie2, model_tweedie3, test = "Chisq")



model_tweedie4 <- glmmTMB(Slug_biomass
                          ~  
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) # looks ok, some issues

# Likelihood ratio test
anova(model_tweedie2, model_tweedie4, test = "Chisq")

ggplot(data=data_field_mean_T1to3 , aes(y=Slug_biomass, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 5))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Mass (g)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  


#----------------------------------------------------------------------------------------
# Snail_biomass
m1<-lmer(Snail_biomass ~  1+
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) #  problems

# The Tweedie distribution is useful when your data are continuous and overdispersed, with support for both continuous (positive) and discrete components.
model_tweedie1 <- glmmTMB(Snail_biomass
                          ~  Treat +
                            Year_since_est+ 
                            Treat*Year_since_est +
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #  looks good

model_tweedie2 <- glmmTMB(Snail_biomass
                          ~  Treat +
                            Year_since_est+ 
                            (1|Site/Year_since_est),                # Mean structure
                          tweedie(link = "log"), # Tweedie family with log link
                          data = data_field_mean_T1to3)
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) # 

# Likelihood ratio test
anova(model_tweedie1, model_tweedie2, test = "Chisq")
summary(model_tweedie2)


ggplot(data=data_field_mean_T1to3 , aes(y=Snail_biomass, x=Treat, fill=Treat)) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey") +
  coord_cartesian(ylim = c(0, 0.2))  +# Adjusts view without removing data
  geom_boxplot(outlier.shape = NA) +  # Avoid duplicate outliers
  theme_bw()+   #  white background
  labs(
    x = NULL  ,                                     # X-axis legend
    y = "Mass (g)"                                        # Y-axis legend
  )  +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
        axis.text = element_text(size = 16),                      # Adjust axis text size
        axis.title = element_text(size = 16)                      # Adjust axis title size
  )  


#####################################################################################################
# Supplementary material additional analysis
#####################################################################################################
#####################################################################################################
#    Environmental predictors of key ES providing taxa

# example model showing tested fixed effects
#Year_since_est =  years since the treatments were implemented
#FM_sown_forb_SR_ave   =  the species richness of the sown plant component of the field margins  -  measuring success of establishment
#Drop_disk_ave  =  A measure of sward structure in the field 
#ratio_FM_to_crop	= The ratio of the area of field margins to the total area of the field.
#Ratio_IFS_to_crop  = the ratio of the area of in-field strips to the total area of the field.
#FYM  =  The application of farm yard manure in the previous year
#CoverCrop_inthepreceedingyear  =  Were cover crops grown in the previous year (technical over the winter period, but before monitoring)
# Shannon_landuse =  Shannon diversity of land use surrounding the experimetnal fields.


m1<-lmer(          ~ Year_since_est+
                     FM_sown_forb_SR_ave  +
                     Drop_disk_ave  +
                     ratio_FM_to_crop	+
                     Ratio_IFS_to_crop  +
                     FYM  +
                     CoverCrop_inthepreceedingyear +
                     Shannon_landuse+
                     (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)




#Worms_Anecic_ave

m1<-lmer(log(Worms_Anecic_ave+1)         ~
           Year_since_est+
           ratio_FM_to_crop	+
           FYM  +
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) #   Mostly ok  - some evidence of quantile devisations but acceptable



# Carabidae_Pitfall_N_mean_ln

m1<-lmer(log(Carabidae_Pitfall_N_mean+1) ~
           1 +
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) #  No problems


#staph_Pitfall_N_mean_ln

m1<-lmer(   log(1+staph_Pitfall_N_mean)
            ~
              1 +
              (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) # no problems


#spider_Pitfall_N_mean_ln

m1<-lmer(log(1+spider_Pitfall_N_mean) ~
           FM_sown_forb_SR_ave  +
           (1|Site/Year_since_est), data=data_field_mean_T1to3,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) #  No problems



#----------------------------------------------------------------------------------------
# Predators_hand_search_sum_ln

#  Tweedie distributed as data is overdispersed

data_field_mean_T1to3$metric<-data_field_mean_T1to3$Predaors_hand_search_sum

m1int<-glmmTMB(   metric       ~ 1+
                    (1|Site/Year_since_est),   
                  tweedie(link = "log"), # Tweedie family with log link
                  data=data_field_mean_T1to3 ) 

model_tweedie1<-glmmTMB(   metric       ~ Year_since_est+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie1, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie1, plot=F)
plot(simulationOutput) #

model_tweedie1<-glmmTMB(   metric       ~    Drop_disk_ave  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie1, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie1, plot=F)
plot(simulationOutput) #

model_tweedie2<-glmmTMB(   metric       ~     ratio_FM_to_crop	+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie2, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #
summary(model_tweedie2)

model_tweedie3<-glmmTMB(   metric       ~  Ratio_IFS_to_crop  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie3, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) #
summary(model_tweedie3)


model_tweedie4<-glmmTMB(   metric       ~   FYM  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie4, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #


model_tweedie5<-glmmTMB(   metric       ~ 
                             CoverCrop_inthepreceedingyear +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie5, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie5, plot=F)
plot(simulationOutput) #

model_tweedie6<-glmmTMB(   metric       ~ 
                             Shannon_landuse+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie6, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie6, plot=F)
plot(simulationOutput) #

model_tweedie7<-glmmTMB(   metric       ~ 
                             FM_sown_forb_SR_ave+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie7, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie7, plot=F)
plot(simulationOutput) #
summary(model_tweedie7)

#----------------------------------------------------------------------------------------
# Cereal_parasitoids

m1<-lmer(   log(1+Parasitoid_Cereal_N)     ~ 
              CoverCrop_inthepreceedingyear +
              (1|Site/Year_since_est), data=Cereal_parasitoids,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) #  no problems

#----------------------------------------------------------------------------------------
# OSR_parasitoids

m1<-lmer(   log(1+Parasitoid_OSR_N)     ~ Year_since_est+
              (1|Site/Year_since_est), data=OSR_parasitoids,  na.action = "na.fail", REML = FALSE) 
drop1(m1,test="F",sumFun=KRSumFun)
summary(m1)
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput) #  No problems

#---------------------------------------------------------------------------------------------
# margins
# Bee_N_margin_m2_ln	

#  Tweedie distributed as data is overdispersed

data_field_mean_T1to3$metric<-data_field_mean_T1to3$Bee_N_margin_m2

m1int<-glmmTMB(   metric       ~ 1+
                    (1|Site/Year_since_est),   
                  tweedie(link = "log"), # Tweedie family with log link
                  data=data_field_mean_T1to3 ) 

model_tweedie1<-glmmTMB(   metric       ~ Year_since_est+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie1, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie1, plot=F)
plot(simulationOutput) #

model_tweedie1<-glmmTMB(   metric       ~    Drop_disk_ave  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie1, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie1, plot=F)
plot(simulationOutput) #
summary(model_tweedie1)

model_tweedie2<-glmmTMB(   metric       ~     ratio_FM_to_crop	+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie2, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #
summary(model_tweedie2)


model_tweedie3<-glmmTMB(   metric       ~  Ratio_IFS_to_crop  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie3, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) #
summary(model_tweedie3)

model_tweedie4<-glmmTMB(   metric       ~   FYM  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie4, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #


model_tweedie5<-glmmTMB(   metric       ~ 
                             CoverCrop_inthepreceedingyear +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie5, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie5, plot=F)
plot(simulationOutput) #

model_tweedie6<-glmmTMB(   metric       ~ 
                             Shannon_landuse+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie6, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie6, plot=F)
plot(simulationOutput) #

model_tweedie7<-glmmTMB(   metric       ~ 
                             FM_sown_forb_SR_ave+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie7, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie7, plot=F)
plot(simulationOutput) #
summary(model_tweedie7)


#--------------------------------------------------------------------------------
# Hoverfly_N_margin_m2_ln

#  Tweedie distributed as data is overdispersed

data_field_mean_T1to3$metric<-data_field_mean_T1to3$Hoverfly_N_margin_m2

m1int<-glmmTMB(   metric       ~ 1+
                    (1|Site/Year_since_est),   
                  tweedie(link = "log"), # Tweedie family with log link
                  data=data_field_mean_T1to3 ) 

model_tweedie1<-glmmTMB(   metric       ~ Year_since_est+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie1, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie1, plot=F)
plot(simulationOutput) #

model_tweedie1<-glmmTMB(   metric       ~    Drop_disk_ave  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie1, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie1, plot=F)
plot(simulationOutput) #
summary(model_tweedie1)

model_tweedie2<-glmmTMB(   metric       ~     ratio_FM_to_crop	+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie2, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #
summary(model_tweedie2)

model_tweedie3<-glmmTMB(   metric       ~  Ratio_IFS_to_crop  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie3, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) #
summary(model_tweedie3)

model_tweedie4<-glmmTMB(   metric       ~   FYM  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie4, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #
summary(model_tweedie4)


model_tweedie5<-glmmTMB(   metric       ~ 
                             CoverCrop_inthepreceedingyear +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie5, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie5, plot=F)
plot(simulationOutput) #
summary(model_tweedie5)

model_tweedie6<-glmmTMB(   metric       ~ 
                             Shannon_landuse+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie6, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie6, plot=F)
plot(simulationOutput) #

model_tweedie7<-glmmTMB(   metric       ~ 
                             FM_sown_forb_SR_ave+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=data_field_mean_T1to3 ) 
anova(m1int, model_tweedie7, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie7, plot=F)
plot(simulationOutput) #



#------------------------------------------------------------------------------------
#Parasitoid_N_margin_m2_ln	
# oulier in EH_T1_2019 removed

MArgin_parasitoids_no_outlier_EH_T1_2019 <- data_field_mean_T1to3[!is.na(data_field_mean_T1to3$Parasitoid_N_margin_m2), ]   #  exclude missing value linked with outlier

MArgin_parasitoids_no_outlier_EH_T1_2019$metric<-MArgin_parasitoids_no_outlier_EH_T1_2019$Parasitoid_N_margin_m2

m1int<-glmmTMB(   metric       ~ 1+
                    (1|Site/Year_since_est),   
                  tweedie(link = "log"), # Tweedie family with log link
                  data=MArgin_parasitoids_no_outlier_EH_T1_2019 ) 

model_tweedie1<-glmmTMB(   metric       ~ Year_since_est+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=MArgin_parasitoids_no_outlier_EH_T1_2019) 
anova(m1int, model_tweedie1, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie1, plot=F)
plot(simulationOutput) #

model_tweedie1<-glmmTMB(   metric       ~    Drop_disk_ave  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=MArgin_parasitoids_no_outlier_EH_T1_2019 ) 
anova(m1int, model_tweedie1, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie1, plot=F)
plot(simulationOutput) #
summary(model_tweedie1)

model_tweedie2<-glmmTMB(   metric       ~     ratio_FM_to_crop	+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=MArgin_parasitoids_no_outlier_EH_T1_2019 ) 
anova(m1int, model_tweedie2, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie2, plot=F)
plot(simulationOutput) #
summary(model_tweedie2)

model_tweedie3<-glmmTMB(   metric       ~  Ratio_IFS_to_crop  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=MArgin_parasitoids_no_outlier_EH_T1_2019 ) 
anova(m1int, model_tweedie3, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie3, plot=F)
plot(simulationOutput) #
summary(model_tweedie3)

model_tweedie4<-glmmTMB(   metric       ~   FYM  +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=MArgin_parasitoids_no_outlier_EH_T1_2019 ) 
anova(m1int, model_tweedie4, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie4, plot=F)
plot(simulationOutput) #
summary(model_tweedie4)


model_tweedie5<-glmmTMB(   metric       ~ 
                             CoverCrop_inthepreceedingyear +
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=MArgin_parasitoids_no_outlier_EH_T1_2019) 
anova(m1int, model_tweedie5, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie5, plot=F)
plot(simulationOutput) #
summary(model_tweedie5)

model_tweedie6<-glmmTMB(   metric       ~ 
                             Shannon_landuse+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=MArgin_parasitoids_no_outlier_EH_T1_2019) 
anova(m1int, model_tweedie6, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie6, plot=F)
plot(simulationOutput) #

model_tweedie7<-glmmTMB(   metric       ~ 
                             FM_sown_forb_SR_ave+
                             (1|Site/Year_since_est),   
                           tweedie(link = "log"), # Tweedie family with log link
                           data=MArgin_parasitoids_no_outlier_EH_T1_2019) 
anova(m1int, model_tweedie7, test = "Chisq")
simulationOutput<-simulateResiduals(fittedModel=model_tweedie7, plot=F)
plot(simulationOutput) #




########################################################################################
#  Green infastrucure effects on yield and profit


#  Yield_SMD

yield_quadrat_SMD<-lm(SMD_Yield_Tonnes_ha_corrected     ~
                        ratio_cGI_to_crop  , data=data_field_mean_T1to3) 
drop1(yield_quadrat_SMD,test="F",sumFun=KRSumFun)

summary(yield_quadrat_SMD)

ggplot(data_field_mean_T1to3, aes(x = ratio_cGI_to_crop, y = SMD_Yield_Tonnes_ha_corrected , color=Treat_name)) +
  geom_point(size=4) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +  # Regression line with ribbon
  labs( x = "Ratio by area of GI:crop", y = "Yield (SMD)", color = "System") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 18, face = "bold"),  # Change axis titles font size
    axis.text = element_text(size = 14),   # Change axis labels font size
    legend.title = element_text(size = 14),  # Change legend title font size
    legend.text = element_text(size = 14)    # Change legend text font size
  )




#  Pofit_no subsidies

profit_SMD<-lm(SMD_Ave_value_crop_ha_noAESpayment    ~
                 ratio_cGI_to_crop   , data=data_field_mean_T1to3) 
drop1(profit_SMD,test="F",sumFun=KRSumFun)

summary(profit_SMD)

# model checks
simulationOutput<-simulateResiduals(fittedModel=profit_SMD, plot=F)
plot(simulationOutput)

ggplot(data_field_mean_T1to3, aes(x = ratio_cGI_to_crop, y = SMD_Ave_value_crop_ha_noAESpayment, color=Treat_name)) +
  geom_point(size=4) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +  # Regression line with ribbon
  labs( x = "Ratio by area of GI:crop", y = "Profit (SMD)", color = "System") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 18, face = "bold"),  # Change axis titles font size
    axis.text = element_text(size = 14),   # Change axis labels font size
    legend.title = element_text(size = 14),  # Change legend title font size
    legend.text = element_text(size = 14)    # Change legend text font size
  )




#----------------------------------------------------------------
#  Profit with subsidies

m1<-lm(SMD_Ave_value_crop_ha_withAESpayment    ~
         ratio_cGI_to_crop   , data=data_field_mean_T1to3) 
drop1(profit_SMD,test="F",sumFun=KRSumFun)


summary(m1)

# model checks
simulationOutput<-simulateResiduals(fittedModel=m1, plot=F)
plot(simulationOutput)

ggplot(data_field_mean_T1to3, aes(x = ratio_cGI_to_crop, y = SMD_Ave_value_crop_ha_withAESpayment , color=Treat_name)) +
  geom_point(size=4) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +  # Regression line with ribbon
  labs( x = "Ratio by area of GI:crop", y = "Profit (SMD) + AES subsidies", color = "System") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 18, face = "bold"),  # Change axis titles font size
    axis.text = element_text(size = 14),   # Change axis labels font size
    legend.title = element_text(size = 14),  # Change legend title font size
    legend.text = element_text(size = 14)    # Change legend text font size
  )



########################################################
# field sizes


field_size <- read.csv("field_size.csv")

# description of variables in field_size.csv
# Site	= site descriptor
# Treat	= treatment number, where 1=BAU control, 2= Enhancing-ES, 3= Maximising-ES
# Field_total_Area_ha = field total area in ha


m1<-lm(Field_total_Area_ha ~  Treat
       
       , data=field_size)
drop1(m1, test="F",sumFun=KRSumFun)
# No significant difference in field size between treatments
summary(m1)
