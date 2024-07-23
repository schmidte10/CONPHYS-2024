---
title: "Absolute scope - Oxygen Consumption"
author: "Elliott Schmidt"
date: "06 December, 2023"
output:
  html_document:
    keep_md: yes
    code_folding: show
    collapse: no
    df_print: paged
    fig_caption: yes
    fig_height: 4
    fig_width: 6
    highlight: monochrome
    theme: flatly
    latex_engine: xelatex
    toc: yes
    toc_float: yes
    css: styles.css
  pdf_document:
    df_print: default
    fig_caption: yes
    fig_height: 4
    fig_width: 4
    highlight: tango
    latex_engine: xelatex
    number_sections: yes
    toc_depth: 2
documentclass: article
fontsize: 12pt
mainfont: Arial
mathfont: LiberationMono
classoption: a4paper
--- 

# Scenario 

For initial details on the experiment performed please read the **ReadMe** file. In breif, _Acanthochromis polyacanthus_ from two different regions on the Great Barrier Reef (GBR) were tested for metabolic performance at four different temperatures, 27$^\circ$C, 28.5$^\circ$C, 30$^\circ$C, and 31.5$^\circ$C. Fish used in this study were collected from two different regions, low- (i.e. Cairns) and high-latitude (i.e., Mackay), within each region fish were collected from a total of three different populations. Individuals were tested at each temperature, resting oxygen consumption, maximum oxygen consumption. Absolute aerboic scope was calculated by using the following formula: 

Absolute aerobic scope = (maximum oxygen consumption - resting oxygen consumption)

Individuals were first tested at 27$^\circ$C. Water temperature was then increased at a rate of 0.5$^\circ$C Day^-1 until the next temperature was reached. Fish were then provided with an additional 5 day to adjust to the new temperature before aerobic physiology was tested again. 

Three traits are included within the aerobic physiology analysis, resting oxygen consumption, maximum oxygen consumption, and absoulte aerboic scope. Data for each metric was collect from respiratory experiments that had data recorded via a combination of programs including, AquaResp and PyroScience. Slopes (i.e., resting and maximum oxygen consumption values) were then calculated via the **RespR** [https://januarharianto.github.io/respR/articles/respR.html] package.  


# Read in the data

Before beginning always make sure that you are working in the correct directory 




```r
knitr::opts_knit$set(root.dir=working.dir)
```

Lets start by loading the packages that are needed 

## Load packages 


```r
library(tidyverse) # data manipulation
library(plyr) # data manipulation
library(dplyr) # data manipulation
library(lubridate) # data manipulation - specifically time data
library(ggplot2) # plotting figures
library(glmmTMB) # running models
library(performance) # model validation
library(car) # dependent
library(DHARMa) # model validation
library(MuMIn) # model validation
library(kableExtra) # creating tables
library(broom) # dependent
library(emmeans) # post-hoc analysis
library(ggeffects) # plotting models/model validation
library(vtable) # creating tables
library(modelr) # model validation
library(kableExtra) # formatting output tables
library(sjPlot) # plotting models
```



Now we can import that data. Replace import data with the PATH to your data file. I have secretly labelled my PATH import.data (i.e. import.data = "PATH TO MY FILE")

## Load data 



```r
resp <- import.data
```

# Data manipulation 

Before the data can be analysed it is important to clean up the data file. Below a number of adjustments are made, primarily making sure that columns are being treated appropriately as either factors, numeric, or as time, as well as the renaming of some columns. Once these changes are made the data is being saved into a new dataframe called **resp2** 


```{.r .style}
resp2 = resp %>% 
  dplyr::rename(EXP_FISH_ID = FISH_ID) %>%
  separate(EXP_FISH_ID, c("FISH_ID"), remove = FALSE) %>%
  mutate(FISH_ID = factor(FISH_ID), 
         POPULATION = factor(POPULATION), 
         REGION = factor(REGION), 
         TEMPERATURE = as.numeric(TEMPERATURE), #run with temperature as a factor
         RESTING_DATE = factor(RESTING_DATE), 
         RESTING_CHAMBER = factor(RESTING_CHAMBER), 
         RESTING_SYSTEM = factor(RESTING_SYSTEM), 
         RESTING_SUMP = factor(RESTING_SUMP), 
         RESTING_AM_PM = factor(RESTING_AM_PM), 
         RESTING_START_TIME = hms(RESTING_START_TIME),
         RESTING_END_TIME = hms(RESTING_ENDTIME),
         MAX_DATE = factor(MAX_DATE), 
         MAX_CHAMBER = factor(MAX_CHAMBER), 
         MAX_SYSTEM = factor(MAX_SYSTEM), 
         MAX_SUMP = factor(MAX_SUMP), 
         MAX_AM_PM = factor(MAX_AM_PM), 
         MAX_START_TIME = hms(MAX_START_TIME), 
         Swim.performance = factor(Swim.performance), 
         NAS = as.numeric(NAS), 
         FAS = as.numeric(FAS), 
         MgO2.hr_Net = as.numeric(MgO2.hr_Net), 
         RESTING_RUNTIME_SECONDS = as.numeric(hms(RESTING_RUNTIME))) %>% 
  dplyr::rename(MASS = WEIGHT) %>% 
  mutate(MASS_CENTERED = scale(MASS, scale = FALSE, center = TRUE))
```

Next select data points will be removed. Beside the removed data points I have provided reasoning for their exclusion, such as fish died during the experiment, or data quailty was poor - which likely indicated that there was an issue with the equipment during the trial. 


```r
resp3 <- resp2 %>% 
  subset(  
    EXP_FISH_ID !="LCHA127_27" & # deceased during experiment
      EXP_FISH_ID !="LCHA132_27" & # deceased during experiment
      EXP_FISH_ID !="LKES168_27" # poor data quality
  ) 
```

So far the analysis has been the same as the protocol outlined in the **aerobic physiology resting** data. One additional data removal step will take place in the maximum oxygen consumption analysis for samples where fish swam poorly and therefore their maximum oxygen consumption data is thought to be unreliable. This step is done before any data analysis has taken place. 


```r
resp4 <- resp3 %>% 
  subset(
    EXP_FISH_ID !="CSUD008_27" &  # poor swim
      EXP_FISH_ID !="CSUD008_30" &  # poor swim 
      EXP_FISH_ID !="CSUD008_28.5" & # poor swim
      EXP_FISH_ID !="CSUD018_31.5" & # poor swim 
      EXP_FISH_ID !="CSUD026_30" & # max. value low 
      EXP_FISH_ID !="CSUD074_28.5" & # fas value low 
      EXP_FISH_ID !="CSUD079_30" &
      EXP_FISH_ID !="CVLA052_27" & #nas value low 
      EXP_FISH_ID !="CVLA054_28.5" & # low max value? 
      EXP_FISH_ID !="LCHA113_27" & # poor data quality 
      EXP_FISH_ID !="LCHA113_30" & # poor swim 
      EXP_FISH_ID !="LCHA127_27" # deceased during experiment
  ) 
```

# Exploratory data analysis {.tabset}

## Mass v Rest


```r
ggplot(resp4, aes(MASS, MgO2.hr_Net)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_classic()
```

![](aerobic_physiology_absolute_scope_files/figure-html/eda-1-1.png)<!-- -->

## Mass v REST (LATITUDE)

```r
ggplot(resp4, aes(MASS, MgO2.hr_Net, color = REGION)) + 
  geom_point() +
  theme_classic() + 
  geom_smooth(method = "lm")
```

![](aerobic_physiology_absolute_scope_files/figure-html/eda-2-1.png)<!-- -->

## TEMPERTURE v REST (LATITUDE)

```r
ggplot(resp4, aes(TEMPERATURE, MgO2.hr_Net, color = REGION)) + 
  geom_point() +
  theme_classic()
```

![](aerobic_physiology_absolute_scope_files/figure-html/eda-3-1.png)<!-- -->

## {-}

# Fit the model 

The model was fit using the **glm** and later **glmmTMB** package in R. A number of different models were tested to determine which hypothesis and associated variables best predicted resting oxygen consumption. Model fit was examined using AICc, BIC, and r-squared values. Additional model were examined via the validation diagonistics provided by the **performance** and **dHARMA** packages in R. 

The first set of models tested looked at three different hypotheses including 1) that mass has a major impact of resting oxygen consumption of fish (this has been documented in the literature), 2) if variables related to time have an important impact on the resting oxygen consumption of fish. 

## Fixed factors (linear regression models)

### model 1

```r
#--- base model ---#
nas.1 <- glm(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4)  
```
#### summary
<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
   <th style="text-align:right;"> t value </th>
   <th style="text-align:right;"> Pr(&gt;|t|) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -4.6505476 </td>
   <td style="text-align:right;"> 4.5988210 </td>
   <td style="text-align:right;"> -1.011248 </td>
   <td style="text-align:right;"> 0.3133266 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> 17.8617376 </td>
   <td style="text-align:right;"> 6.8020030 </td>
   <td style="text-align:right;"> 2.625953 </td>
   <td style="text-align:right;"> 0.0094246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TEMPERATURE </td>
   <td style="text-align:right;"> 0.4921098 </td>
   <td style="text-align:right;"> 0.1581257 </td>
   <td style="text-align:right;"> 3.112143 </td>
   <td style="text-align:right;"> 0.0021773 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 167.6469322 </td>
   <td style="text-align:right;"> 24.9382811 </td>
   <td style="text-align:right;"> 6.722473 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:TEMPERATURE </td>
   <td style="text-align:right;"> -0.6707332 </td>
   <td style="text-align:right;"> 0.2336280 </td>
   <td style="text-align:right;"> -2.870945 </td>
   <td style="text-align:right;"> 0.0046100 </td>
  </tr>
</tbody>
</table>

### model 2

```r
#--- experimental rmr equipment hypothesis ---#
nas.2 <- glm(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_SUMP + RESTING_RUNTIME_SECONDS + 
                   RESTING_AM_PM, 
                 family=gaussian(),
                 data = resp4) 
```

#### summary
<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
   <th style="text-align:right;"> t value </th>
   <th style="text-align:right;"> Pr(&gt;|t|) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -4.1166901 </td>
   <td style="text-align:right;"> 4.6900539 </td>
   <td style="text-align:right;"> -0.8777490 </td>
   <td style="text-align:right;"> 0.3813335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> 17.7878857 </td>
   <td style="text-align:right;"> 6.8578591 </td>
   <td style="text-align:right;"> 2.5937958 </td>
   <td style="text-align:right;"> 0.0103304 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TEMPERATURE </td>
   <td style="text-align:right;"> 0.4314945 </td>
   <td style="text-align:right;"> 0.1732534 </td>
   <td style="text-align:right;"> 2.4905397 </td>
   <td style="text-align:right;"> 0.0137265 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 167.2794598 </td>
   <td style="text-align:right;"> 25.2897560 </td>
   <td style="text-align:right;"> 6.6145146 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RESTING_SUMP2 </td>
   <td style="text-align:right;"> 0.0078043 </td>
   <td style="text-align:right;"> 0.4108966 </td>
   <td style="text-align:right;"> 0.0189932 </td>
   <td style="text-align:right;"> 0.9848690 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RESTING_RUNTIME_SECONDS </td>
   <td style="text-align:right;"> 0.0000818 </td>
   <td style="text-align:right;"> 0.0000936 </td>
   <td style="text-align:right;"> 0.8735816 </td>
   <td style="text-align:right;"> 0.3835932 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RESTING_AM_PMPM </td>
   <td style="text-align:right;"> -0.1069842 </td>
   <td style="text-align:right;"> 0.3960028 </td>
   <td style="text-align:right;"> -0.2701603 </td>
   <td style="text-align:right;"> 0.7873685 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:TEMPERATURE </td>
   <td style="text-align:right;"> -0.6678804 </td>
   <td style="text-align:right;"> 0.2355090 </td>
   <td style="text-align:right;"> -2.8359020 </td>
   <td style="text-align:right;"> 0.0051315 </td>
  </tr>
</tbody>
</table>

### model comparison table
<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> model </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> AICc </th>
   <th style="text-align:right;"> BIC </th>
   <th style="text-align:right;"> r2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> nas.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 839.8550 </td>
   <td style="text-align:right;"> 858.3808 </td>
   <td style="text-align:right;"> 0.440634 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nas.2 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 845.5166 </td>
   <td style="text-align:right;"> 872.9666 </td>
   <td style="text-align:right;"> 0.439167 </td>
  </tr>
</tbody>
</table>

The model that contains **MASS_CENTERED** seems to do better than the model that incorporates variables that are associated with the time that experiments are performed.This is demonstrated by the lower AIC and BIC scores, as well as higher r-squared value.  

## Polynomials 

### polynomial models 

Note that the linear model has already been created via model _rmr.3_ in the previous section.


```r
#--- second order polynomial ---# 
nas.1.p2 <- glm(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE,2) + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4) 

#--- third order polynomial ---#
nas.1.p3 <- glm(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4) 
```

#### polynomial model comparisons
<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> model </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> AICc </th>
   <th style="text-align:right;"> BIC </th>
   <th style="text-align:right;"> r2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> nas.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 839.8550 </td>
   <td style="text-align:right;"> 858.3808 </td>
   <td style="text-align:right;"> 0.4463407 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nas.1.p2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 840.4431 </td>
   <td style="text-align:right;"> 864.9447 </td>
   <td style="text-align:right;"> 0.4580961 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nas.1.p3 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 844.9023 </td>
   <td style="text-align:right;"> 875.2738 </td>
   <td style="text-align:right;"> 0.4581327 </td>
  </tr>
</tbody>
</table>

From our model comparison we can see the there is no additional benefit to the model by including temperature as a 2^nd^ or 3^rd^ order polynomial. However, the linear and quadratic model both perform well. 

## Random factors 

Fish were repeatedly sampled over four different temperatures, therefore repeated sampling needs to be accounted for. To do this random factors will be included within the model. There are a number of options that can be used for random factors including 1) accounting for repeated sampling of individuals, 2) accounting for repeated sampling of individuals nested within population, 3) account for repeated sampling of individuals and populations without nesting. All three models will be run a compaired. 

### random factor models


```r
nas.1a <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE) 

nas.1b <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|POPULATION/FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

nas.1c <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID) + (REGION|POPULATION), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE) # convergence problem
```

#### random factor model comparisons 

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> model </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> AICc </th>
   <th style="text-align:right;"> BIC </th>
   <th style="text-align:right;"> r2m </th>
   <th style="text-align:right;"> r2c </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> nas.1a </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 828.7221 </td>
   <td style="text-align:right;"> 850.2488 </td>
   <td style="text-align:right;"> 0.4525497 </td>
   <td style="text-align:right;"> 0.5901940 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nas.1b </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 830.8106 </td>
   <td style="text-align:right;"> 855.3122 </td>
   <td style="text-align:right;"> 0.4542611 </td>
   <td style="text-align:right;"> 0.5928177 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nas.1c </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.4560622 </td>
   <td style="text-align:right;"> 0.5944637 </td>
  </tr>
</tbody>
</table>

Model _nas.3a_ appears to be the best model, however, there seems to be little difference in how the models change depending on how the random factors are arranged.

# Model validation {.tabset .tabset-faded}

## performance {.tabset .tabset-faded}

### rmr.3a (linear)
![](aerobic_physiology_absolute_scope_files/figure-html/model-valid-1-1.png)<!-- -->

The _nas.1a_ model performs well, however, in the model validation performed by the **performance** model it looks like there are two variables that are highly correlated. If we expand the figure we can see that the highly correlated variables are REGION and REGION:TEMPERATURE. Perhaps this is unsurprising  but lets see what happens when we run the quadratic (2^nd^ polynomial) model to see if this helps deal with the high correlation between these two variables, as it performed very similarly to _nas.1a_, and even had a higher r2 value. 

### nas.1.p2a (quadratic)

First we need to update the model by adding in the missing random factor

```r
nas.1.p2a <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE,2) + MASS_CENTERED + (1|FISH_ID), 
                 family=gaussian(),
                 data = resp4, 
                 REML=TRUE) 
```

![](aerobic_physiology_absolute_scope_files/figure-html/model-valid-1.2-1.png)<!-- -->

## DHARMa residuals {.tabset .tabset-faded}

### nas.1a (linear)

```r
nas.1a %>% simulateResiduals(plot=TRUE)
```

![](aerobic_physiology_absolute_scope_files/figure-html/model-valid-2-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.212 0.244 0.556 0.532 0.74 0.728 1 0.712 0.596 0.868 0.788 0.364 0.28 0.416 0.744 0.124 0.676 0.996 0.788 0.864 ...
```

```r
nas.1a %>% DHARMa::testResiduals(plot=TRUE)
```

![](aerobic_physiology_absolute_scope_files/figure-html/model-valid-2-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.063636, p-value = 0.4741
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.96215, p-value = 0.744
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 176, p-value = 0.4095
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001379164 0.040444565
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01136364
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.063636, p-value = 0.4741
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.96215, p-value = 0.744
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 176, p-value = 0.4095
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001379164 0.040444565
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01136364
```

### nas.1.p2 (quadratic)

First we need to update the model by adding in the missing random factor

```r
nas.1.p2a <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE,2) + MASS_CENTERED + (1|FISH_ID), 
                 family=gaussian(),
                 data = resp4, 
                 REML=TRUE)  
```


```r
nas.1.p2a %>% simulateResiduals(plot=TRUE) 
```

![](aerobic_physiology_absolute_scope_files/figure-html/model-valid-2.2-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.268 0.2 0.648 0.472 0.664 0.8 1 0.756 0.488 0.812 0.844 0.4 0.224 0.352 0.792 0.172 0.596 0.996 0.876 0.828 ...
```

```r
nas.1.p2a %>% DHARMa::testResiduals(plot=TRUE)
```

![](aerobic_physiology_absolute_scope_files/figure-html/model-valid-2.2-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.060045, p-value = 0.5497
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.95158, p-value = 0.64
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 176, p-value = 0.4095
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001379164 0.040444565
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01136364
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.060045, p-value = 0.5497
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.95158, p-value = 0.64
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 176, p-value = 0.4095
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001379164 0.040444565
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01136364
```

## {-}

# {-}

It looks like the model that treats temperature as a second order polynomial does a better job at avoiding high levels of collinearity within the model. The quadratic model will be used moving forward because it: 

* The **quadratic model** performs just as well as the linear model based on the model validation scores (i.e., AIC, BIC, and r2) 
* The **quadratic model** does a **better** job at dealing with collinearity that appeared in the model 

# Partial plots {.tabset .tabset-faded}

## ggemmeans 

![](aerobic_physiology_absolute_scope_files/figure-html/partial-plots-1-1.png)<!-- -->

## plot_model 

![](aerobic_physiology_absolute_scope_files/figure-html/partial-plots-2-1.png)<!-- -->

# {-} 

# Model investigation {.tabset .tabset-faded}

## summary 
<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> StdError </th>
   <th style="text-align:right;"> Zvalue </th>
   <th style="text-align:right;"> Pvalue </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> 9.676637 </td>
   <td style="text-align:right;"> 0.3717003 </td>
   <td style="text-align:right;"> 26.033437 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -1.632353 </td>
   <td style="text-align:right;"> 0.6124654 </td>
   <td style="text-align:right;"> -2.665216 </td>
   <td style="text-align:right;"> 0.0076939 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> 10.324322 </td>
   <td style="text-align:right;"> 3.0798654 </td>
   <td style="text-align:right;"> 3.352199 </td>
   <td style="text-align:right;"> 0.0008017 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -6.546341 </td>
   <td style="text-align:right;"> 3.0573845 </td>
   <td style="text-align:right;"> -2.141157 </td>
   <td style="text-align:right;"> 0.0322613 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 179.593117 </td>
   <td style="text-align:right;"> 32.6469695 </td>
   <td style="text-align:right;"> 5.501066 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> -14.500327 </td>
   <td style="text-align:right;"> 4.5325854 </td>
   <td style="text-align:right;"> -3.199129 </td>
   <td style="text-align:right;"> 0.0013784 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> 4.897260 </td>
   <td style="text-align:right;"> 4.4894335 </td>
   <td style="text-align:right;"> 1.090841 </td>
   <td style="text-align:right;"> 0.2753427 </td>
  </tr>
</tbody>
</table>

## Anova 
<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Chisq </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> Pr(&gt;Chisq) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> REGION </td>
   <td style="text-align:right;"> 6.537818 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0105605 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2) </td>
   <td style="text-align:right;"> 6.515869 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0384678 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 30.261721 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGION:poly(TEMPERATURE, 2) </td>
   <td style="text-align:right;"> 11.460373 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0032465 </td>
  </tr>
</tbody>
</table>

## confint 
<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> 2.5 % </th>
   <th style="text-align:right;"> 97.5 % </th>
   <th style="text-align:right;"> Estimate </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> 8.9481173 </td>
   <td style="text-align:right;"> 10.4051557 </td>
   <td style="text-align:right;"> 9.676637 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -2.8327627 </td>
   <td style="text-align:right;"> -0.4319423 </td>
   <td style="text-align:right;"> -1.632353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> 4.2878974 </td>
   <td style="text-align:right;"> 16.3607477 </td>
   <td style="text-align:right;"> 10.324322 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -12.5387048 </td>
   <td style="text-align:right;"> -0.5539778 </td>
   <td style="text-align:right;"> -6.546341 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 115.6062324 </td>
   <td style="text-align:right;"> 243.5800011 </td>
   <td style="text-align:right;"> 179.593117 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> -23.3840309 </td>
   <td style="text-align:right;"> -5.6166225 </td>
   <td style="text-align:right;"> -14.500327 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -3.9018683 </td>
   <td style="text-align:right;"> 13.6963875 </td>
   <td style="text-align:right;"> 4.897260 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Std.Dev.(Intercept)|FISH_ID </td>
   <td style="text-align:right;"> 0.8790124 </td>
   <td style="text-align:right;"> 1.9917277 </td>
   <td style="text-align:right;"> 1.323160 </td>
  </tr>
</tbody>
</table>

## r-squared
<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:right;"> R2_conditional </th>
   <th style="text-align:right;"> R2_marginal </th>
   <th style="text-align:left;"> optional </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.6031021 </td>
   <td style="text-align:right;"> 0.4622269 </td>
   <td style="text-align:left;"> FALSE </td>
  </tr>
</tbody>
</table>

# {-} 

# Pairwise comparisons {.tabset .tabset-faded} 

## emtrends [latitudes]



```r
nas.1.p2a %>% emtrends(var = "TEMPERATURE", type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"(Core MASS_CENTERED-0.000127610645933016) - (Leading MASS_CENTERED-0.000127610645933016)","2":"29.09659","3":"0.690616","4":"0.2070788","5":"174","6":"0.2819063","7":"1.099326","8":"3.335039","9":"0.00104241","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
SCROLL TO THE RIGHT -->

The numbers in the left most column in the table just mention that the slopes are assuming mean **MASS_CENTERED** values when looking at differences between latitudinal slopes.

## emmeans [latitudes]

```r
nas.1.p2a %>% emmeans(pairwise ~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"Core - Leading","2":"29.09659","3":"2.090783","4":"0.7384902","5":"174","6":"0.633231","7":"3.548335","8":"2.831158","9":"0.005184537","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## temperature 

```r
nas.1.p2a %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["TEMPERATURE"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["REGION"],"name":[2],"type":["fct"],"align":["left"]},{"label":["emmean"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"29.09659","2":"Core","3":"10.266519","4":"0.4706466","5":"174","6":"9.337607","7":"11.195430","8":"21.81365","9":"1.152724e-51","_rn_":"1"},{"1":"29.09659","2":"Leading","3":"8.175736","4":"0.5182601","5":"174","6":"7.152850","7":"9.198621","8":"15.77535","9":"2.202822e-35","_rn_":"2"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


## Means - f(temperature)

```r
nas.1.p2a %>% update(.~1+ REGION * as.factor(TEMPERATURE) + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID)) %>% 
  emmeans(~REGION*TEMPERATURE, type = "response") %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["REGION"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["emmean"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"Core","2":"27.0","3":"8.319363","4":"0.5374013","5":"174","6":"7.258698","7":"9.380027","8":"15.48073","9":"1.507435e-34","_rn_":"1"},{"1":"Leading","2":"27.0","3":"8.368538","4":"0.6002530","5":"174","6":"7.183824","7":"9.553252","8":"13.94169","9":"3.790369e-30","_rn_":"2"},{"1":"Core","2":"28.5","3":"9.901579","4":"0.5337707","5":"174","6":"8.848080","7":"10.955078","8":"18.55025","9":"4.381380e-43","_rn_":"3"},{"1":"Leading","2":"28.5","3":"8.374571","4":"0.5959145","5":"174","6":"7.198420","7":"9.550722","8":"14.05331","9":"1.810995e-30","_rn_":"4"},{"1":"Core","2":"30.0","3":"10.465447","4":"0.5811305","5":"174","6":"9.318475","7":"11.612420","8":"18.00877","9":"1.311792e-41","_rn_":"5"},{"1":"Leading","2":"30.0","3":"7.830956","4":"0.6367256","5":"174","6":"6.574256","7":"9.087656","8":"12.29879","9":"2.033216e-25","_rn_":"6"},{"1":"Core","2":"31.5","3":"10.171357","4":"0.5949223","5":"174","6":"8.997164","7":"11.345550","8":"17.09695","9":"4.309548e-39","_rn_":"7"},{"1":"Leading","2":"31.5","3":"7.377908","4":"0.6439594","5":"174","6":"6.106931","7":"8.648886","8":"11.45710","9":"5.271775e-23","_rn_":"8"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## Abs. diff - f(temperature)

```r
nas.1.p2a %>% update(.~1+ REGION * as.factor(TEMPERATURE) + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID)) %>% 
  emmeans(~REGION*TEMPERATURE, type = "response") %>% pairs(by ="REGION") %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["REGION"],"name":[2],"type":["fct"],"align":["left"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"TEMPERATURE27 - TEMPERATURE28.5","2":"Core","3":"-1.582216316","4":"0.6413915","5":"174","6":"-3.2460061","7":"0.081573433","8":"-2.466849333","9":"0.06880331","_rn_":"1"},{"1":"TEMPERATURE27 - TEMPERATURE30","2":"Core","3":"-2.146084736","4":"0.7011379","5":"174","6":"-3.9648585","7":"-0.327310983","8":"-3.060859522","9":"0.01351335","_rn_":"2"},{"1":"TEMPERATURE27 - TEMPERATURE31.5","2":"Core","3":"-1.851993894","4":"0.7106950","5":"174","6":"-3.6955590","7":"-0.008428825","8":"-2.605891266","9":"0.04849296","_rn_":"3"},{"1":"TEMPERATURE28.5 - TEMPERATURE30","2":"Core","3":"-0.563868420","4":"0.6790505","5":"174","6":"-2.3253466","7":"1.197609751","8":"-0.830377775","9":"0.83997841","_rn_":"4"},{"1":"TEMPERATURE28.5 - TEMPERATURE31.5","2":"Core","3":"-0.269777578","4":"0.6885327","5":"174","6":"-2.0558529","7":"1.516297771","8":"-0.391815205","9":"0.97952846","_rn_":"5"},{"1":"TEMPERATURE30 - TEMPERATURE31.5","2":"Core","3":"0.294090842","4":"0.6955968","5":"174","6":"-1.5103090","7":"2.098490676","8":"0.422789246","9":"0.97452238","_rn_":"6"},{"1":"TEMPERATURE27 - TEMPERATURE28.5","2":"Leading","3":"-0.006033198","4":"0.6806629","5":"174","6":"-1.7716941","7":"1.759627718","8":"-0.008863709","9":"0.99999975","_rn_":"7"},{"1":"TEMPERATURE27 - TEMPERATURE30","2":"Leading","3":"0.537581507","4":"0.7521353","5":"174","6":"-1.4134809","7":"2.488643910","8":"0.714740466","9":"0.89114214","_rn_":"8"},{"1":"TEMPERATURE27 - TEMPERATURE31.5","2":"Leading","3":"0.990629420","4":"0.7523973","5":"174","6":"-0.9611127","7":"2.942371519","8":"1.316630775","9":"0.55360508","_rn_":"9"},{"1":"TEMPERATURE28.5 - TEMPERATURE30","2":"Leading","3":"0.543614705","4":"0.7490150","5":"174","6":"-1.3993537","7":"2.486583068","8":"0.725772783","9":"0.88668113","_rn_":"10"},{"1":"TEMPERATURE28.5 - TEMPERATURE31.5","2":"Leading","3":"0.996662618","4":"0.7497366","5":"174","6":"-0.9481776","7":"2.941502840","8":"1.329350343","9":"0.54553766","_rn_":"11"},{"1":"TEMPERATURE30 - TEMPERATURE31.5","2":"Leading","3":"0.453047912","4":"0.7406739","5":"174","6":"-1.4682834","7":"2.374379199","8":"0.611669872","9":"0.92826253","_rn_":"12"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
# {-}

# Summary figure 


```
## Warning: Removed 4 rows containing missing values (`geom_point()`).
```

![](aerobic_physiology_absolute_scope_files/figure-html/sum-fig-1.png)<!-- -->

# Conclusion 

* In conclusion while maximum oxygen consumption is **significantly** positively correlated with temperature and fish from low latitudes have **significantly** higher maximum consumption at elevated temperatures compared to fish from high latitudes.
