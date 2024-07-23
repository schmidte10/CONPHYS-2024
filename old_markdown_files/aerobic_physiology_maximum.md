---
title: "Maximum Oxygen Consumption"
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
ggplot(resp4, aes(MASS, MAX_MgO2.hr_RESPR)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_classic()
```

![](aerobic_physiology_maximum_files/figure-html/eda-1-1.png)<!-- -->

## Mass v REST (LATITUDE)

```r
ggplot(resp4, aes(MASS, MAX_MgO2.hr_RESPR, color = REGION)) + 
  geom_point() +
  theme_classic() + 
  geom_smooth(method = "lm")
```

![](aerobic_physiology_maximum_files/figure-html/eda-2-1.png)<!-- -->

## TEMPERTURE v REST (LATITUDE)

```r
ggplot(resp4, aes(TEMPERATURE, MAX_MgO2.hr_RESPR, color = REGION)) + 
  geom_point() +
  theme_classic()
```

![](aerobic_physiology_maximum_files/figure-html/eda-3-1.png)<!-- -->

## {-}

# Fit the model 

The model was fit using the **glm** and later **glmmTMB** package in R. A number of different models were tested to determine which hypothesis and associated variables best predicted resting oxygen consumption. Model fit was examined using AICc, BIC, and r-squared values. Additional model were examined via the validation diagonistics provided by the **performance** and **dHARMA** packages in R. 

The first set of models tested looked at three different hypotheses including 1) that mass has a major impact of resting oxygen consumption of fish (this has been documented in the literature), 2) if variables related to time have an important impact on the resting oxygen consumption of fish. 

## Fixed factors (linear regression models)

### model 1

```r
#--- base model ---#
mmr.1 <- glm(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
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
   <td style="text-align:right;"> -4.5748341 </td>
   <td style="text-align:right;"> 4.4522838 </td>
   <td style="text-align:right;"> -1.027525 </td>
   <td style="text-align:right;"> 0.3056245 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> 17.6468408 </td>
   <td style="text-align:right;"> 6.5852635 </td>
   <td style="text-align:right;"> 2.679747 </td>
   <td style="text-align:right;"> 0.0080880 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TEMPERATURE </td>
   <td style="text-align:right;"> 0.6944978 </td>
   <td style="text-align:right;"> 0.1530872 </td>
   <td style="text-align:right;"> 4.536615 </td>
   <td style="text-align:right;"> 0.0000107 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 298.1596890 </td>
   <td style="text-align:right;"> 24.1436460 </td>
   <td style="text-align:right;"> 12.349406 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:TEMPERATURE </td>
   <td style="text-align:right;"> -0.6609056 </td>
   <td style="text-align:right;"> 0.2261837 </td>
   <td style="text-align:right;"> -2.921986 </td>
   <td style="text-align:right;"> 0.0039482 </td>
  </tr>
</tbody>
</table>

### model 2

```r
#--- experimental rmr equipment hypothesis ---#
mmr.2 <- glm(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MAX_SUMP + MAX_CHAMBER + 
                   MAX_AM_PM, 
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
   <td style="text-align:right;"> -6.2242856 </td>
   <td style="text-align:right;"> 6.4786135 </td>
   <td style="text-align:right;"> -0.9607435 </td>
   <td style="text-align:right;"> 0.3380701 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> 16.0665387 </td>
   <td style="text-align:right;"> 9.1812310 </td>
   <td style="text-align:right;"> 1.7499330 </td>
   <td style="text-align:right;"> 0.0819666 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TEMPERATURE </td>
   <td style="text-align:right;"> 0.7876835 </td>
   <td style="text-align:right;"> 0.2182937 </td>
   <td style="text-align:right;"> 3.6083665 </td>
   <td style="text-align:right;"> 0.0004069 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAX_SUMP2 </td>
   <td style="text-align:right;"> -0.1739027 </td>
   <td style="text-align:right;"> 0.5317501 </td>
   <td style="text-align:right;"> -0.3270383 </td>
   <td style="text-align:right;"> 0.7440484 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAX_CHAMBER2 </td>
   <td style="text-align:right;"> 0.9203132 </td>
   <td style="text-align:right;"> 0.7549758 </td>
   <td style="text-align:right;"> 1.2189970 </td>
   <td style="text-align:right;"> 0.2245644 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAX_CHAMBER3 </td>
   <td style="text-align:right;"> 1.1525268 </td>
   <td style="text-align:right;"> 0.7667736 </td>
   <td style="text-align:right;"> 1.5030861 </td>
   <td style="text-align:right;"> 0.1347055 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAX_CHAMBER4 </td>
   <td style="text-align:right;"> 0.5061699 </td>
   <td style="text-align:right;"> 0.8313336 </td>
   <td style="text-align:right;"> 0.6088649 </td>
   <td style="text-align:right;"> 0.5434413 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAX_AM_PMPM </td>
   <td style="text-align:right;"> -0.4757712 </td>
   <td style="text-align:right;"> 0.5382845 </td>
   <td style="text-align:right;"> -0.8838656 </td>
   <td style="text-align:right;"> 0.3780394 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:TEMPERATURE </td>
   <td style="text-align:right;"> -0.7177722 </td>
   <td style="text-align:right;"> 0.3148211 </td>
   <td style="text-align:right;"> -2.2799368 </td>
   <td style="text-align:right;"> 0.0238764 </td>
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
   <td style="text-align:left;"> mmr.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 828.4562 </td>
   <td style="text-align:right;"> 846.9821 </td>
   <td style="text-align:right;"> 0.6622066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmr.2 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 945.6551 </td>
   <td style="text-align:right;"> 976.0266 </td>
   <td style="text-align:right;"> 0.3732903 </td>
  </tr>
</tbody>
</table>

The model that contains **MASS_CENTERED** seems to do better than the model that incorporates variables that are associated with the time that experiments are performed.This is demonstrated by the lower AIC and BIC scores, as well as higher r-squared value.  

## Polynomials 

### polynomial models 

Note that the linear model has already been created via model _rmr.3_ in the previous section.


```r
mmr.1.p2 <- glm(MAX_MgO2.hr_RESPR ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4)

mmr.1.p3 <- glm(MAX_MgO2.hr_RESPR ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED, 
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
   <td style="text-align:left;"> mmr.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 828.4562 </td>
   <td style="text-align:right;"> 846.9821 </td>
   <td style="text-align:right;"> 0.6673593 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmr.1.p2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 832.3856 </td>
   <td style="text-align:right;"> 856.8872 </td>
   <td style="text-align:right;"> 0.6681820 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmr.1.p3 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 836.8419 </td>
   <td style="text-align:right;"> 867.2134 </td>
   <td style="text-align:right;"> 0.6682098 </td>
  </tr>
</tbody>
</table>

From our model comparison we can see the there is no additional benefit to the model by including temperature as a 2^nd^ or 3^rd^ order polynomial. However, the linear and quadratic model both perform well. 

## Random factors 

Fish were repeatedly sampled over four different temperatures, therefore repeated sampling needs to be accounted for. To do this random factors will be included within the model. There are a number of options that can be used for random factors including 1) accounting for repeated sampling of individuals, 2) accounting for repeated sampling of individuals nested within population, 3) account for repeated sampling of individuals and populations without nesting. All three models will be run a compaired. 

### random factor models


```r
mmr.1a <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE) 

mmr.1b <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|POPULATION/FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

mmr.1c <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID) + (REGION|POPULATION), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE) # Convergence problem
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
   <td style="text-align:left;"> mmr.1a </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 809.5847 </td>
   <td style="text-align:right;"> 831.1115 </td>
   <td style="text-align:right;"> 0.6693399 </td>
   <td style="text-align:right;"> 0.7836732 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmr.1b </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 811.5475 </td>
   <td style="text-align:right;"> 836.0491 </td>
   <td style="text-align:right;"> 0.6694043 </td>
   <td style="text-align:right;"> 0.7848247 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmr.1c </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.6697666 </td>
   <td style="text-align:right;"> 0.7848011 </td>
  </tr>
</tbody>
</table>

Model _rmr.3a_ appears to be the best model, however, there seems to be little difference in how the models change depending on how the random factors are arranged.

# Model validation {.tabset .tabset-faded}

## performance {.tabset .tabset-faded}

### rmr.3a (linear)
![](aerobic_physiology_maximum_files/figure-html/model-valid-1-1.png)<!-- -->

The _rmr.3a_ model performs well, however, in the model validation performed by the **performance** model it looks like there are two variables that are highly correlated. If we expand the figure we can see that the highly correlated variables are REGION and REGION:TEMPERATURE. Perhaps this is unsurprising  but lets see what happens when we run the quadratic (2^nd^ polynomial) model to see if this helps deal with the high correlation between these two variables, as it performed very similarly to _rmr.3a_, and even had a higher r2 value. 

### rmr.3.p2a (quadratic)

First we need to update the model by adding in the missing random factor

```r
mmr.1.p2a <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED + (1|FISH_ID), 
                 family=gaussian(),
                 data = resp4, 
                 REML = TRUE) 
```

![](aerobic_physiology_maximum_files/figure-html/model-valid-1.2-1.png)<!-- -->

## DHARMa residuals {.tabset .tabset-faded}

### rmr.3a (linear)

```r
mmr.1a %>% simulateResiduals(plot=TRUE)
```

![](aerobic_physiology_maximum_files/figure-html/model-valid-2-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.212 0.328 0.456 0.744 0.812 0.852 1 0.7 0.372 0.684 0.692 0.032 0.216 0.204 0.612 0.084 0.452 0.988 0.936 0.848 ...
```

```r
mmr.1a %>% DHARMa::testResiduals(plot=TRUE)
```

![](aerobic_physiology_maximum_files/figure-html/model-valid-2-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.0715, p-value = 0.3293
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.94899, p-value = 0.68
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 3, observations = 176, p-value = 0.1665
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.003529072 0.049003686
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01704545
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.0715, p-value = 0.3293
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.94899, p-value = 0.68
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 3, observations = 176, p-value = 0.1665
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.003529072 0.049003686
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01704545
```

### rmr.3.p2 (quadratic)

First we need to update the model by adding in the missing random factor

```r
mmr.1.p2a <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED + (1|FISH_ID), 
                 family=gaussian(),
                 data = resp4, 
                 REML = TRUE) 
```


```r
mmr.1.p2a %>% simulateResiduals(plot=TRUE) 
```

![](aerobic_physiology_maximum_files/figure-html/model-valid-2.2-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.228 0.308 0.492 0.728 0.788 0.876 1 0.728 0.356 0.664 0.724 0.04 0.188 0.192 0.644 0.092 0.42 0.98 0.936 0.836 ...
```

```r
mmr.1.p2a %>% DHARMa::testResiduals(plot=TRUE)
```

![](aerobic_physiology_maximum_files/figure-html/model-valid-2.2-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.076136, p-value = 0.2594
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.94042, p-value = 0.608
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 3, observations = 176, p-value = 0.1665
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.003529072 0.049003686
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01704545
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.076136, p-value = 0.2594
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.94042, p-value = 0.608
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 3, observations = 176, p-value = 0.1665
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.003529072 0.049003686
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01704545
```

## {-}

# {-}

It looks like the model that treats temperature as a second order polynomial does a better job at avoiding high levels of collinearity within the model. The quadratic model will be used moving forward because it: 

* The **quadratic model** performs just as well as the linear model based on the model validation scores (i.e., AIC, BIC, and r2) 
* The **quadratic model** does a **better** job at dealing with collinearity that appeared in the model 

# Partial plots {.tabset .tabset-faded}

## ggemmeans 

![](aerobic_physiology_maximum_files/figure-html/partial-plots-1-1.png)<!-- -->

## plot_model 

![](aerobic_physiology_maximum_files/figure-html/partial-plots-2-1.png)<!-- -->

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
   <td style="text-align:right;"> 15.665104 </td>
   <td style="text-align:right;"> 0.3860064 </td>
   <td style="text-align:right;"> 40.5824948 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -1.540213 </td>
   <td style="text-align:right;"> 0.6378566 </td>
   <td style="text-align:right;"> -2.4146702 </td>
   <td style="text-align:right;"> 0.0157495 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> 14.695118 </td>
   <td style="text-align:right;"> 2.8646793 </td>
   <td style="text-align:right;"> 5.1297602 </td>
   <td style="text-align:right;"> 0.0000003 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -2.507947 </td>
   <td style="text-align:right;"> 2.8404686 </td>
   <td style="text-align:right;"> -0.8829342 </td>
   <td style="text-align:right;"> 0.3772718 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 313.864413 </td>
   <td style="text-align:right;"> 34.0009493 </td>
   <td style="text-align:right;"> 9.2310485 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> -13.979779 </td>
   <td style="text-align:right;"> 4.2142722 </td>
   <td style="text-align:right;"> -3.3172464 </td>
   <td style="text-align:right;"> 0.0009091 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> 1.659481 </td>
   <td style="text-align:right;"> 4.1665787 </td>
   <td style="text-align:right;"> 0.3982838 </td>
   <td style="text-align:right;"> 0.6904210 </td>
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
   <td style="text-align:right;"> 5.244328 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0220184 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2) </td>
   <td style="text-align:right;"> 16.284260 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0002910 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 85.212257 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGION:poly(TEMPERATURE, 2) </td>
   <td style="text-align:right;"> 11.182851 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0037297 </td>
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
   <td style="text-align:right;"> 14.908545 </td>
   <td style="text-align:right;"> 16.4216626 </td>
   <td style="text-align:right;"> 15.665104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -2.790389 </td>
   <td style="text-align:right;"> -0.2900373 </td>
   <td style="text-align:right;"> -1.540213 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> 9.080450 </td>
   <td style="text-align:right;"> 20.3097864 </td>
   <td style="text-align:right;"> 14.695118 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -8.075163 </td>
   <td style="text-align:right;"> 3.0592692 </td>
   <td style="text-align:right;"> -2.507947 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS_CENTERED </td>
   <td style="text-align:right;"> 247.223777 </td>
   <td style="text-align:right;"> 380.5050495 </td>
   <td style="text-align:right;"> 313.864413 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> -22.239601 </td>
   <td style="text-align:right;"> -5.7199575 </td>
   <td style="text-align:right;"> -13.979779 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -6.506863 </td>
   <td style="text-align:right;"> 9.8258252 </td>
   <td style="text-align:right;"> 1.659481 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Std.Dev.(Intercept)|FISH_ID </td>
   <td style="text-align:right;"> 1.054358 </td>
   <td style="text-align:right;"> 2.1168357 </td>
   <td style="text-align:right;"> 1.493956 </td>
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
   <td style="text-align:right;"> 0.7829506 </td>
   <td style="text-align:right;"> 0.6683685 </td>
   <td style="text-align:left;"> FALSE </td>
  </tr>
</tbody>
</table>

# {-} 

# Pairwise comparisons {.tabset .tabset-faded} 

## emtrends [latitudes]



```r
mmr.1.p2a %>% emtrends(var = "TEMPERATURE", type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"(Core MASS_CENTERED-0.000127610645933016) - (Leading MASS_CENTERED-0.000127610645933016)","2":"29.09659","3":"0.6432331","4":"0.1924589","5":"174","6":"0.2633786","7":"1.023088","8":"3.342184","9":"0.001017531","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
SCROLL TO THE RIGHT -->

The numbers in the left most column in the table just mention that the slopes are assuming mean **MASS_CENTERED** values when looking at differences between latitudinal slopes.

## emmeans [latitudes]

```r
mmr.1.p2a %>% emmeans(pairwise ~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"Core - Leading","2":"29.09659","3":"1.695557","4":"0.7435233","5":"174","6":"0.2280712","7":"3.163042","8":"2.280435","9":"0.02379491","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## temperature 

```r
mmr.1.p2a %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["TEMPERATURE"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["REGION"],"name":[2],"type":["fct"],"align":["left"]},{"label":["emmean"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"29.09659","2":"Core","3":"15.85982","4":"0.4704955","5":"174","6":"14.93121","7":"16.78843","8":"33.70876","9":"3.378353e-78","_rn_":"1"},{"1":"29.09659","2":"Leading","3":"14.16426","4":"0.5214344","5":"174","6":"13.13511","7":"15.19341","8":"27.16403","9":"1.735592e-64","_rn_":"2"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


## Means - f(temperature)

```r
mmr.1.p2a %>% update(.~1+ REGION * as.factor(TEMPERATURE) + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID)) %>% 
  emmeans(~REGION*TEMPERATURE, type = "response") %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["REGION"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["emmean"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"Core","2":"27.0","3":"13.90434","4":"0.5275303","5":"174","6":"12.86315","7":"14.94552","8":"26.35742","9":"1.187092e-62","_rn_":"1"},{"1":"Leading","2":"27.0","3":"13.83242","4":"0.5923828","5":"174","6":"12.66324","7":"15.00160","8":"23.35048","9":"1.659212e-55","_rn_":"2"},{"1":"Core","2":"28.5","3":"15.44277","4":"0.5243370","5":"174","6":"14.40789","7":"16.47765","8":"29.45200","9":"1.638937e-69","_rn_":"3"},{"1":"Leading","2":"28.5","3":"14.06788","4":"0.5882491","5":"174","6":"12.90686","7":"15.22891","8":"23.91484","9":"6.943035e-57","_rn_":"4"},{"1":"Core","2":"30.0","3":"16.46500","4":"0.5667728","5":"174","6":"15.34637","7":"17.58364","8":"29.05045","9":"1.195713e-68","_rn_":"5"},{"1":"Leading","2":"30.0","3":"14.25629","4":"0.6251032","5":"174","6":"13.02253","7":"15.49005","8":"22.80629","9":"3.676862e-54","_rn_":"6"},{"1":"Core","2":"31.5","3":"17.12324","4":"0.5798435","5":"174","6":"15.97880","7":"18.26767","8":"29.53079","9":"1.112103e-69","_rn_":"7"},{"1":"Leading","2":"31.5","3":"14.20198","4":"0.6317663","5":"174","6":"12.95507","7":"15.44889","8":"22.47980","9":"2.402005e-53","_rn_":"8"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## Abs. diff - f(temperature)

```r
mmr.1.p2a %>% update(.~1+ REGION * as.factor(TEMPERATURE) + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID)) %>% 
  emmeans(~REGION*TEMPERATURE, type = "response") %>% pairs(by ="REGION") %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["REGION"],"name":[2],"type":["fct"],"align":["left"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"TEMPERATURE27 - TEMPERATURE28.5","2":"Core","3":"-1.53843734","4":"0.5938799","5":"174","6":"-3.078980","7":"0.00210572","8":"-2.59048559","9":"5.045819e-02","_rn_":"1"},{"1":"TEMPERATURE27 - TEMPERATURE30","2":"Core","3":"-2.56066721","4":"0.6492325","5":"174","6":"-4.244797","7":"-0.87653775","8":"-3.94414504","9":"6.654002e-04","_rn_":"2"},{"1":"TEMPERATURE27 - TEMPERATURE31.5","2":"Core","3":"-3.21889902","4":"0.6599391","5":"174","6":"-4.930802","7":"-1.50699640","8":"-4.87756996","9":"1.427813e-05","_rn_":"3"},{"1":"TEMPERATURE28.5 - TEMPERATURE30","2":"Core","3":"-1.02222987","4":"0.6287983","5":"174","6":"-2.653352","7":"0.60889261","8":"-1.62568801","9":"3.669065e-01","_rn_":"4"},{"1":"TEMPERATURE28.5 - TEMPERATURE31.5","2":"Core","3":"-1.68046168","4":"0.6392790","5":"174","6":"-3.338772","7":"-0.02215185","8":"-2.62868259","9":"4.570512e-02","_rn_":"5"},{"1":"TEMPERATURE30 - TEMPERATURE31.5","2":"Core","3":"-0.65823181","4":"0.6452530","5":"174","6":"-2.332038","7":"1.01557474","8":"-1.02011425","9":"7.378999e-01","_rn_":"6"},{"1":"TEMPERATURE27 - TEMPERATURE28.5","2":"Leading","3":"-0.23546008","4":"0.6295418","5":"174","6":"-1.868511","7":"1.39759118","8":"-0.37401815","9":"9.821016e-01","_rn_":"7"},{"1":"TEMPERATURE27 - TEMPERATURE30","2":"Leading","3":"-0.42386353","4":"0.6960953","5":"174","6":"-2.229557","7":"1.38182947","8":"-0.60891595","9":"9.291393e-01","_rn_":"8"},{"1":"TEMPERATURE27 - TEMPERATURE31.5","2":"Leading","3":"-0.36955401","4":"0.6966401","5":"174","6":"-2.176660","7":"1.43755227","8":"-0.53048052","9":"9.515685e-01","_rn_":"9"},{"1":"TEMPERATURE28.5 - TEMPERATURE30","2":"Leading","3":"-0.18840346","4":"0.6933597","5":"174","6":"-1.987000","7":"1.61019341","8":"-0.27172541","9":"9.929633e-01","_rn_":"10"},{"1":"TEMPERATURE28.5 - TEMPERATURE31.5","2":"Leading","3":"-0.13409394","4":"0.6943276","5":"174","6":"-1.935202","7":"1.66701369","8":"-0.19312775","9":"9.974383e-01","_rn_":"11"},{"1":"TEMPERATURE30 - TEMPERATURE31.5","2":"Leading","3":"0.05430952","4":"0.6839527","5":"174","6":"-1.719885","7":"1.82850415","8":"0.07940538","9":"9.998198e-01","_rn_":"12"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
# {-}

# Summary figure 


```
## Warning: Removed 1 rows containing missing values (`geom_point()`).
```

![](aerobic_physiology_maximum_files/figure-html/sum-fig-1.png)<!-- -->

# Conclusion 

* In conclusion while maximum oxygen consumption is **significantly** positively correlated with temperature and fish from low latitudes have **significantly** higher maximum consumption at elevated temperatures compared to fish from high latitudes.
