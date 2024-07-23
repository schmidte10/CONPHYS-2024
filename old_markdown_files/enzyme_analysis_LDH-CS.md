---
title: "Citrate Synthase (CS)"
author: "Elliott Schmidt"
date: "07 December, 2023"
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

For initial details on the experiment performed please read the **ReadMe** file. In brief, _Acanthochromis polyacanthus_ from two different regions on the Great Barrier Reef (GBR) were tested for metabolic performance at four different temperatures, 27$^\circ$C, 28.5$^\circ$C, 30$^\circ$C, and 31.5$^\circ$C. Fish used in this study were collected from two different regions, low- (i.e. Cairns) and high-latitude (i.e., Mackay), within each region fish were collected from a total of three different populations. After metabolic performance was tested blood and tissue samples were collected. White muscle tissue samples were used to look at the relationship between activity and temperature in two different enzymes, Lactate Dehydrogenase (LDH; anaerobic) and Citrate Synthase (CS: aerobic). Enzyme activity was measured over four different temperatures including 20$^\circ$C, 30$^\circ$C, 40$^\circ$C, and 50$^\circ$C. Enzyme activity was measured using a spectophotometer and wavelength absoprtion levels were recorded using the software program LabX. 

This analysis is looking at the LDH:CS ratio. Methods used to determine LDH and CS activity level have been previously outlined in separate markdown sheets for each enzyme, and therefore will be repeated here. Within this worksheet LDH and CS data frames will be imported from saved files in a directory.

# Load packages 

Lets start by loading the packages that are needed 

```r
library(tidyverse) # data manipulation
library(janitor) # data manipulation
library(plyr) # data manipulation
library(dplyr) # data manipulation
library(lubridate) # data manipulation - specifically time data
library(ggplot2) # plotting figures
library(glmmTMB) # running models
library(performance) # model validation
library(chron) # data manipulation - specifically time data
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
library(car) # used for Anova function
```

## Read in the data

Before beginning always make sure that you are working in the correct directory 




```r
knitr::opts_knit$set(root.dir=working.dir)
```

Now we can import that data. Two different data frames are being imported. The first has all the enzyme wave length absorption data for each sample and the tissue.mass data file contained information pertaining to the tissue samples that was used for each sample. Later on these two data frames will be merged. 

## Load data 


# Data manipulation 

Both LDH and CS dataframes have been previously cleaned and manipulated, therefore, the only remaining step is to join the data frames together and then make another column that has the LDH:CS ratio


```r
#--- data preparation/manipulation ---# 
ldh.cs.data <- ldh.data %>% 
  inner_join(select(cs.data, c("UNIQUE_SAMPLE_ID","CS_ACTIVITY")), by = "UNIQUE_SAMPLE_ID") %>% 
  mutate(LCr = LDH_ACTIVITY/CS_ACTIVITY)
```

# Exploratory data analysis {.tabset}

## LDH v TEMPERATURE [LATITUDE]

```r
ggplot(ldh.cs.data, aes(x =as.numeric(temperature), y= LCr, color = REGION)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)
```

![](enzyme_analysis_LDH-CS_files/figure-html/eda-1-1.png)<!-- -->

## LDH V TEMPERATURE [DENSITY]

```r
ggplot(ldh.cs.data, aes(x = LCr)) + 
  geom_density(alpha =0.5, position = "identity") 
```

![](enzyme_analysis_LDH-CS_files/figure-html/eda-2-1.png)<!-- -->

## LDH v TISSUE MASS (LATITUDE)

```r
ggplot(ldh.cs.data, aes(x =TISSUE_MASS_CENTERED, y= LCr, color = REGION)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)
```

![](enzyme_analysis_LDH-CS_files/figure-html/eda-3-1.png)<!-- -->


## {-}

# Fit the model 

The model was fit using the **glm** and later **glmmTMB** package in R. A number of different models were tested to determine which hypothesis and associated variables best predicted resting oxygen consumption. Model fit was examined using AICc, BIC, and r-squared values. Additional model were examined via the validation diagonistics provided by the **performance** and **dHARMA** packages in R. 

## Fixed factors (linear regression models)

### model 1

```r
#--- base model ---#
ldh.cs.model.1 <- glm(LCr~ 1 + REGION*temperature + TISSUE_MASS_CENTERED, 
                       family=gaussian(), 
                       data = ldh.cs.data)  
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
   <td style="text-align:right;"> 15.5825978 </td>
   <td style="text-align:right;"> 5.0803058 </td>
   <td style="text-align:right;"> 3.0672558 </td>
   <td style="text-align:right;"> 0.0026407 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -3.9697588 </td>
   <td style="text-align:right;"> 7.4743176 </td>
   <td style="text-align:right;"> -0.5311199 </td>
   <td style="text-align:right;"> 0.5962633 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> temperature </td>
   <td style="text-align:right;"> 0.4263426 </td>
   <td style="text-align:right;"> 0.1403591 </td>
   <td style="text-align:right;"> 3.0375120 </td>
   <td style="text-align:right;"> 0.0028959 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TISSUE_MASS_CENTERED </td>
   <td style="text-align:right;"> -0.5452691 </td>
   <td style="text-align:right;"> 0.2308206 </td>
   <td style="text-align:right;"> -2.3623065 </td>
   <td style="text-align:right;"> 0.0196809 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:temperature </td>
   <td style="text-align:right;"> -0.0073350 </td>
   <td style="text-align:right;"> 0.2047640 </td>
   <td style="text-align:right;"> -0.0358219 </td>
   <td style="text-align:right;"> 0.9714806 </td>
  </tr>
</tbody>
</table>

### model 2

```r
ldh.cs.model.2 <- glm(LCr ~ 1 + REGION*temperature, 
                       family=gaussian(), 
                       data = ldh.cs.data) 
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
   <td style="text-align:right;"> 14.9398547 </td>
   <td style="text-align:right;"> 5.1629861 </td>
   <td style="text-align:right;"> 2.8936461 </td>
   <td style="text-align:right;"> 0.0044775 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -2.9607306 </td>
   <td style="text-align:right;"> 7.5944452 </td>
   <td style="text-align:right;"> -0.3898548 </td>
   <td style="text-align:right;"> 0.6972919 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> temperature </td>
   <td style="text-align:right;"> 0.4346427 </td>
   <td style="text-align:right;"> 0.1428037 </td>
   <td style="text-align:right;"> 3.0436375 </td>
   <td style="text-align:right;"> 0.0028374 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:temperature </td>
   <td style="text-align:right;"> -0.0117352 </td>
   <td style="text-align:right;"> 0.2083869 </td>
   <td style="text-align:right;"> -0.0563143 </td>
   <td style="text-align:right;"> 0.9551792 </td>
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
   <td style="text-align:left;"> ldh.cs.model.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1058.448 </td>
   <td style="text-align:right;"> 1075.073 </td>
   <td style="text-align:right;"> 0.1609504 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ldh.cs.model.2 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1061.928 </td>
   <td style="text-align:right;"> 1075.866 </td>
   <td style="text-align:right;"> 0.1259275 </td>
  </tr>
</tbody>
</table>

The model that contains **TISSUE_MASS_CENTERED** seems to do better than the model that leaves TISSUE_MASS_CENTERED out. Therefore we will move ahead with the model that contains **TISSUE_MASS_CENTERED** as a co-variate.  

## Polynomials 

### polynomial models 

Note that the linear model has already been created via model _ldh.cs.model.1_ in the previous section.


```r
#--- second order polynomial ---# 
ldh.cs.model.1.p2 <- glm(LCr ~ 1 + REGION*poly(temperature, 2) + TISSUE_MASS_CENTERED, 
                      family=gaussian(), 
                      data = ldh.cs.data)  

#--- third order polynomial ---# 
ldh.cs.model.1.p3 <- glm(LCr ~ 1 + REGION*poly(temperature, 3) + TISSUE_MASS_CENTERED, 
                      family=gaussian(), 
                      data = ldh.cs.data)  
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
   <td style="text-align:left;"> ldh.cs.model.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1058.448 </td>
   <td style="text-align:right;"> 1075.073 </td>
   <td style="text-align:right;"> 0.1651823 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ldh.cs.model.1.p2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 1062.071 </td>
   <td style="text-align:right;"> 1083.963 </td>
   <td style="text-align:right;"> 0.1707019 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ldh.cs.model.1.p3 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 1064.193 </td>
   <td style="text-align:right;"> 1091.203 </td>
   <td style="text-align:right;"> 0.1864147 </td>
  </tr>
</tbody>
</table>

From our model comparison we can see that the model that runs temperature as a linear model performs the best. Therefore, moving forward we will use the linear model. 

## Random factors 

Fish were repeatedly sampled over four different temperatures, therefore repeated sampling needs to be accounted for. To do this random factors will be included within the model. There are a number of options that can be used for random factors including 1) accounting for repeated sampling of individuals, 2) accounting for repeated sampling of individuals nested within population, 3) account for repeated sampling of individuals and populations without nesting. All three models will be run a compaired. 

### random factor models


```r
ldh.cs.model.1a <- glmmTMB(LCr ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|fish_id), 
                       family=gaussian(), 
                       data = ldh.cs.data, 
                       REML = TRUE) 

ldh.cs.model.1b <- glmmTMB(LCr ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|POPULATION/fish_id), 
                       family=gaussian(), 
                       data = ldh.cs.data,
                       REML = TRUE) 

ldh.cs.model.1c <- glmmTMB(LCr ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|fish_id) + (1 + REGION|POPULATION), 
                       family=gaussian(), 
                       data = ldh.cs.data,
                       REML = TRUE) # convergnece problem
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
   <td style="text-align:left;"> ldh.cs.model.1a </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 994.5888 </td>
   <td style="text-align:right;"> 1013.865 </td>
   <td style="text-align:right;"> 0.1483105 </td>
   <td style="text-align:right;"> 0.1483105 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ldh.cs.model.1b </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 996.8563 </td>
   <td style="text-align:right;"> 1018.748 </td>
   <td style="text-align:right;"> 0.1483098 </td>
   <td style="text-align:right;"> 0.1483098 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ldh.cs.model.1c </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 1001.5037 </td>
   <td style="text-align:right;"> 1028.514 </td>
   <td style="text-align:right;"> 0.1483106 </td>
   <td style="text-align:right;"> 0.1483106 </td>
  </tr>
</tbody>
</table>

Model _ldh.cs.model.1a_ appears to be the best model, however, there seems to be little difference in how the models change depending on how the random factors are arranged.

# Model validation {.tabset .tabset-faded}

## performance {.tabset .tabset-faded}

### rmr.3a (linear)
![](enzyme_analysis_LDH-CS_files/figure-html/model-valid-1-1.png)<!-- -->

## DHARMa residuals {.tabset .tabset-faded}

### nas.1a (linear)

```r
ldh.cs.model.1a %>% simulateResiduals(plot=TRUE)
```

![](enzyme_analysis_LDH-CS_files/figure-html/model-valid-2-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.072 0.052 0.076 0.152 0.944 0.512 0.9 0.864 0.708 0.692 0.436 0.66 0.764 0.928 0.984 0.272 0.136 0.388 1 0.868 ...
```

```r
ldh.cs.model.1a %>% DHARMa::testResiduals(plot=TRUE)
```

![](enzyme_analysis_LDH-CS_files/figure-html/model-valid-2-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.098, p-value = 0.1584
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.91404, p-value = 0.624
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 132, p-value = 0.2834
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001840215 0.053659931
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01515152
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.098, p-value = 0.1584
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.91404, p-value = 0.624
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 132, p-value = 0.2834
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001840215 0.053659931
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01515152
```

## {-}

# {-}

The _ldh.cs.model.1a_ model looks good, and there seem to be no major violations of assumptions. 

# Partial plots {.tabset .tabset-faded}

## ggemmeans 

![](enzyme_analysis_LDH-CS_files/figure-html/partial-plots-1-1.png)<!-- -->

## plot_model 

![](enzyme_analysis_LDH-CS_files/figure-html/partial-plots-2-1.png)<!-- -->

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
   <td style="text-align:right;"> 14.6974700 </td>
   <td style="text-align:right;"> 3.9311000 </td>
   <td style="text-align:right;"> 3.7387677 </td>
   <td style="text-align:right;"> 0.0001849 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -3.3896394 </td>
   <td style="text-align:right;"> 5.8100316 </td>
   <td style="text-align:right;"> -0.5834115 </td>
   <td style="text-align:right;"> 0.5596163 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> temperature </td>
   <td style="text-align:right;"> 0.4426061 </td>
   <td style="text-align:right;"> 0.0836021 </td>
   <td style="text-align:right;"> 5.2941995 </td>
   <td style="text-align:right;"> 0.0000001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TISSUE_MASS_CENTERED </td>
   <td style="text-align:right;"> -0.4573831 </td>
   <td style="text-align:right;"> 0.3990843 </td>
   <td style="text-align:right;"> -1.1460814 </td>
   <td style="text-align:right;"> 0.2517615 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:temperature </td>
   <td style="text-align:right;"> -0.0132587 </td>
   <td style="text-align:right;"> 0.1215687 </td>
   <td style="text-align:right;"> -0.1090631 </td>
   <td style="text-align:right;"> 0.9131524 </td>
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
   <td style="text-align:right;"> 0.9282535 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.3353172 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> temperature </td>
   <td style="text-align:right;"> 51.6720557 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TISSUE_MASS_CENTERED </td>
   <td style="text-align:right;"> 1.3135026 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.2517615 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGION:temperature </td>
   <td style="text-align:right;"> 0.0118948 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.9131524 </td>
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
   <td style="text-align:right;"> 6.9926555 </td>
   <td style="text-align:right;"> 22.4022845 </td>
   <td style="text-align:right;"> 14.6974700 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -14.7770922 </td>
   <td style="text-align:right;"> 7.9978133 </td>
   <td style="text-align:right;"> -3.3896394 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> temperature </td>
   <td style="text-align:right;"> 0.2787490 </td>
   <td style="text-align:right;"> 0.6064632 </td>
   <td style="text-align:right;"> 0.4426061 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TISSUE_MASS_CENTERED </td>
   <td style="text-align:right;"> -1.2395740 </td>
   <td style="text-align:right;"> 0.3248078 </td>
   <td style="text-align:right;"> -0.4573831 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:temperature </td>
   <td style="text-align:right;"> -0.2515289 </td>
   <td style="text-align:right;"> 0.2250116 </td>
   <td style="text-align:right;"> -0.0132587 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Std.Dev.(Intercept)|fish_id </td>
   <td style="text-align:right;"> 8.3076257 </td>
   <td style="text-align:right;"> 14.5083432 </td>
   <td style="text-align:right;"> 10.9786103 </td>
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
   <td style="text-align:right;"> 0.7239123 </td>
   <td style="text-align:right;"> 0.1483105 </td>
   <td style="text-align:left;"> FALSE </td>
  </tr>
</tbody>
</table>

# {-} 

# Pairwise comparisons {.tabset .tabset-faded} 

## emtrends [latitudes]



```r
ldh.cs.model.1a  %>% emtrends(var = "temperature", type = "response") %>% pairs(by = "temperature") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["temperature"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"(Core TISSUE_MASS_CENTERED-0.0859139162323838) - (Leading TISSUE_MASS_CENTERED-0.0859139162323838)","2":"34.69697","3":"0.01325866","4":"0.1215687","5":"130","6":"-0.2272504","7":"0.2537678","8":"0.1090631","9":"0.9133206","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
SCROLL TO THE RIGHT -->

The numbers in the left most column in the table just mention that the slopes are assuming mean **TISSUE_MASS_CENTERED** values when looking at differences between latitudinal slopes.

## emmeans [latitudes]

```r
ldh.cs.model.1a  %>% emmeans(pairwise ~ temperature*REGION, type = "response") %>% pairs(by = "temperature") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["temperature"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"Core - Leading","2":"34.69697","3":"3.849675","4":"3.995668","5":"130","6":"-4.055276","7":"11.75463","8":"0.9634622","9":"0.3371044","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## TEMPERATURE 

```r
ldh.cs.model.1a  %>% emmeans(~ temperature*REGION, type = "response")  %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["temperature"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["REGION"],"name":[2],"type":["fct"],"align":["left"]},{"label":["emmean"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"34.69697","2":"Core","3":"30.09386","4":"2.688083","5":"130","6":"24.77581","7":"35.41191","8":"11.195285","9":"8.699907e-21","_rn_":"1"},{"1":"34.69697","2":"Leading","3":"26.24418","4":"2.931895","5":"130","6":"20.44378","7":"32.04459","8":"8.951268","9":"3.115589e-15","_rn_":"2"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


## Means - f(TEMPERATURE)

```r
ldh.cs.model.1a  %>% update(.~1+ REGION * as.factor(temperature) + TISSUE_MASS_CENTERED + (1|fish_id)) %>% 
  emmeans(~REGION*temperature, type = "response") %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["REGION"],"name":[1],"type":["fct"],"align":["left"]},{"label":["temperature"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["emmean"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"Core","2":"20","3":"25.99800","4":"3.088150","5":"130","6":"19.88847","7":"32.10754","8":"8.418634","9":"6.037693e-14","_rn_":"1"},{"1":"Leading","2":"20","3":"21.05130","4":"3.380730","5":"130","6":"14.36293","7":"27.73967","8":"6.226850","9":"6.089996e-09","_rn_":"2"},{"1":"Core","2":"30","3":"24.43566","4":"3.051838","5":"130","6":"18.39796","7":"30.47335","8":"8.006865","9":"5.774001e-13","_rn_":"3"},{"1":"Leading","2":"30","3":"22.21330","4":"3.337702","5":"130","6":"15.61005","7":"28.81654","8":"6.655267","9":"7.171644e-10","_rn_":"4"},{"1":"Core","2":"40","3":"32.91752","4":"3.125331","5":"130","6":"26.73442","7":"39.10061","8":"10.532489","9":"3.906587e-19","_rn_":"5"},{"1":"Leading","2":"40","3":"29.39876","4":"3.337702","5":"130","6":"22.79552","7":"36.00200","8":"8.808086","9":"6.941184e-15","_rn_":"6"},{"1":"Core","2":"50","3":"37.99695","4":"3.165301","5":"130","6":"31.73478","7":"44.25912","8":"12.004214","9":"8.374284e-23","_rn_":"7"},{"1":"Leading","2":"50","3":"32.88505","4":"3.382059","5":"130","6":"26.19405","7":"39.57605","8":"9.723380","9":"3.974312e-17","_rn_":"8"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## Abs. diff - f(TEMPERATURE)

```r
ldh.cs.model.1a  %>% update(.~1+ REGION * as.factor(temperature) + TISSUE_MASS_CENTERED + (1|fish_id)) %>% 
  emmeans(~REGION*temperature, type = "response") %>% pairs(by ="REGION") %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["REGION"],"name":[2],"type":["fct"],"align":["left"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"temperature20 - temperature30","2":"Core","3":"1.562348","4":"2.455105","5":"130","6":"-4.827244","7":"7.9519400","8":"0.6363669","9":"9.200786e-01","_rn_":"1"},{"1":"temperature20 - temperature40","2":"Core","3":"-6.919511","4":"2.546509","5":"130","6":"-13.546987","7":"-0.2920357","8":"-2.7172543","9":"3.709612e-02","_rn_":"2"},{"1":"temperature20 - temperature50","2":"Core","3":"-11.998947","4":"2.599422","5":"130","6":"-18.764134","7":"-5.2337597","8":"-4.6160055","9":"5.440284e-05","_rn_":"3"},{"1":"temperature30 - temperature40","2":"Core","3":"-8.481859","4":"2.502843","5":"130","6":"-14.995691","7":"-1.9680271","8":"-3.3888903","9":"5.085203e-03","_rn_":"4"},{"1":"temperature30 - temperature50","2":"Core","3":"-13.561295","4":"2.554217","5":"130","6":"-20.208833","7":"-6.9137558","8":"-5.3093734","9":"2.744597e-06","_rn_":"5"},{"1":"temperature40 - temperature50","2":"Core","3":"-5.079436","4":"2.650034","5":"130","6":"-11.976345","7":"1.8174735","8":"-1.9167433","9":"2.261013e-01","_rn_":"6"},{"1":"temperature20 - temperature30","2":"Leading","3":"-1.161996","4":"2.684328","5":"130","6":"-8.148156","7":"5.8241641","8":"-0.4328815","9":"9.727137e-01","_rn_":"7"},{"1":"temperature20 - temperature40","2":"Leading","3":"-8.347462","4":"2.684328","5":"130","6":"-15.333622","7":"-1.3613020","8":"-3.1097031","9":"1.218348e-02","_rn_":"8"},{"1":"temperature20 - temperature50","2":"Leading","3":"-11.833749","4":"2.741209","5":"130","6":"-18.967948","7":"-4.6995505","8":"-4.3169812","9":"1.807901e-04","_rn_":"9"},{"1":"temperature30 - temperature40","2":"Leading","3":"-7.185466","4":"2.628383","5":"130","6":"-14.026027","7":"-0.3449050","8":"-2.7337967","9":"3.549609e-02","_rn_":"10"},{"1":"temperature30 - temperature50","2":"Leading","3":"-10.671753","4":"2.684294","5":"130","6":"-17.657826","7":"-3.6856806","8":"-3.9756274","9":"6.630035e-04","_rn_":"11"},{"1":"temperature40 - temperature50","2":"Leading","3":"-3.486287","4":"2.684294","5":"130","6":"-10.472360","7":"3.4997854","8":"-1.2987725","9":"5.653739e-01","_rn_":"12"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
# {-}

# Summary figure 

![](enzyme_analysis_LDH-CS_files/figure-html/sum-fig-1.png)<!-- -->

# Conclusion 

* In conclusion the LDH:CS ratio has a **significantly** positively correlated with temperature, however, there is no significant difference in the relationship between temperature and LDH:CS when comparing fish from low- and high-latitudes.


