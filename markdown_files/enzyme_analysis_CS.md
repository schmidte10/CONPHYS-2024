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

```r
cs <- read_delim("./enzymes/CS_LocalAdapt6.txt", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(...21 = col_skip(), 
                                              ...22 = col_skip()), trim_ws = TRUE) %>% 
  clean_names() %>% 
  mutate(creation_time = as.POSIXct(creation_time, format = "%d/%m/%Y %H:%M:%S"))
```

```
## New names:
## • `` -> `...21`
## • `` -> `...22`
```

```r
tissue.mass <- read.delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/enzymes/tissue_mass.txt") %>% 
  dplyr::rename(FISH_ID = fish_id)
```

# Data manipulation 

Before the data can be analysed it is important to clean up the data file. I won't explain step, that can be figured out by examining the different functions. The main steps that are occurring below are columns being broken up to make new columns (this is the _separate_ function), or columns being combined to make a unique_sample_Id value. Time data can also be tricky to deal with in R, so there are a number of data manipulation steps being used to make sure that time is being read properly.


```r
#--- data preparation/manipulation ---# 
cs2 <- cs %>%
  mutate(muscle_type = str_replace(muscle_type, " ", ".")) %>%
  unite("UNIQUE_SAMPLE_ID", c(fish_id,temperature,sample_index), sep="_", remove = FALSE) %>% 
  separate(creation_time, into=c('DATE','TIME'), sep = " ", remove = FALSE) %>% 
  arrange(sample_id_1, DATE, TIME) 

cs3 <- cs2 %>% 
  mutate(DATE = as.Date(creation_time), 
         TIME = format(creation_time, "%H:%M:%S")) %>%
  mutate(TIME = hms(cs2$TIME)) %>% 
  mutate(TIME = chron(times=cs2$TIME)) %>% 
  arrange(TIME) %>%
  group_by(UNIQUE_SAMPLE_ID, sample_id_1) %>% 
  mutate(TIME_DIFF = TIME - first(TIME)) %>% 
  filter(TIME != first(TIME)) %>%
  ungroup() %>% 
  mutate(TIME_DIFF_SECS = period_to_seconds(hms(TIME_DIFF))) %>% 
  mutate(MINUTES = TIME_DIFF_SECS/60) %>% 
  mutate(MINUTES = round(MINUTES, digits = 2)) %>% 
  dplyr::rename(CUVETTE = sample_id_1) %>% 
  mutate(REGION = substr(fish_id, 1, 1 ), 
         POPULATION = substr(fish_id, 2, 4), 
         SAMPLE_NO = substr(fish_id, 5, 7)) %>% 
  mutate(REGION = case_when( REGION =="L"~ "Leading", 
                             REGION == "C" ~ "Core", 
                             TRUE ~ "na"))   
```

## Data cleaning

Next select data points will be removed. Data points have been checked using previously written app script that allows the user to look at the plotted points of samples that were run. Because we are interested in the slope, points that plateaued were removed from the analysis, as it signifies that the reaction ran out of 'fuel' to use. For LDH 'fuel' would refer to NADH, for CS 'fuel' would refer to Oxaloacetic acid. Samples that were removed were placed into one of five different groups. No reactions reached plateau for CS, however there were a number of runs and/or cuvettes where data quailty was poor. 


```r
cs3.filtered <- cs3 %>% 
  dplyr::rename(TEMPERATURE = temperature, 
         FISH_ID = fish_id) %>%
  filter(!(TEMPERATURE == "50" & FISH_ID == "LCKM158")) %>% 
  filter(!(TEMPERATURE == "50" & FISH_ID == "CSUD010")) %>% 
  filter(!(TEMPERATURE == "40" & FISH_ID == "CSUD018")) %>% 
  filter(!(TEMPERATURE == "20" & FISH_ID == "CTON061")) %>% 
  filter(!(TEMPERATURE == "50" & FISH_ID == "CTON065")) %>%
  filter(!(FISH_ID == "CTON069")) %>% 
  filter(!(UNIQUE_SAMPLE_ID %in% c("CTON060_50_3", 
                                   "CTON061_30_3", 
                                   "CTON061_40_3", 
                                   "CTON061_50_2", 
                                   "CTON061_50_3")))%>% 
  mutate(UNIQUE_SAMPLE_ID = str_sub(UNIQUE_SAMPLE_ID, end = -3))
```

Great! Now we have all the data points that we want to keep. However, the data needs to manipulated in a way that we can obtain and pull out slopes from the absorption readings, and the calculate enzyme activity based on these slopes. This will involve a number of steps. 

## Data calculations

### Step1: Extract slopes 

Step1 will produce a data frame that provides you with the slope that was obtained for cuvettes 1-3 for each sample run at each experimental temperature

```r
CS_activity <- cs3.filtered %>% 
  dplyr::group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
  do({
    mod = lm(result ~ MINUTES, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2], 
               r2 = summary(mod)$adj.r.squared)
  }) %>%
  dplyr::ungroup() %>%
  filter(CUVETTE != ("6"))%>% 
  filter(CUVETTE != ("4"))%>% 
  filter(CUVETTE != ("5"))
```

### Step2: Slope means

Step2 will calculate the mean slope for cuvette 1-3. 

```r
CS_activity_means <- CS_activity %>%
  dplyr::group_by(UNIQUE_SAMPLE_ID) %>% 
  dplyr::mutate(Mean = mean(Slope))
```

### Step3: Background activity level

Step3 will calculate background activity level by measuring the slope from cuvette 5 (postive control)

```r
CS_background <- cs3.filtered %>% 
  group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
  do({
    mod = lm(result ~ MINUTES, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2], 
               r2 = summary(mod)$adj.r.squared)
  }) %>%
  ungroup() %>%
  filter(CUVETTE == ("5")) %>% 
  dplyr::rename(Background = Slope)
```

### Step4: Merging dataframes

Step4 will merge the data frames that you created with the mean slopes and the background slopes.

```r
final_table <- CS_activity %>% 
  full_join(distinct(CS_activity_means[,c(1,6)]), by = "UNIQUE_SAMPLE_ID") %>% 
  full_join(CS_background[,c(1,4)], by = "UNIQUE_SAMPLE_ID") 
final_table$Mean[duplicated(final_table$Mean)] <- ""
final_table$Background[duplicated(final_table$Background)] <- ""
final_table <- final_table %>% 
  mutate(Mean = as.numeric(Mean), 
         Background = as.numeric(Background), 
         Background_perc = Background/Mean) 
```

### Step5: Enzyme activity levels 

Step5 is where enzyme activity levels are calculated. See further details in manuscript (doi: xxx). Within this step background activity level is taken into account and subtracted from slopes where background activity was >5% or more of the sample slope. 

```r
CS.data <- final_table %>% 
  select(c(UNIQUE_SAMPLE_ID, Mean, Background, Background_perc)) %>% 
  mutate(Mean = as.numeric(Mean), 
         Background = as.numeric(Background), 
         Background_perc = as.numeric(Background_perc)) %>% 
  mutate(Background2 = case_when(Background_perc <= 0.05 ~ 0, 
                                 TRUE ~ Background), 
         CS_ABSORBANCE = Mean - Background2) %>%
  inner_join(select(cs3.filtered, c(UNIQUE_SAMPLE_ID, REGION, POPULATION, TEMPERATURE, FISH_ID)), by ="UNIQUE_SAMPLE_ID") %>% 
  inner_join(tissue.mass, by = "FISH_ID") %>% 
  mutate(TISSUE_MASS_CENTERED = scale(TISSUE_MASS, center = TRUE, scale = FALSE)) %>%
  distinct(UNIQUE_SAMPLE_ID, REGION, POPULATION, .keep_all = TRUE) %>% 
  mutate(REGION = factor(REGION),
         PATH_LENGTH = 1, 
         EXTINCTION_COEFFICIENT = 13.6, 
         TISSUE_CONCENTRATION = 0.01, 
         ASSAY_VOL = 0.930, 
         SAMPLE_VOL = 0.020, 
         CS_ACTIVITY = ((CS_ABSORBANCE/(PATH_LENGTH*EXTINCTION_COEFFICIENT*TISSUE_CONCENTRATION))*(ASSAY_VOL/SAMPLE_VOL)))
```

By the end of this stage you should have a data frame that included a column called **LDH_ACTIVITY** along with necessary metadata - this data frame will be used to perform the statistical analysis. 

# Exploratory data analysis {.tabset}

## LDH v TEMPERATURE [LATITUDE]

```r
ggplot(CS.data, aes(x =as.numeric(TEMPERATURE), y= CS_ACTIVITY, color = REGION)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)
```

![](enzyme_analysis_CS_files/figure-html/eda-1-1.png)<!-- -->

## LDH V TEMPERATURE [DENSITY]

```r
ggplot(CS.data, aes(x = CS_ACTIVITY, fill = TEMPERATURE, color = TEMPERATURE)) + 
  geom_density(alpha =0.5, position = "identity") 
```

![](enzyme_analysis_CS_files/figure-html/eda-2-1.png)<!-- -->

## LDH v TISSUE MASS (LATITUDE)

```r
ggplot(CS.data, aes(x =TISSUE_MASS_CENTERED, y= CS_ACTIVITY, color = REGION)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)
```

![](enzyme_analysis_CS_files/figure-html/eda-3-1.png)<!-- -->


## {-}

# Fit the model 

The model was fit using the **glm** and later **glmmTMB** package in R. A number of different models were tested to determine which hypothesis and associated variables best predicted resting oxygen consumption. Model fit was examined using AICc, BIC, and r-squared values. Additional model were examined via the validation diagonistics provided by the **performance** and **dHARMA** packages in R. 

## Fixed factors (linear regression models)

### model 1

```r
#--- base model ---#
cs.model.1 <- glm(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED, 
                       family=gaussian(), 
                       data = CS.data)  
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
   <td style="text-align:right;"> -1.3911899 </td>
   <td style="text-align:right;"> 0.5554433 </td>
   <td style="text-align:right;"> -2.5046481 </td>
   <td style="text-align:right;"> 0.0135037 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -0.7102694 </td>
   <td style="text-align:right;"> 0.8088557 </td>
   <td style="text-align:right;"> -0.8781164 </td>
   <td style="text-align:right;"> 0.3815132 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TEMPERATURE </td>
   <td style="text-align:right;"> 0.1582013 </td>
   <td style="text-align:right;"> 0.0153213 </td>
   <td style="text-align:right;"> 10.3255866 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TISSUE_MASS_CENTERED </td>
   <td style="text-align:right;"> -0.0552681 </td>
   <td style="text-align:right;"> 0.0247763 </td>
   <td style="text-align:right;"> -2.2306875 </td>
   <td style="text-align:right;"> 0.0274296 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:TEMPERATURE </td>
   <td style="text-align:right;"> 0.0281609 </td>
   <td style="text-align:right;"> 0.0222012 </td>
   <td style="text-align:right;"> 1.2684382 </td>
   <td style="text-align:right;"> 0.2069258 </td>
  </tr>
</tbody>
</table>

### model 2

```r
cs.model.2 <- glm(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE, 
                       family=gaussian(), 
                       data = CS.data) 
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
   <td style="text-align:right;"> -1.4649587 </td>
   <td style="text-align:right;"> 0.5628730 </td>
   <td style="text-align:right;"> -2.6026454 </td>
   <td style="text-align:right;"> 0.0103250 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -0.5879236 </td>
   <td style="text-align:right;"> 0.8192418 </td>
   <td style="text-align:right;"> -0.7176435 </td>
   <td style="text-align:right;"> 0.4742647 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TEMPERATURE </td>
   <td style="text-align:right;"> 0.1594114 </td>
   <td style="text-align:right;"> 0.0155441 </td>
   <td style="text-align:right;"> 10.2554460 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:TEMPERATURE </td>
   <td style="text-align:right;"> 0.0269880 </td>
   <td style="text-align:right;"> 0.0225318 </td>
   <td style="text-align:right;"> 1.1977743 </td>
   <td style="text-align:right;"> 0.2331842 </td>
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
   <td style="text-align:left;"> cs.model.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 481.1473 </td>
   <td style="text-align:right;"> 497.8730 </td>
   <td style="text-align:right;"> 0.6531851 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs.model.2 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 484.0263 </td>
   <td style="text-align:right;"> 498.0468 </td>
   <td style="text-align:right;"> 0.6417294 </td>
  </tr>
</tbody>
</table>

The model that contains **TISSUE_MASS_CENTERED** seems to do better than the model that leaves TISSUE_MASS_CENTERED out. Therefore we will move ahead with the model that contains **TISSUE_MASS_CENTERED** as a co-variate.  

## Polynomials 

### polynomial models 

Note that the linear model has already been created via model _ldh.model.1_ in the previous section.


```r
#--- second order polynomial ---# 
cs.model.1.p2 <- glm(CS_ACTIVITY ~ 1 + REGION*poly(TEMPERATURE, 2) + TISSUE_MASS_CENTERED, 
                      family=gaussian(), 
                      data = CS.data)  

#--- third order polynomial ---# 
cs.model.1.p3 <- glm(CS_ACTIVITY ~ 1 + REGION*poly(TEMPERATURE, 3) + TISSUE_MASS_CENTERED, 
                      family=gaussian(), 
                      data = CS.data)  
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
   <td style="text-align:left;"> cs.model.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 481.1473 </td>
   <td style="text-align:right;"> 497.8730 </td>
   <td style="text-align:right;"> 0.6600699 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs.model.1.p2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 484.7950 </td>
   <td style="text-align:right;"> 506.8257 </td>
   <td style="text-align:right;"> 0.6622016 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs.model.1.p3 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 489.1785 </td>
   <td style="text-align:right;"> 516.3682 </td>
   <td style="text-align:right;"> 0.6628391 </td>
  </tr>
</tbody>
</table>

From our model comparison we can see that the model that runs temperature as a linear model performs the best. Therefore, moving forward we will use the linear model. 

## Random factors 

Fish were repeatedly sampled over four different temperatures, therefore repeated sampling needs to be accounted for. To do this random factors will be included within the model. There are a number of options that can be used for random factors including 1) accounting for repeated sampling of individuals, 2) accounting for repeated sampling of individuals nested within population, 3) account for repeated sampling of individuals and populations without nesting. All three models will be run a compaired. 

### random factor models


```r
cs.model.1a <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                       family=gaussian(), 
                       data = CS.data, 
                       REML = TRUE) 

cs.model.1b <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED + (1|POPULATION/FISH_ID), 
                       family=gaussian(), 
                       data = CS.data,
                       REML = TRUE) 

cs.model.1c <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED + (1|FISH_ID) + (1 + REGION|POPULATION), 
                       family=gaussian(), 
                       data = CS.data,
                       REML = TRUE)
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
   <td style="text-align:left;"> cs.model.1a </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 415.8362 </td>
   <td style="text-align:right;"> 435.2322 </td>
   <td style="text-align:right;"> 0.6471845 </td>
   <td style="text-align:right;"> 0.6471845 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs.model.1b </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 416.4375 </td>
   <td style="text-align:right;"> 438.4682 </td>
   <td style="text-align:right;"> 0.6304244 </td>
   <td style="text-align:right;"> 0.6304244 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs.model.1c </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 420.9443 </td>
   <td style="text-align:right;"> 448.1341 </td>
   <td style="text-align:right;"> 0.6305698 </td>
   <td style="text-align:right;"> 0.6305698 </td>
  </tr>
</tbody>
</table>

Model _ldh.model.1a_ appears to be the best model, however, there seems to be little difference in how the models change depending on how the random factors are arranged.

# Model validation {.tabset .tabset-faded}

## performance {.tabset .tabset-faded}

### rmr.3a (linear)
![](enzyme_analysis_CS_files/figure-html/model-valid-1-1.png)<!-- -->

## DHARMa residuals {.tabset .tabset-faded}

### nas.1a (linear)

```r
cs.model.1a %>% simulateResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-2-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.628 0.612 0.56 0.52 0.22 0.188 0.092 0.088 0.276 0.348 0.576 0.684 0.504 0.36 0.056 0.552 0.736 0.592 0.24 0.18 ...
```

```r
cs.model.1a %>% DHARMa::testResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-2-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.10788, p-value = 0.08839
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.97124, p-value = 0.92
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 134, p-value = 0.2892
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001812671 0.052874653
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01492537
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.10788, p-value = 0.08839
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.97124, p-value = 0.92
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 134, p-value = 0.2892
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001812671 0.052874653
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01492537
```

## {-}

# {-}

The _cs.model.1a_ model looks okay....lets play around with a link transformation to see if we can get any improvement 

# Fit the model (link transformations)


```r
cs.model.1a <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED + (1|FISH_ID), #deafult
                       family=gaussian(link = "identity"), 
                       data = CS.data, 
                       REML = TRUE) 

cs.model.1a.log <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                       family=gaussian(link="log"), 
                       data = CS.data, 
                       REML = TRUE)  

cs.model.1a.inv <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                       family=gaussian(link="inverse"), 
                       data = CS.data, 
                       REML = TRUE) 
```

# Model re-validation {.tabset .tabset-faded}

## performance {.tabset .tabset-faded}

### Gaussian (identity)

![](enzyme_analysis_CS_files/figure-html/model-valid-2.2a-1.png)<!-- -->

### Gaussian (log)
![](enzyme_analysis_CS_files/figure-html/model-valid-2.2b-1.png)<!-- -->

### Gaussian (inverse)
![](enzyme_analysis_CS_files/figure-html/model-valid-2.2c-1.png)<!-- -->

## DHARMa {.tabset .tabset-faded}


```r
cs.model.1a %>% simulateResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-2.3a-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.628 0.612 0.56 0.52 0.22 0.188 0.092 0.088 0.276 0.348 0.576 0.684 0.504 0.36 0.056 0.552 0.736 0.592 0.24 0.18 ...
```

```r
cs.model.1a %>% DHARMa::testResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-2.3a-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.10788, p-value = 0.08839
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.97124, p-value = 0.92
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 134, p-value = 0.2892
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001812671 0.052874653
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01492537
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.10788, p-value = 0.08839
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.97124, p-value = 0.92
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 2, observations = 134, p-value = 0.2892
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.001812671 0.052874653
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                             0.01492537
```


```r
cs.model.1a.log %>% simulateResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-2.3b-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.6 0.744 0.68 0.516 0.012 0.136 0.124 0.132 0.056 0.4 0.704 0.636 0.576 0.464 0.076 0.28 0.856 0.588 0.016 0.112 ...
```

```r
cs.model.1a.log %>% DHARMa::testResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-2.3b-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.085672, p-value = 0.279
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 1.0782, p-value = 0.68
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 1, observations = 134, p-value = 1
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.000188921 0.040877038
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                            0.007462687
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.085672, p-value = 0.279
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 1.0782, p-value = 0.68
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 1, observations = 134, p-value = 1
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.000188921 0.040877038
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                            0.007462687
```


```r
cs.model.1a.inv %>% simulateResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-2.3c-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.356 0.692 0.716 0.548 0.008 0.104 0.124 0.124 0.032 0.34 0.784 0.38 0.492 0.488 0.092 0.156 0.756 0.576 0.004 0.08 ...
```

```r
cs.model.1a.inv %>% DHARMa::testResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-2.3c-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.14233, p-value = 0.008775
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.10791, p-value = 0.8
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 0, observations = 134, p-value = 0.6316
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.00000000 0.02715348
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                                      0
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.14233, p-value = 0.008775
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 0.10791, p-value = 0.8
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 0, observations = 134, p-value = 0.6316
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.00000000 0.02715348
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                                      0
```
## {-}

# {-}

From looking at the different models it looks like the model with the log-link function performs the best. In the DHARMa validation test we can see that one of quantile deviations is violated. Because the model passes all the other data validations realtively well we could move on with the log-link model. However, previously we showed that the 2^nd^ and 3^rd^ order polynomials also performed quite well, and we know the LDH model was not linear. So before we choose out final model, lets see what the 2^nd^ and 3^rd^ order polynomials look like with a log-link. 

# Fit model - polynomials and link functions 


```r
cs.model.1a.log <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                       family=gaussian(link="log"), 
                       data = CS.data, 
                       REML = TRUE) 

cs.model.1a.log.p2 <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*poly(TEMPERATURE, 2) + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                       family=gaussian(link="log"), 
                       data = CS.data, 
                       REML = TRUE)

cs.model.1a.log.p3 <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*poly(TEMPERATURE, 3) + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                       family=gaussian(link="log"), 
                       data = CS.data, 
                       REML = TRUE)
```

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
   <td style="text-align:left;"> cs.model.1a.log </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 381.7131 </td>
   <td style="text-align:right;"> 401.1091 </td>
   <td style="text-align:right;"> 0.3191810 </td>
   <td style="text-align:right;"> 0.3191810 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs.model.1a.log.p2 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 336.8235 </td>
   <td style="text-align:right;"> 361.4524 </td>
   <td style="text-align:right;"> 0.4330873 </td>
   <td style="text-align:right;"> 0.4330873 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs.model.1a.log.p3 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 343.1807 </td>
   <td style="text-align:right;"> 372.8930 </td>
   <td style="text-align:right;"> 0.4330877 </td>
   <td style="text-align:right;"> 0.4330877 </td>
  </tr>
</tbody>
</table>

From this model comparison we can see that the 2^nd^ order polynomial with the log-link seems to be the best model. Let's look at our model validations 

# Model re-re-validation 

![](enzyme_analysis_CS_files/figure-html/model-valid-3.2a-1.png)<!-- -->

## DHARMa 


```r
cs.model.1a.log.p2 %>% simulateResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-3.2b-1.png)<!-- -->

```
## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
##  
## Scaled residual values: 0.756 0.704 0.592 0.564 0.024 0.096 0.06 0.152 0.108 0.36 0.64 0.796 0.548 0.388 0.088 0.4 0.844 0.652 0.032 0.096 ...
```

```r
cs.model.1a.log.p2 %>% DHARMa::testResiduals(plot=TRUE)
```

![](enzyme_analysis_CS_files/figure-html/model-valid-3.2b-2.png)<!-- -->

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.071522, p-value = 0.4994
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 1.0848, p-value = 0.696
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 1, observations = 134, p-value = 1
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.000188921 0.040877038
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                            0.007462687
```

```
## $uniformity
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.071522, p-value = 0.4994
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs.
## 	simulated
## 
## data:  simulationOutput
## dispersion = 1.0848, p-value = 0.696
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
## 	DHARMa outlier test based on exact binomial test with approximate
## 	expectations
## 
## data:  simulationOutput
## outliers at both margin(s) = 1, observations = 134, p-value = 1
## alternative hypothesis: true probability of success is not equal to 0.007968127
## 95 percent confidence interval:
##  0.000188921 0.040877038
## sample estimates:
## frequency of outliers (expected: 0.00796812749003984 ) 
##                                            0.007462687
```

Validations look great! Moving ahead with the quadratic log-link model. 

# Partial plots {.tabset .tabset-faded}

## ggemmeans 

![](enzyme_analysis_CS_files/figure-html/partial-plots-1-1.png)<!-- -->

## plot_model 

![](enzyme_analysis_CS_files/figure-html/partial-plots-2-1.png)<!-- -->

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
   <td style="text-align:right;"> 1.2606077 </td>
   <td style="text-align:right;"> 0.0644433 </td>
   <td style="text-align:right;"> 19.5614886 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> 0.0672739 </td>
   <td style="text-align:right;"> 0.0956223 </td>
   <td style="text-align:right;"> 0.7035383 </td>
   <td style="text-align:right;"> 0.4817204 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> 5.3478569 </td>
   <td style="text-align:right;"> 0.2433255 </td>
   <td style="text-align:right;"> 21.9781996 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -0.8924033 </td>
   <td style="text-align:right;"> 0.1895792 </td>
   <td style="text-align:right;"> -4.7072860 </td>
   <td style="text-align:right;"> 0.0000025 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TISSUE_MASS_CENTERED </td>
   <td style="text-align:right;"> -0.0081022 </td>
   <td style="text-align:right;"> 0.0092873 </td>
   <td style="text-align:right;"> -0.8723992 </td>
   <td style="text-align:right;"> 0.3829906 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> 0.5404937 </td>
   <td style="text-align:right;"> 0.3578765 </td>
   <td style="text-align:right;"> 1.5102800 </td>
   <td style="text-align:right;"> 0.1309720 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> 0.1512415 </td>
   <td style="text-align:right;"> 0.2733172 </td>
   <td style="text-align:right;"> 0.5533553 </td>
   <td style="text-align:right;"> 0.5800201 </td>
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
   <td style="text-align:right;"> 1.2058605 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.2721535 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2) </td>
   <td style="text-align:right;"> 1312.4238850 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TISSUE_MASS_CENTERED </td>
   <td style="text-align:right;"> 0.7610803 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.3829906 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGION:poly(TEMPERATURE, 2) </td>
   <td style="text-align:right;"> 6.0796922 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0478423 </td>
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
   <td style="text-align:right;"> 1.1343011 </td>
   <td style="text-align:right;"> 1.3869144 </td>
   <td style="text-align:right;"> 1.2606077 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading </td>
   <td style="text-align:right;"> -0.1201423 </td>
   <td style="text-align:right;"> 0.2546901 </td>
   <td style="text-align:right;"> 0.0672739 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> 4.8709476 </td>
   <td style="text-align:right;"> 5.8247661 </td>
   <td style="text-align:right;"> 5.3478569 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -1.2639717 </td>
   <td style="text-align:right;"> -0.5208350 </td>
   <td style="text-align:right;"> -0.8924033 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TISSUE_MASS_CENTERED </td>
   <td style="text-align:right;"> -0.0263050 </td>
   <td style="text-align:right;"> 0.0101005 </td>
   <td style="text-align:right;"> -0.0081022 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)1 </td>
   <td style="text-align:right;"> -0.1609313 </td>
   <td style="text-align:right;"> 1.2419187 </td>
   <td style="text-align:right;"> 0.5404937 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REGIONLeading:poly(TEMPERATURE, 2)2 </td>
   <td style="text-align:right;"> -0.3844504 </td>
   <td style="text-align:right;"> 0.6869335 </td>
   <td style="text-align:right;"> 0.1512415 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Std.Dev.(Intercept)|FISH_ID </td>
   <td style="text-align:right;"> 0.2053835 </td>
   <td style="text-align:right;"> 0.3437883 </td>
   <td style="text-align:right;"> 0.2657225 </td>
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
   <td style="text-align:right;"> 0.5570967 </td>
   <td style="text-align:right;"> 0.4330873 </td>
   <td style="text-align:left;"> FALSE </td>
  </tr>
</tbody>
</table>

# {-} 

# Pairwise comparisons {.tabset .tabset-faded} 

## emtrends [latitudes]



```r
cs.model.1a.log.p2  %>% emtrends(var = "TEMPERATURE", type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["estimate"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[9],"type":["dbl"],"align":["right"]}],"data":[{"1":"(Core TISSUE_MASS_CENTERED-0.117338329414431) - (Leading TISSUE_MASS_CENTERED-0.117338329414431)","2":"34.62687","3":"-0.00416897","4":"0.002870779","5":"132","6":"-0.009847655","7":"0.001509715","8":"-1.452209","9":"0.1488165","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
SCROLL TO THE RIGHT -->

The numbers in the left most column in the table just mention that the slopes are assuming mean **TISSUE_MASS_CENTERED** values when looking at differences between latitudinal slopes.

## emmeans [latitudes]

```r
cs.model.1a.log.p2  %>% emmeans(pairwise ~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["fct"],"align":["left"]},{"label":["ratio"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["null"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[10],"type":["dbl"],"align":["right"]}],"data":[{"1":"Core / Leading","2":"34.6268656716418","3":"0.9499583","4":"0.091741","5":"132","6":"0.7847662","7":"1.149923","8":"1","9":"-0.5315852","10":"0.5959068","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## TEMPERATURE 

```r
cs.model.1a.log.p2  %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["TEMPERATURE"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["REGION"],"name":[2],"type":["fct"],"align":["left"]},{"label":["response"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["null"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[10],"type":["dbl"],"align":["right"]}],"data":[{"1":"34.62687","2":"Core","3":"3.879063","4":"0.2534753","5":"132","6":"3.408716","7":"4.414310","8":"1","9":"20.74535","10":"2.261412e-43","_rn_":"1"},{"1":"34.62687","2":"Leading","3":"4.083403","4":"0.2883748","5":"132","6":"3.551021","7":"4.695603","8":"1","9":"19.92222","10":"1.311044e-41","_rn_":"2"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


## Means - f(TEMPERATURE)

```r
cs.model.1a.log.p2  %>% update(.~1+ REGION * as.factor(TEMPERATURE) + TISSUE_MASS_CENTERED + (1|FISH_ID)) %>% 
  emmeans(~REGION*TEMPERATURE, type = "response") %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["REGION"],"name":[1],"type":["fct"],"align":["left"]},{"label":["TEMPERATURE"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["response"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["null"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[10],"type":["dbl"],"align":["right"]}],"data":[{"1":"Core","2":"20","3":"1.729170","4":"0.1518268","5":"132","6":"1.453476","7":"2.057157","8":"1","9":"6.237139","10":"5.602971e-09","_rn_":"1"},{"1":"Leading","2":"20","3":"1.844954","4":"0.1716901","5":"132","6":"1.534760","7":"2.217843","8":"1","9":"6.581338","10":"1.002373e-09","_rn_":"2"},{"1":"Core","2":"30","3":"3.210342","4":"0.2244893","5":"132","6":"2.795624","7":"3.686582","8":"1","9":"16.679954","10":"2.644407e-34","_rn_":"3"},{"1":"Leading","2":"30","3":"3.185232","4":"0.2445628","5":"132","6":"2.736408","7":"3.707671","8":"1","9":"15.088848","10":"1.607698e-30","_rn_":"4"},{"1":"Core","2":"40","3":"4.708608","4":"0.3108166","5":"132","6":"4.132231","7":"5.365380","8":"1","9":"23.471980","10":"5.964682e-49","_rn_":"5"},{"1":"Leading","2":"40","3":"5.184144","4":"0.3675634","5":"132","6":"4.505751","7":"5.964677","8":"1","9":"23.209744","10":"1.972310e-48","_rn_":"6"},{"1":"Core","2":"50","3":"6.202633","4":"0.4019846","5":"132","6":"5.456327","7":"7.051017","8":"1","9":"28.159393","10":"1.168851e-57","_rn_":"7"},{"1":"Leading","2":"50","3":"7.150521","4":"0.4960251","5":"132","6":"6.233676","7":"8.202214","8":"1","9":"28.358238","10":"5.262215e-58","_rn_":"8"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## Abs. diff - f(TEMPERATURE)

```r
cs.model.1a.log.p2  %>% update(.~1+ REGION * as.factor(TEMPERATURE) + TISSUE_MASS_CENTERED + (1|FISH_ID)) %>% 
  emmeans(~REGION*TEMPERATURE, type = "response") %>% pairs(by ="REGION") %>% summary(infer=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["contrast"],"name":[1],"type":["fct"],"align":["left"]},{"label":["REGION"],"name":[2],"type":["fct"],"align":["left"]},{"label":["ratio"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["SE"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["df"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["lower.CL"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["upper.CL"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["null"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["t.ratio"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[10],"type":["dbl"],"align":["right"]}],"data":[{"1":"TEMPERATURE20 / TEMPERATURE30","2":"Core","3":"0.5386247","4":"0.03811747","5":"132","6":"0.4480371","7":"0.6475280","8":"1","9":"-8.743147","10":"6.894485e-14","_rn_":"1"},{"1":"TEMPERATURE20 / TEMPERATURE40","2":"Core","3":"0.3672359","4":"0.02454860","5":"132","6":"0.3086057","7":"0.4370048","8":"1","9":"-14.985739","10":"0.000000e+00","_rn_":"2"},{"1":"TEMPERATURE20 / TEMPERATURE50","2":"Core","3":"0.2787799","4":"0.01831343","5":"132","6":"0.2349775","7":"0.3307476","8":"1","9":"-19.444461","10":"0.000000e+00","_rn_":"3"},{"1":"TEMPERATURE30 / TEMPERATURE40","2":"Core","3":"0.6818028","4":"0.02762984","5":"132","6":"0.6135691","7":"0.7576247","8":"1","9":"-9.451396","10":"1.487699e-14","_rn_":"4"},{"1":"TEMPERATURE30 / TEMPERATURE50","2":"Core","3":"0.5175773","4":"0.01990454","5":"132","6":"0.4682916","7":"0.5720502","8":"1","9":"-17.125471","10":"0.000000e+00","_rn_":"5"},{"1":"TEMPERATURE40 / TEMPERATURE50","2":"Core","3":"0.7591305","4":"0.02312775","5":"132","6":"0.7012742","7":"0.8217599","8":"1","9":"-9.045515","10":"2.320366e-14","_rn_":"6"},{"1":"TEMPERATURE20 / TEMPERATURE30","2":"Leading","3":"0.5792214","4":"0.04324379","5":"132","6":"0.4769535","7":"0.7034176","8":"1","9":"-7.314245","10":"1.344181e-10","_rn_":"7"},{"1":"TEMPERATURE20 / TEMPERATURE40","2":"Leading","3":"0.3558841","4":"0.02441492","5":"132","6":"0.2977025","7":"0.4254365","8":"1","9":"-15.059715","10":"0.000000e+00","_rn_":"8"},{"1":"TEMPERATURE20 / TEMPERATURE50","2":"Leading","3":"0.2580168","4":"0.01729113","5":"132","6":"0.2167287","7":"0.3071705","8":"1","9":"-20.215181","10":"0.000000e+00","_rn_":"9"},{"1":"TEMPERATURE30 / TEMPERATURE40","2":"Leading","3":"0.6144180","4":"0.02698986","5":"132","6":"0.5480538","7":"0.6888183","8":"1","9":"-11.088261","10":"0.000000e+00","_rn_":"10"},{"1":"TEMPERATURE30 / TEMPERATURE50","2":"Leading","3":"0.4454545","4":"0.01843643","5":"132","6":"0.3999747","7":"0.4961057","8":"1","9":"-19.538566","10":"0.000000e+00","_rn_":"11"},{"1":"TEMPERATURE40 / TEMPERATURE50","2":"Leading","3":"0.7250023","4":"0.02098010","5":"132","6":"0.6724155","7":"0.7817017","8":"1","9":"-11.112749","10":"0.000000e+00","_rn_":"12"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
# {-}

# Summary figure 

![](enzyme_analysis_CS_files/figure-html/sum-fig-1.png)<!-- -->

# Conclusion 

* In conclusion CS enzyme activity has a **significantly** positively correlated with temperature, however, there is no significant difference in the relationship between temperature and CS activity when comparing fish from low- and high-latitudes.


