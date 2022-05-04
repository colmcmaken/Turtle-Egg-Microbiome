###Determine SSand Grain Size and Sorting Coefficient
library("G2Sd")
#SAND Sieve
#Import Data
SandSieveMetadata <- read.csv("~/Grad School/Thesis/Data/Sand Sieve Data 2021.csv", row.names=1)
#Run Stats
SandStats <- granstat(SandSieveMetadata,statistic="all",aggr=TRUE,modes=FALSE)
#Save as Excel file
write.csv(SandStats, "SandStats2021.csv")


###Environmental Data Statistics
#Import Dataset
Nest_Differences <- read.delim("~/Grad School/Thesis/Data/Turtle_Metadata_2021_simplified.txt", row.names=1)

#RQ: How do species and beach affect the hatching success?

#Hatching Success = Response (DV), proportion between 0 and 1
#Species = categorical predictor (IV1)
    #2-levels: CC, CM
#Beach = categorical predictor (IV2)
    #2-levels: FT, H

#Statistical Test
#Binomial GLM
#Dependent variable is set to proportions

#Null Hypotheses
#H0: Species does not affect hatching success of turtles
#H0: Beach does not affect hatching success of turtles
#H0: Species and beach do not interact

#Alternative Hypotheses
#H1: Species affects hatching success of turtles
#H2: Beach affects hatching success of turtles
#H3: Species and beach interact


#Fit the full (most complex)
#2 IVs and DV is counts and already expressed as proportion
mod.1 <- glm(Hatch.Success ~ Species * Beach,
             weights = Clutch.Size,
             family = 'binomial', data = Nest_Differences) 

#Determine if there is over-dispersion in the data
par(mfrow=c(2,2))
plot(mod.1) > mod.1$deviance / mod.1$df.residual
mod.1$deviance / mod.1$df.residual
#15.64805
#Result is >1, it suggests over-dispersion in the data, try refitting the model

#Refit data with a quasibinomial model
mod.2 <- glm(Hatch.Success ~ Species * Beach,
             weights = Clutch.Size,
             family = 'quasibinomial', data = Nest_Differences) 

#More than 1 IV, Term selection for quasibinomial model
anova(mod.2, test='Chisq')
#               Df Deviance Resid. Df Resid. Dev Pr(>Chi)  
# NULL                             26     466.60           
# Species        1   90.305        25     376.29   0.0225 *
# Beach          1   12.598        24     363.69   0.3940  
# Species:Beach  1    3.787        23     359.91   0.6403  

###Only Species is significant

mod.3 <- glm(Hatch.Success ~ Species,
             weights = Clutch.Size,
             family = 'quasibinomial', data = Nest_Differences) 
anova(mod.3, test='Chisq')
#         Df Deviance Resid. Df Resid. Dev Pr(>Chi)  
# NULL                       26     466.60           
# Species  1   90.305        25     376.29  0.02392 *

#Validate the fitted model
par(mfrow=c(2,2))
plot(mod.3)
par(mfrow=c(1,1))
plot(residuals(mod.3)~ fitted(mod.3))
plot(residuals(mod.3)~ factor(Nest_Differences$Species))

#Estimate approximate R^2
(mod.3$null.deviance-mod.3$deviance)/mod.3$null.deviance
#0.1935399

#Produce a summary plot
predictions <- (predict(mod.3, newdata = Nest_Differences, type='response'))
all.data <- cbind(Nest_Differences, predictions)
par(mfrow=c(1,1)) 
#If you have a single, categorical IV:
boxplot(predictions ~ Nest_Differences$Species)



###Descriptive Stats & T-Test for Species Differences
#Import Dataset
Nest_Differences <- read.delim("~/Grad School/Thesis/Data/Turtle_Metadata_2021_simplified.txt", row.names=1)

#####Descriptive Statistics
tapply(Nest_Differences$Hatch.Success, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.370   0.845   0.890   0.835   0.930   1.000 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9300  0.9500  0.9600  0.9543  0.9600  0.9700 

tapply(Nest_Differences$High.Tide.Distance, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   18.00   27.50   32.05   35.75  149.00 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.00   37.50   44.00   55.29   53.00  128.00 
tapply(Nest_Differences$Dune.Distance, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    13.5    31.0   197.3   363.2   613.0 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0       0       3       0      21
tapply(Nest_Differences$Incubation.Length, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 50.00   52.75   54.00   53.70   55.00   57.00 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 50.0    52.5    55.0    54.0    55.5    57.0
tapply(Nest_Differences$Clutch.Size, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 85.00   98.75  112.50  113.30  129.00  139.00 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 84     115     116     114     121     126 
tapply(Nest_Differences$Chamber.Depth, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 41.00   51.75   54.00   56.95   63.00   75.00 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 60.0    67.5    71.0    71.0    75.5    80.0 
tapply(Nest_Differences$Temperature.Average, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 29.75   31.54   32.25   32.02   32.77   33.40 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 28.85   29.00   31.05   30.58   31.62   32.90
tapply(Nest_Differences$pH.Average, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.715   7.564   7.960   7.817   8.246   8.705 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.585   6.843   7.495   7.386   7.925   8.090
tapply(Nest_Differences$Conductivity.Average, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.50   19.88   27.50  106.88   47.50 1121.00 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13.00   23.00   35.50   47.79   56.25  127.50
tapply(Nest_Differences$Sand.Grain.Size, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 176.8   369.1   423.7   421.2   460.9   580.4 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 341.5   393.5   435.7   428.0   467.7   496.8
tapply(Nest_Differences$Sorting.Coefficient, Nest_Differences$Species, summary)
# $CC
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.189   1.281   1.413   1.471   1.668   1.844 
# $CM
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.273   1.292   1.306   1.433   1.550   1.772 


#####T-Test
class(Nest_Differences$Species)
Nest_Differences$Species <-factor(Nest_Differences$Species)

attach(Nest_Differences)
levels(Species)
#"CC" "CM"

#Check Normality
boxplot(Hatch.Success~Species)
boxplot(High.Tide.Distance~Species)
boxplot(Dune.Distance~Species)
boxplot(Incubation.Length~Species)
boxplot(Clutch.Size~Species)
boxplot(Chamber.Depth~Species)
boxplot(Temperature.Average~Species)
boxplot(Conductivity.Average~Species)
boxplot(pH.Average~Species)
boxplot(Sand.Grain.Size~Species)
boxplot(Sorting.Coefficient~Species)

#Confirm using the Shapiro-Wilk test (H0: data is normal):
shapiro.test(subset(Nest_Differences,Species=="CC")$Hatch.Success)
#W = 0.72904, p-value = 9.006e-05
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences,Species=="CC")$Hatch.Success))
#W = 0.65789, p-value = 1.247e-05
shapiro.test(sqrt(subset(Nest_Differences,Species=="CC")$Hatch.Success))
#W = 0.69378, p-value = 3.277e-05
shapiro.test(subset(Nest_Differences,Species=="CC")$High.Tide.Distance)
#W = 0.65977, p-value = 1.31e-05
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences,Species=="CC")$High.Tide.Distance))
#W = NaN, p-value = NA
shapiro.test(sqrt(subset(Nest_Differences,Species=="CC")$High.Tide.Distance))
#W = 0.88889, p-value = 0.02566
shapiro.test(subset(Nest_Differences,Species=="CC")$Dune.Distance)
#W = 0.79645, p-value = 0.0007651
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences,Species=="CC")$Dune.Distance))
#W = NaN, p-value = NA
shapiro.test(sqrt(subset(Nest_Differences,Species=="CC")$Dune.Distance))
#W = 0.84215, p-value = 0.003949
shapiro.test(subset(Nest_Differences,Species=="CC")$Incubation.Length)
#W = 0.97455, p-value = 0.8464
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CM")$Incubation.Length)
#W = 0.9526, p-value = 0.7533
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CC")$Clutch.Size)
#W = 0.93067, p-value = 0.1591
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CM")$Clutch.Size)
#W = 0.73116, p-value = 0.008033
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences,Species=="CM")$Clutch.Size))
#W = 0.69202, p-value = 0.003032
shapiro.test(sqrt(subset(Nest_Differences,Species=="CM")$Clutch.Size))
#W = 0.71138, p-value = 0.004929
shapiro.test(subset(Nest_Differences,Species=="CC")$Chamber.Depth)
#W = 0.94659, p-value = 0.3182
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CM")$Chamber.Depth)
#W = 0.98149, p-value = 0.9665
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CC")$Temperature.Average)
#W = 0.92821, p-value = 0.1426
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CM")$Temperature.Average)
#W = 0.8743, p-value = 0.2023
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CC")$Conductivity.Average)
#W = 0.40588, p-value = 4.975e-08
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences,Species=="CC")$Conductivity.Average))
#W = 0.87194, p-value = 0.01271
shapiro.test(sqrt(subset(Nest_Differences,Species=="CC")$Conductivity.Average))
#W = 0.6161, p-value = 4.343e-06
shapiro.test(subset(Nest_Differences,Species=="CC")$pH.Average)
#W = 0.8495, p-value = 0.005229
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences,Species=="CC")$pH.Average))
#W = 0.83161, p-value = 0.002662
shapiro.test(sqrt(subset(Nest_Differences,Species=="CC")$pH.Average))
#W = 0.84066, p-value = 0.003732
shapiro.test(subset(Nest_Differences,Species=="CC")$Sand.Grain.Size)
#W = 0.93241, p-value = 0.1718
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CM")$Sand.Grain.Size)
#W = 0.96866, p-value = 0.8887
##Normally distributed
shapiro.test(subset(Nest_Differences,Species=="CC")$Sorting.Coefficient)
#W = 0.88455, p-value = 0.02138
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences,Species=="CC")$Sorting.Coefficient))
#W = 0.89392, p-value = 0.03176
shapiro.test(sqrt(subset(Nest_Differences,Species=="CC")$Sorting.Coefficient))
#W = 0.88955, p-value = 0.02639

#Homogeneity of Variances
bartlett.test(Incubation.Length~Beach)
#Bartlett's K-squared = 0.053113, df = 1, p-value = 0.8177
##Homogenous, parametric
bartlett.test(Chamber.Depth~Beach)
#Bartlett's K-squared = 0.23743, df = 1, p-value = 0.6261
##Homogenous, parametric
bartlett.test(Temperature.Average~Beach)
#Bartlett's K-squared = 0.20518, df = 1, p-value = 0.6506
##Homogenous, parametric
bartlett.test(Sand.Grain.Size~Beach)
#Bartlett's K-squared = 3.5405, df = 1, p-value = 0.05989
##Homogenous, parametric


#Non-parametric
#CC<CM
wilcox.test(Hatch.Success~Species,alternative='less')
#W = 16, p-value = 0.00148
wilcox.test(High.Tide.Distance~Species,alternative='less')
#W = 24, p-value = 0.005872
wilcox.test(Clutch.Size~Species,alternative='less')
#W = 71.5, p-value = 0.5441

#CC>CM
wilcox.test(Dune.Distance~Species,alternative='greater')
#W = 130, p-value = 0.0004492
wilcox.test(Conductivity.Average~Species,alternative='greater')
#W = 66, p-value = 0.5983
wilcox.test(pH.Average~Species,alternative='greater')
#W = 101, p-value = 0.04573
wilcox.test(Sorting.Coefficient~Species,alternative='greater')
#W = 72, p-value = 0.4677


#Parametric
#CC<CM
t.test(Incubation.Length ~ Species, alternative='less')
#t = -0.29769, df = 8.333, p-value = 0.3866
#95 percent confidence interval: -Inf 1.564302
#Sample estimates:
#Mean in group CC mean in group CM 
          # 53.7             54.0 
t.test(Chamber.Depth ~ Species, alternative='less')
#t = -4.3579, df = 13.855, p-value = 0.0003357
#95 percent confidence interval: -Inf -8.367254
#Sample estimates:
#Mean in group CC mean in group CM 
          # 56.95            71.00 
t.test(Sand.Grain.Size ~ Species, alternative='less')
#t = -0.23594, df = 17.181, p-value = 0.4081
#95 percent confidence interval: -Inf 43.39921
#Sample estimates:
#Mean in group CC mean in group CM 
        # 421.2102         428.0247

#CC>CM
t.test(Temperature.Average ~ Species, alternative='greater')
#t = 2.206, df = 7.7769, p-value = 0.02969
#95 percent confidence interval: 0.2226272  Inf
#Sample estimates:
#Mean in group CC mean in group CM 
        # 32.02500         30.57857 


###95% Confidence Intervals
wilcox.test(Hatch.Success ~ Species, conf.int=TRUE)
#-0.11998032 -0.02006407
wilcox.test(High.Tide.Distance ~ Species, conf.int=TRUE)
#-35.00005  -5.99996
wilcox.test(Dune.Distance ~ Species, conf.int=TRUE)
#11 361
t.test(Incubation.Length ~ Species)
#-2.607816  2.007816
t.test(Chamber.Depth ~ Species)
#-20.971684  -7.128316
wilcox.test(Clutch.Size ~ Species, conf.int=TRUE)
#-19.99998  13.00002
t.test(Temperature.Average ~ Species)
#-0.07313405  2.96599120
wilcox.test(pH.Average ~ Species, conf.int=TRUE)
#-0.09492955  1.18003369
wilcox.test(Conductivity.Average ~ Species, conf.int=TRUE)
#-29.00000  25.00004
t.test(Sand.Grain.Size ~ Species)
#-67.70227  54.07334
wilcox.test(Sorting.Coefficient ~ Species, conf.int=TRUE)
#-0.117  0.241

detach(Nest_Differences)


###T-Test for Beach Differences
#Import Dataset
Nest_Differences_CC <- droplevels(Nest_Differences[!Nest_Differences$Species == 'CM',])

#####Descriptive Statistics
tapply(Nest_Differences_CC$Hatch.Success, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5400  0.8400  0.8900  0.8573  0.9300  0.9600 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3700  0.8600  0.8900  0.8078  0.9300  1.0000 

tapply(Nest_Differences_CC$High.Tide.Distance, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   17.00   28.00   36.82   37.50  149.00 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.00   19.00   27.00   26.22   34.00   46.00
tapply(Nest_Differences_CC$Dune.Distance, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   276.0   361.0   345.8   515.5   613.0 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00   11.00   14.00   15.78   25.00   30.00
tapply(Nest_Differences_CC$Incubation.Length, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 50.00   52.00   53.00   53.27   55.00   57.00 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.00   54.00   54.00   54.22   55.00   56.00 
tapply(Nest_Differences_CC$Clutch.Size, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 85.0    99.5   113.0   112.1   127.5   134.0 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 89.0    98.0   112.0   114.8   129.0   139.0 
tapply(Nest_Differences_CC$Chamber.Depth, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 41.00   50.00   54.00   53.36   57.00   63.00 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 48.00   53.00   63.00   61.33   69.00   75.00
tapply(Nest_Differences_CC$Temperature.Average, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 29.75   30.68   32.00   31.59   32.60   32.85 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 31.70   31.95   32.50   32.56   33.05   33.40
tapply(Nest_Differences_CC$pH.Average, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.715   6.782   7.820   7.486   8.012   8.665 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.820   7.975   8.225   8.221   8.325   8.705 
tapply(Nest_Differences_CC$Conductivity.Average, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.00   27.25   41.00  175.91  165.50 1121.00 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.5    19.0    20.0    22.5    27.0    38.0 
tapply(Nest_Differences_CC$Sand.Grain.Size, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 176.8   408.2   431.0   430.1   466.6   580.4 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 344.8   368.3   399.8   410.3   443.0   538.5 
tapply(Nest_Differences_CC$Sorting.Coefficient, Nest_Differences_CC$Beach, summary)
# $`Fort Lauderdale`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.189   1.442   1.556   1.573   1.762   1.844 
# $Hillsboro
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.269   1.277   1.286   1.345   1.315   1.662 


#####T-Test
class(Nest_Differences_CC$Beach)
Nest_Differences_CC$Beach <-factor(Nest_Differences_CC$Beach)

attach(Nest_Differences_CC)
levels(Beach)
#"Fort Lauderdale" "Hillsboro" 

#Check Normality
boxplot(Hatch.Success~Beach)
boxplot(High.Tide.Distance~Beach)
boxplot(Dune.Distance~Beach)
boxplot(Incubation.Length~Beach)
boxplot(Clutch.Size~Beach)
boxplot(Chamber.Depth~Beach)
boxplot(Temperature.Average~Beach)
boxplot(Conductivity.Average~Beach)
boxplot(pH.Average~Beach)
boxplot(Sand.Grain.Size~Beach)
boxplot(Sorting.Coefficient~Beach)

#Confirm using the Shapiro-Wilk test (H0: data is normal):
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Hatch.Success)
#W = 0.76248, p-value = 0.003026
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Hatch.Success))
#W = 0.6936, p-value = 0.0003871
shapiro.test(sqrt(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Hatch.Success))
#W = 0.72829, p-value = 0.001085
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$High.Tide.Distance)
#W = 0.70007, p-value = 0.0004688
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$High.Tide.Distance))
#W = NaN, p-value = NA
shapiro.test(sqrt(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$High.Tide.Distance))
#W = 0.91206, p-value = 0.258
shapiro.test(sqrt(subset(Nest_Differences_CC,Beach=="Hillsboro")$High.Tide.Distance))
#W = 0.97019, p-value = 0.8963
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Dune.Distance)
#W = 0.92389, p-value = 0.3524
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Hillsboro")$Dune.Distance)
#W = 0.91697, p-value = 0.3678
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Incubation.Length)
#W = 0.95715, p-value = 0.7357
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Hillsboro")$Incubation.Length)
#W = 0.917, p-value = 0.3679
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Clutch.Size)
#W = 0.91514, p-value = 0.2802
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Hillsboro")$Clutch.Size)
#W = 0.93564, p-value = 0.5368
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Chamber.Depth)
#W = 0.95712, p-value = 0.7353
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Hillsboro")$Chamber.Depth)
#W = 0.92204, p-value = 0.4094
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Temperature.Average)
#W = 0.89358, p-value = 0.1541
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Hillsboro")$Temperature.Average)
#W = 0.91756, p-value = 0.3724
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Conductivity.Average)
#W = 0.53726, p-value = 4.291e-06
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Conductivity.Average))
#W = 0.92963, p-value = 0.4071
##Normally distributed
shapiro.test(log(subset(Nest_Differences_CC,Beach=="Hillsboro")$Conductivity.Average))
#W = 0.93974, p-value = 0.5792
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$pH.Average)
#W = 0.83169, p-value = 0.02457
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$pH.Average))
#W = 0.82439, p-value = 0.0197
shapiro.test(sqrt(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$pH.Average))
#W = 0.82814, p-value = 0.02207
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Sand.Grain.Size)
#W = 0.89486, p-value = 0.1598
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Hillsboro")$Sand.Grain.Size)
#W = 0.90808, p-value = 0.3027
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Fort Lauderdale")$Sorting.Coefficient)
#W = 0.92696, p-value = 0.3809
##Normally distributed
shapiro.test(subset(Nest_Differences_CC,Beach=="Hillsboro")$Sorting.Coefficient)
#W = 0.63806, p-value = 0.0002674
##Not normally distributed >> transform data
shapiro.test(log(subset(Nest_Differences_CC,Beach=="Hillsboro")$Sorting.Coefficient))
#W = 0.6542, p-value = 0.0004133
shapiro.test(sqrt(subset(Nest_Differences_CC,Beach=="Hillsboro")$Sorting.Coefficient))
#W = 0.64618, p-value = 0.000333


#Homogeneity of Variances
bartlett.test(sqrt(High.Tide.Distance)~Beach)
#Bartlett's K-squared = 4.7866, df = 1, p-value = 0.02868
##Not homogenous, non-parametric
bartlett.test(Incubation.Length~Beach)
#Bartlett's K-squared = 1.6395, df = 1, p-value = 0.2004
##Homogenous, parametric
bartlett.test(Dune.Distance~Beach)
#Bartlett's K-squared = 35.466, df = 1, p-value = 2.596e-09
##Not homogenous, non-parametric
bartlett.test(Clutch.Size~Beach)
#Bartlett's K-squared = 0.048941, df = 1, p-value = 0.8249
##Homogenous, parametric
bartlett.test(Chamber.Depth~Beach)
#Bartlett's K-squared = 1.8083, df = 1, p-value = 0.1787
##Homogenous, parametric
bartlett.test(Temperature.Average~Beach)
#Bartlett's K-squared = 2.3481, df = 1, p-value = 0.1254
##Homogenous, parametric
bartlett.test(log(Conductivity.Average)~Beach)
#Bartlett's K-squared = 7.1513, df = 1, p-value = 0.007491
##Not homogenous, non-parametric
bartlett.test(Sand.Grain.Size~Beach)
#Bartlett's K-squared = 2.1702, df = 1, p-value = 0.1407
##Homogenous, parametric


#Non-parametric
#FT<H
wilcox.test(pH.Average~Beach,alternative='less')
#W = 16.5, p-value = 0.006754

#FT>H
wilcox.test(Hatch.Success~Beach,alternative='greater')
#W = 49, p-value = 0.5303
wilcox.test(High.Tide.Distance~Beach,alternative='greater')
#W = 53.5, p-value = 0.395
wilcox.test(Dune.Distance~Beach,alternative='greater')
#W = 90, p-value = 0.0005775
wilcox.test(Conductivity.Average~Beach,alternative='greater')
#W = 82.5, p-value = 0.006754
wilcox.test(Sorting.Coefficient~Beach,alternative='greater')
#W = 76, p-value = 0.02323


#Parametric
#FT<H
t.test(Incubation.Length ~ Beach, alternative='less')
#t = -1.2557, df = 17.083, p-value = 0.1131
#95 percent confidence interval: -Inf 0.3655195
#Sample estimates:
#Mean in group Fort Lauderdale       mean in group Hillsboro 
                    # 53.27273                      54.22222
t.test(Chamber.Depth ~ Beach, alternative='less')
#t = -2.0966, df = 12.931, p-value = 0.02813
#95 percent confidence interval:-Inf -1.235349
#Sample estimates:
#Mean in group Fort Lauderdale       mean in group Hillsboro 
                    # 53.36364                      61.33333 
t.test(Temperature.Average ~ Beach, alternative='less')
#t = -2.454, df = 16.389, p-value = 0.01283
#95 percent confidence interval: -Inf -0.2822822
#Sample estimates:
#Mean in group Fort Lauderdale       mean in group Hillsboro 
                    # 31.58636                      32.56111 

#FT>H
t.test(Clutch.Size ~ Beach, alternative='greater')
#t = -0.33728, df = 16.626, p-value = 0.6299
#95 percent confidence interval: -16.56304   Inf
#Sample estimates:
#Mean in group Fort Lauderdale       mean in group Hillsboro 
                    # 112.0909                      114.7778 
t.test(Sand.Grain.Size ~ Beach, alternative='greater')
#t = 0.51297, df = 16.562, p-value = 0.3074
#95 percent confidence interval:-47.54909    Inf
#Sample estimates:
#Mean in group Fort Lauderdale       mean in group Hillsboro 
                    # 430.1391                      410.2972 


###95% Confidence Intervals using all Two-Tailed Tests
wilcox.test(Hatch.Success ~ Beach, conf.int=TRUE)
#-0.07998733  0.08995845
wilcox.test(High.Tide.Distance ~ Beach, conf.int=TRUE)
#-13.00003  18.00004
wilcox.test(Dune.Distance ~ Beach, conf.int=TRUE)
#252 504
t.test(Incubation.Length ~ Beach)
#-2.5442162  0.6452263
t.test(Chamber.Depth ~ Beach)
#-16.1860518   0.2466579
t.test(Clutch.Size ~ Beach)
#-19.52313  14.14939
t.test(Temperature.Average ~ Beach)
#-1.8151628 -0.1343321
wilcox.test(pH.Average ~ Beach, conf.int=TRUE)
#-1.44005094 -0.09006377
wilcox.test(Conductivity.Average ~ Beach, conf.int=TRUE)
#4.999986 150.500016
t.test(Sand.Grain.Size ~ Beach)
#-61.93105 101.61479
wilcox.test(Sorting.Coefficient ~ Beach, conf.int=TRUE)
#0.024 0.464

detach(Nest_Differences_CC)
                