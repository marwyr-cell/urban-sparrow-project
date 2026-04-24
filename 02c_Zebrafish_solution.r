## ##################################################### Zebrafish 15-04-2026 ##
##                                                                            ##
## Two research questions we could try to answer:                             ##
## - After controlling for differences in initial wound size, does wound size ##
##   at day 14 differ between the two treatment groups on the one hand, and   ##
##   the control group on the other?                                          ##
## - After controlling for differences in initial wound size, does wound size ##
##   at day 14 differ between the nanoparticle GF and standard GF treatment?  ##
##                                                                            ##
## ########################################################################## ##

# 00. Prepare session, import and correct the data -----------------------------

## Load packages -----
  library(pastecs)
  library(ggpubr); theme_set(theme_bw())
  library(GGally)
  library(performance)
  library(DHARMa)
  library(car)
  library(ggeffects)
  library(ggplot2)
  library(carData)
# 01. Import, screen and correct the data --------------------------------------
  fish<- read.csv("~/UZH/BME321/data/raw/Zebrafish.csv", stringsAsFactors= T)
#characters are automatically used as factors
  str(fish)
  summary(fish)

  # We could decide to reorder the levels of the treatment variable
  fish$treatment<- factor(fish$treatment,
                          levels= levels(fish$treatment)[c(1, 3, 2)])

## Set contrasts -----
## (to make sure our model actually answers our two research questions)
  contrasts(fish$treatment)<- cbind(Control.vs.GF= c(2, -1, -1),
                                    Standard.vs.NanoP= c(0, 1, -1))

# 02. Explore the data ---------------------------------------------------------

## Numerically -----

  # Let's get descriptive statistics for the three different treatment groups
  by(fish[, -1], fish$treatment,
     function(x){round(stat.desc(x, norm= T), 3)})

  # Correlation between continuous independent variable and dependent variable
  cor(fish$baseline, fish$day14)

## Graphically -----
  ggpairs(fish, aes(col= treatment, alpha= .4))

# Plot to capture the research question
  ggarrange(
    ggplot(fish, aes(treatment, day14)) +
      geom_point(aes(col= treatment), size= 3, alpha= .4,
                 position= position_jitter(height= 0, width= .3)) +
      geom_boxplot(fill= NA, outliers= F),
    ggplot(fish, aes(baseline, day14)) +
      geom_point(aes(col= treatment), size= 3, alpha= .4) +
      geom_smooth(method= lm),
    ncol= 2, nrow= 1, common.legend= T)

# 03. Fit model ----------------------------------------------------------------

## Formulate candidate models -----
  mod.00<- lm(day14~ 1, fish)
  mod.01<- lm(day14~ treatment + baseline, fish)
  mod.02<- lm(day14~ treatment * baseline, fish)

## Compare models -----
  compare_performance(mod.00, mod.01, mod.02)# Use AICc to decide: small sample!

# 04. Validate preferred model -------------------------------------------------
  check_model(mod.01)

# 05. Interpret preferred model ------------------------------------------------

## Numerically -----
  summary(mod.01)
  model_performance(mod.01)

## Graphically -----
  (pred.treat<- predict_response(mod.01, "treatment"))
  (pred.base<- predict_response(mod.01, "baseline"))

  ggarrange(
    ggplot(pred.treat, aes(x, predicted)) +
      geom_point(data= fish, aes(treatment, day14, col= treatment), size= 3,
                 alpha= .2, position= position_jitter(height= 0, width= .1)) +
      geom_errorbar(aes(ymin= conf.low, ymax= conf.high), width= .2) +
      geom_point(size= 3) +
      labs(x= "Treatment",
           y= expression("Wound size at day 14 (mm"^2*")"),
           col= "Treatment"),
    ggplot(pred.base, aes(x, predicted)) +
      geom_point(data= fish, aes(baseline, day14, col= treatment), size= 3,
                 alpha= .2) +
      geom_ribbon(aes(ymin= conf.low, ymax= conf.high), alpha= .2) +
      geom_line() +
      labs(x= expression("Wound size at day 0 (mm"^2*")"),
           y= expression("Wound size at day 14 (mm"^2*")"),
           col= "Treatment"),
    ncol= 2, common.legend= T, labels= c("a)", "b)"))

  # Fig a: difference between treatments at mean baseline wound size (145 mm2)
  # Fig b: the relationship between wound size at day 14 and baseline for the
  #        control group!

  # Alternative plot (both effects in the same graph)
  pred.fish<- predict_response(mod.01, c("baseline", "treatment"))

  ggplot(pred.fish, aes(x, predicted)) +
    geom_point(data= fish, aes(baseline, day14, col= treatment), size= 3,
               alpha= .2) +
    geom_ribbon(aes(ymin= conf.low, ymax= conf.high, fill= group), alpha= .2) +
    geom_line(aes(col= group)) +
    labs(x= expression("Wound size at day 0 (mm"^2*")"),
         y= expression("Wound size at day 14 (mm"^2*")"),
         col= "Treatment", fill= "Treatment")


# 06. Report findings ----------------------------------------------------------
  # A linear model showed that, whilst controlling for differences in wound size
  # at day 14 due to differences at day 0 (b= 0.65, SE= 0.18, t= 3.68, p< 0.001),
  # wounds in the control group were larger than in the two treatment groups
  # (b= 13.4, SE= 2.16, t= 6.21, p< 0.0001). Moreover, wounds in the standard
  # treatment group were larger than those in the nanoparticle treatment group
  # (b= 10.01, SE= 3.71, t= 2.70, p< 0.05). Overall, the model explained 60 % of
  # the variability in wound size on day 14 (R^2= 0.600, F(3, 32= 18.53,
  # p< 0.0001).

# 07. Store session information ------------------------------------------------
  version<- sessionInfo()

# 08. Appendix -----------------------------------------------------------------

## Express healing as a percentage/proportion -----

# Cool idea, takes initial wound size into account as well
  # - We might still want to use:
  #   healing % ~ Treatment + Initial condition; why?!
  # - The outcome is the ratio of two continuous variables, and can take values
  #   between 0 and 1 only -> it is bounded! In order to make sure our model
  #   doesn't make non-sensical predictions (e.g. healing of >100%, negative
  #   CIs), we need to fit a different linear model: Beta glm
  # - Since what we measure in the real world is expressed in mm^2, one could
  #   argue that if we can come up with a model on the same scale, that would be
  #   preferable

## Express healing as a difference -----

# Also interesting idea, as our model is still on mm^2 scale
  # - Again, we would still want to use:
  #   healing mm^2 ~ Treatment + Initial condition; why?!

# As long as you make explicit in your script/report WHY you "transformed" the
# outcome variable, I think there is something to say for all three approaches!

