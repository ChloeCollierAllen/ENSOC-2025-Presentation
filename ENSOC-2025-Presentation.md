ENSOC 2025 Presentation
================
Chloe Collier-Allen
2025-10-28

- [Set Up](#set-up)
- [Linear Model](#linear-model)
- [GLMs](#glms)
- [ANOVA](#anova)
- [Continuation of frequentest stats with
  interactions](#continuation-of-frequentest-stats-with-interactions)
- [Single-group Bayesian SEM](#single-group-bayesian-sem)
- [Multi-group Bayesian SEM (for
  2023–2025)](#multi-group-bayesian-sem-for-20232025)

## Set Up

    #Set working directory
    setwd("/Users/chloecollier/Desktop")

    #Read in the file we will be working with
    library(readr)
    CDC25<- read_csv("Combined-Data-Cat.csv")
    View(CDC25)

    na.omit(CDC25)

    summary(CDC25)


    #Packages
    library(DHARMa)
    library(car)
    library(lme4)
    library(lmerTest)
    library(ggplot2)
    library(cowplot)

    #Short name variables 
    site<-CDC25$'site'
    elev<-CDC25$`elev (ft)`
    date<-CDC25$date
    lat<-CDC25$lat
    long<-CDC25$long
    dbh<-CDC25$`dbh(m)`
    cone<-CDC25$`Cone Count`
    cat<-CDC25$`Caterpillar Count`
    dprec<-CDC25$`Day Prep (cm)`
    mprec<-CDC25$`Monthly Prep (cm)`
    burn<-CDC25$`Year Burned `
    nfplant<-CDC25$`Non-focal tree plants present `
    winter<-CDC25$`Winter Prior Prep (cm)`

## Linear Model

    #Cats and monthly precipitation 
    lm1<-lm(cat~winter, data=CDC25)

    plot(lm1)
    #data is non-normal (right skewed) based on Q-Q plot thus we need to use a GLM 

    #Cats and non-focal plants
    lm2<-lm(cat~nfplant, data=CDC25)

    plot(lm2)
    #data is non-normal (right skewed) based on Q-Q plot thus we need to use a GLM 

    #Cats and cone count 
    lm3<-lm(cat~cone, data=CDC25)

    plot(lm3)
    #data is non-normal (right skewed) based on Q-Q plot thus we need to use a GLM 

    #Cats and dbh 
    lm4<-lm(cat~dbh, data=CDC25)

    plot(lm4)
    #data is non-normal (right skewed) based on Q-Q plot thus we need to use a GLM 

\##Check Model Assumptions

    #Check AIC-lm1
    pineP<-glm(lm1,family=poisson)
    AIC(pineP) #445.8215

    pineG<-glm(lm1, family=gaussian)
    AIC(pineG)#468.7801

    pineGa<-glm(lm1, family=Gamma)
    AIC(pineGa)
    #Gamma cannot even be used due to non-positive values 
    #Thus, poisson is best for lm1

    #Check AIC-lm2
    pineP2<-glm(lm2,family=poisson)
    AIC(pineP2) #425.3668

    pineG2<-glm(lm2, family=gaussian)
    AIC(pineG2)#463.6354

    pineGa2<-glm(lm2, family=Gamma)
    AIC(pineGa2)
    #Gamma cannot even be used due to non-positive values 
    #Thus, poisson is best for lm2

    #Check AIC-lm3
    pineP3<-glm(lm3,family=poisson)
    AIC(pineP3) #482.7881

    pineG3<-glm(lm3, family=gaussian)
    AIC(pineG3)#463.1583

    pineGa3<-glm(lm3, family=Gamma)
    AIC(pineGa3)
    #Gamma cannot even be used due to non-positive values 
    #Thus, gaussian is best for lm3

    #Check AIC-lm4
    pineP4<-glm(lm4,family=poisson)
    AIC(pineP4) #457.1442

    pineG4<-glm(lm4, family=gaussian)
    AIC(pineG4)#466.1625

    pineGa4<-glm(lm4, family=Gamma)
    AIC(pineGa4)
    #Gamma cannot even be used due to non-positive values 
    #Thus, poisson is best for lm4

## GLMs

    length(cat)  # Length of the response variable--> 83
    length(winter)  # Length of the predictor variable--> 83

    sum(is.na(cat))  # Number of NAs in response variable --> 2
    sum(is.na(winter))  # Number of NAs in predictor variable--> 2


    data <- na.omit(data.frame( winter))
    glm1 <- glm(cat ~ winter, family = poisson, data = CDC25)


    #GLM models for lm1-4 using Poisson 

    glm1<-glm(cat~winter, family=poisson)

    plot(glm1)

    summary(glm1)
    #Caterpillar count and annual winter precipitation are significant but right skewed

    glm2<-glm(cat~nfplant, family=poisson)

    plot(glm2)

    summary(glm2)
    #Caterpillar count and non-focal plant count are significant but right skewed

    glm3<-glm(cat~cone, family=gaussian)

    plot(glm3)

    summary(glm3)
    #Caterpillar count and cone count are not significant but right skewed

    glm4<-glm(cat~dbh, family=gaussian)

    plot(glm4)

    summary(glm4)

    #Caterpillar count and dbh are slightly significant but right skewed

## ANOVA

    bcaov<- aov(cat ~ burn,
                  data = CDC25)

    hist(bcaov$residuals)

    # QQ-plot
    library(car)
    qqPlot(bcaov$residuals,
           id = FALSE # id = FALSE to remove point identification
    )
    #Normality not met; still right skewed

    shapiro.test(bcaov$residuals)
    #P-value of the Shapiro-Wilk test on the residuals is smaller than the usual significance 
    #level of α=5%, thus we do reject the hypothesis that residuals follow a normal 
    #distribution

    # Boxplot to visualize equality of the variance 
    boxplot(cat ~ burn,
            data = CDC25)

    #Test 
    # Levene's test
    library(car)

    leveneTest(cat ~ burn,
               data = CDC25)
    #The p-value being larger than the significance level of 0.05 thus, we do not reject the 
    #null hypothesis, so we cannot reject the hypothesis that variances are equal between 
    #burns (p-value = 0.3629).

    #One way ANOVA as we do not assume equal variance between burns 
    oneway.test(cat ~ burn,
                data = CDC25,
                var.equal = FALSE) # assuming unequal variances

    #A higher F-value indicates a greater difference between group means relative to the variation 
    #within the groups.
    #P-value=0.01354 means we cab reject the null hypothesis, which states that all population means 
    #are equal.Therefore, the data provides sufficient evidence to conclude that at least one 
    #of the group means is significantly different from the others. 

    #Lets look at this through a Tukey 
    library(multcomp)

    # Tukey HSD test:
    tukey.test <- TukeyHSD(bcaov)
    tukey.test
    plot(tukey.test)
    #The confidence intervals does cross the zero line, which indicate that all groups 
    #are significantly different.


    #I am having an issue with a duplication of the no-burn category thus I need to fix this:
    # ---- 0) Packages ----
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(forcats)

    # ---- 1) Clean the factor labels (fix the duplicate 'No burn ') ----
    CDC25 <- CDC25 %>%
      mutate(
        burn = burn %>%
          str_trim() %>%     # remove leading/trailing spaces
          str_squish() %>%   # collapse any double spaces
          as.factor()
      )

    # sanity checks
    levels(CDC25$burn)
    table(CDC25$burn, useNA = "ifany")

    # ---- 2) Refit ANOVA and Tukey ----
    fit <- aov(cat ~ burn, data = CDC25)
    tuk  <- TukeyHSD(fit, "burn", conf.level = 0.95)

    # ---- 3) Tidy Tukey results and mark significance ----
    tdf <- as.data.frame(tuk$burn)
    tdf$contrast <- rownames(tdf)
    tdf <- tdf %>%
      rename(diff = diff, lwr = lwr, upr = upr, padj = `p adj`) %>%
      mutate(
        sig = padj < 0.05,
        # nicer ordering: put significant ones at bottom/top for visibility
        contrast = fct_reorder(contrast, diff)
      )

    # ---- 4) ggplot: highlight significant comparisons ----
    ggplot(tdf, aes(x = diff, y = contrast, xmin = lwr, xmax = upr, color = sig)) +
      geom_segment(aes(x = lwr, xend = upr, yend = contrast), linewidth = 0.7) +
      geom_point(size = 2) +
      geom_vline(xintercept = 0, linetype = 2) +
      scale_color_manual(values = c("FALSE" = "grey30", "TRUE" = "firebrick"),
                         labels = c("Not significant", "Significant (p < 0.05)"),
                         name = NULL) +
      labs(
        x = "Difference in means (first - second)",
        y = NULL,
        title = "Tukey HSD: 95% family-wise CIs for pairwise differences",
        subtitle = "Significant contrasts highlighted"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "bottom",
        panel_grid_minor = element_blank()
      )

    require(ggplot2)
    require(ggsignif)

    ggplot(CDC25, aes(x = cat, y = burn)) +
      geom_boxplot(fill = "green", colour = "darkgreen") +
      scale_x_discrete() + xlab("Cat") +
      ylab("Burn") +
      geom_signif(comparisons = list(c("group_1", "group_2")), 
                  map_signif_level=TRUE)
    #We see that no burn has a higher cat count, however, there could be outliers.
    #Burns from 2000, 2006, and 2018 have higher caterpillar counts as well. 

    ggplot(CDC25, aes(x = cat, y = nfplant, color=burn)) +
      geom_boxplot()

    #This suggests that cats does not differ significantly across burn 
    #years vs. no-burn conditions in this dataset.


    set.seed(123)
    data <- data.frame(
      response = rnorm(100),
      factor1  = factor(rep(c("cat","burn"), each = 50)),
      factor2  = factor(rep(c("cat","burn"), times = 50))
    )

    # add an interaction bump when both are "cat"
    data$response <- with(data, response + ifelse(factor1 == "cat" & factor2 == "cat", 1, 0))

    fit <- aov(response ~ factor1 * factor2, data = data)
    summary(fit)
    #Factor1 is slightly significant and factor2 and factor1:factor2 are significant

    with(data, tapply(response, list(factor1, factor2), mean))  # cell means

    interaction.plot(data$factor1, data$factor2, data$response, ylab="Mean response")
    # or flipped order:
    interaction.plot(data$factor2, data$factor1, data$response)

    #Inspect cell means and plot the interaction
    with(data, tapply(response, list(factor1, factor2), mean))

    interaction.plot(data$factor1, data$factor2, data$response,
                     ylab="Mean response", xlab="factor1", trace.label="factor2")


    #Adjusted Tukey
    install.packages("emmeans")   # if needed
    library(emmeans)

    m <- aov(response ~ factor1*factor2, data=data)

    # Simple effects of factor1 at each level of factor2
    emmeans(m, ~ factor1 | factor2) %>% pairs(adjust="tukey")
    #Because the simple effect of factor1 depends on the level of factor2 (near zero at factor2=burn, 
    #sizeable and opposite-signed at factor2=cat), you likely have a factor1 × factor2 interaction. In 
    #a 2×2 this even looks like a crossover (sign flip).


    # Simple effects of factor2 at each level of factor1
    emmeans(m, ~ factor2 | factor1) %>% pairs(adjust="tukey")
    #The effect of factor2 depends on factor1: it’s ~0 when factor1=burn, but large and negative when 
    #factor1=cat → consistent with a factor1 × factor2 interaction (difference-in-differences ≈ −1.105 
    #− 0.024 ≈ −1.13).

    # Estimated marginal means table
    emm <- emmeans(m, ~ factor1*factor2)
    emm

    # Interaction (DiD) – matches your grid order: (burn,burn), (cat,burn), (burn,cat), (cat,cat)
    int_con <- contrast(emm, list("interaction (DiD)" = c(+1, -1, -1, +1)))
    summary(int_con, infer = c(TRUE, TRUE))  # estimate, SE, df, t, p, CI

    # Simple effects with CIs (what you showed, plus intervals)
    confint(pairs(emmeans(m, ~ factor1 | factor2), adjust = "tukey"))
    confint(pairs(emmeans(m, ~ factor2 | factor1), adjust = "tukey"))

    # Nice interaction plot with 95% CIs
    emmip(m, factor2 ~ factor1, CIs = TRUE)  # or swap axes

    plot(emmip(m, factor2 ~ factor1, CIs = TRUE))

    #Effect size 
    install.packages("effectsize")   # if needed
    library(effectsize)

    eta_squared(m, partial = TRUE)  # partial η² for each term
    #A two-way ANOVA showed a significant factor1 × factor2 interaction.Partial η² = 0.09 (Type III), 
    #95% CI [CI_low, CI_high]. Main effects: factor1, ηp² = 0.04; factor2, ηp² = 0.08. Simple-effects 
    #(Tukey-adjusted) indicated that when factor2 = cat, factor1 = burn was ~0.95 units lower than 
    #factor1 = cat (p = 0.0004), but no difference when factor2 = burn (p = 0.93).

    #Assumptions check 
    par(mfrow=c(1,2))
    plot(m, which=1)   # residuals vs fitted (homogeneity)
    plot(m, which=2)   # QQ plot (normality)-it is normal
    par(mfrow=c(1,1))

    #Figure 
    library(ggplot2)
    emm_df <- as.data.frame(emm)  # from your `emm <- emmeans(m, ~ factor1*factor2)`

    ggplot(emm_df, aes(x = factor1, y = emmean, group = factor2)) +
      geom_point(position = position_dodge(width = 0.2), size = 2) +
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                    width = .05, position = position_dodge(width = 0.2)) +
      geom_line(position = position_dodge(width = 0.2)) +
      labs(y = "Estimated marginal mean of response",
           x = "factor1", color = "factor2",
           title = "Interaction plot with 95% CIs (from emmeans)") +
      theme_minimal(base_size = 12)


    #Effects table
    library(emmeans)
    simple_f1 <- contrast(emmeans(m, ~ factor1 | factor2), method = "revpairwise", adjust = "tukey")
    simple_f2 <- contrast(emmeans(m, ~ factor2 | factor1), method = "revpairwise", adjust = "tukey")

    confint(simple_f1)  # estimates + 95% CIs for factor1 within each factor2
    confint(simple_f2)  # estimates + 95% CIs for factor2 within each factor1

    ##These results indicate that the effect of each factor is contingent on the other, driven by a 
    ###pronounced elevation in the cat–cat condition.


    #standardized effect sizes (Hedges’ g or Cohen’s d)
    library(emmeans)
    library(effectsize)
    library(dplyr)

    # Recreate your contrasts
    f1_contrasts <- emmeans(m, ~ factor1 | factor2) %>% pairs(adjust = "tukey")
    f2_contrasts <- emmeans(m, ~ factor2 | factor1) %>% pairs(adjust = "tukey")

    # Convert the 'summary_emm' objects to data frames
    df_f1 <- as.data.frame(f1_contrasts)
    df_f2 <- as.data.frame(f2_contrasts)

    # Compute Hedges’ g (or Cohen’s d) for each simple effect
    df_f1 <- as.data.frame(f1_contrasts)
    df_f1$g <- effectsize::t_to_d(df_f1$t.ratio, df = df_f1$df, hedges.correction = TRUE)
    df_f1$group <- "factor1 | factor2"

    df_f2 <- as.data.frame(f2_contrasts)
    df_f2$g <- effectsize::t_to_d(df_f2$t.ratio, df = df_f2$df, hedges.correction = TRUE)
    df_f2$group <- "factor2 | factor1"


    # Combine your two contrast data frames first
    contr_all <- rbind(df_f1, df_f2)

    # Helper: compute Hedges' g + 95% CI from a single row's t and df
    calc_g_row <- function(tval, dfval) {
      g_obj <- effectsize::t_to_d(tval, df = dfval, hedges.correction = TRUE, ci = 0.95)
      gd <- as.data.frame(g_obj)  # has columns: d, CI_low, CI_high, ...
      c(g = gd$d, g_low = gd$CI_low, g_high = gd$CI_high)
    }

    # Apply to every row
    g_mat <- t(apply(contr_all[, c("t.ratio","df")], 1, function(x) calc_g_row(x[1], x[2])))
    g_df  <- as.data.frame(g_mat)
    rownames(g_df) <- NULL

    # Bind back to a neat table (choose your columns)
    effect_table_ci <- cbind(
      contr_all[, c("group","contrast","estimate","SE","t.ratio","p.value")],
      g_df
    )

    # Print nicely
    print(effect_table_ci, digits = 3)

    library(ggplot2)

    et <- effect_table_ci
    et$label <- c("factor1 | factor2=burn",
                  "factor1 | factor2=cat",
                  "factor2 | factor1=burn",
                  "factor2 | factor1=cat")

    ggplot(et, aes(x = g, y = label)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_errorbarh(aes(xmin = g_low, xmax = g_high), height = 0.12) +
      geom_point(size = 2) +
      labs(x = "Hedges' g (95% CI)", y = NULL,
           title = "Simple-effects standardized differences") +
      theme_minimal(base_size = 12)

    #Simple-effects contrasts from the emmeans model were converted to Hedges’ g (effect size). When 
    ##factor2 = cat, factor1 = cat exceeded factor1 = burn by 0.95 units (g = 0.74, 95% CI [0.33, 1.16], 
    ### Tukey-adjusted p < .001). When factor1 = cat, factor2 = cat exceeded factor2 = burn by 1.10 units
    ####(g = 0.86, 95% CI [0.44, 1.28], Tukey-adjusted p < .001). In contrast, simple effects were 
    #####negligible when the other factor equaled burn (g ≈ 0.02–0.14, 95% CIs spanning 0). These 
    ######results indicate a practically large effect concentrated in the cat–cat condition, 
    #######consistent with the significant interaction.

## Continuation of frequentest stats with interactions

    glm1 <- glm(cat ~ winter * nfplant, family = poisson)

    summary(glm1)

    ggplot(CDC25, aes(x = cat, y = nfplant*CDC25$`Winter Prior Prep (cm)`)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot()  
    #We see  significance between cats and nfplant


    glm2 <- glm(cat ~ nfplant *cone, family = poisson)

    summary(glm2) 

    ggplot(CDC25, aes(x = cat, y = nfplant*cone)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot()    
    #We see  significance between cats and nfplant 

    glm3 <- glm(cat ~ nfplant *elev, family = poisson)

    summary(glm3) 

    ggplot(CDC25, aes(x = cat, y = nfplant*elev)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot()    
    #We see  significance between cats and nfplant 

    glm4 <- glm(cat ~ nfplant *dbh, family = poisson)

    summary(glm4) 

    ggplot(CDC25, aes(x = cat, y = nfplant*dbh)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot()    
    ##Intercept is  significant meaning the expected value of the dependent variable is 
    #significantly different from zero when all independent variables in the model is 0.
    #Cats are significantly tied to nfplant*dbh

    glm5 <- glm(cat ~ winter + nfplant +dbh+ elev, family = poisson)

    summary(glm5)
    ##Intercept is  slightly significant meaning the expected value of the dependent variable is 
    #significantly different from zero when all independent variables in the model is 0.
    #Cats are fairly significantly tied to winter precep and highly tied to dbh and elev.

    ggplot(CDC25, aes(x = cat, y = mprec + nfplant +dbh+ elev)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot()  
    #We see a significance between cats and nfplant
    --------------------------------------------------------------------------------
    library(DHARMa)
    library(car)
    library(lme4)
    library(lmerTest)
    library(ggplot2)
    library(cowplot)

    # Y = Cone Count
    # X = DBH

    # Simple model with DBH predicting cone count
    mod1 <- lm(cat ~ dbh, data = CDC25)
    summary(mod1)
    #Slightly significant 

    # Check residuals of model
    simulationOutput <- simulateResiduals(mod1, plot = T)
    # Skewed data-deviation significant

    shapiro.test(residuals(mod1))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # Interaction with site
    mod2 <- lm(cat ~ dbh*site, data = CDC25)
    summary(mod2)
    #No significance

    # Check residuals of model
    simulationOutput1 <- simulateResiduals(mod2, plot = T)
    #qu = 0.25, log(sigma) = -3.902867 : outer Newton did not converge fully.

    shapiro.test(residuals(mod2))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # Site as a random effect
    mod3 <- lmer(cat ~ dbh + (1|site), data = CDC25)
    summary(mod3)
    #Slightly significant 

    # Check residuals of model
    simulationOutput2 <- simulateResiduals(mod3, plot = T)
    # Skewed data-deviation significant

    shapiro.test(residuals(mod3))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # monthly precipitation as a random effect with cats and elevation
    mod4 <- lmer(cat ~ elev + (1|mprec), data = CDC25)
    summary(mod4)
    #Not significant 

    # Check residuals of model
    simulationOutput3 <- simulateResiduals(mod4, plot = T)
    # Skewed data-deviation significant 

    shapiro.test(residuals(mod4))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # elevation as a random effect with cats and monthly precipitation
    mod5 <- lmer(cat ~ mprec + (1|elev), data = CDC25)
    summary(mod5)
    #intercept is slightly significant meaning the expected value of the dependent variable is 
    #significantly different from zero when all independent variables in the model is 0.

    # Check residuals of model
    simulationOutput4 <- simulateResiduals(mod5, plot = T)
    # Skewed data-deviation significant 

    shapiro.test(residuals(mod5))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # Interaction with site
    mod6 <- lm(cat ~ elev*nfplant, data = CDC25)
    summary(mod6)
    #Cats is slightly significant in relation to non-focal plants

    # Check residuals of model
    simulationOutput5 <- simulateResiduals(mod6, plot = T)
    #Skewed data-deviation significant

    shapiro.test(residuals(mod6))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # Simple model with elev predicting cone count
    mod7 <- lm(cat ~ nfplant, data = CDC25)
    summary(mod7)
    #confirms significance between cat and non-focal plant species 

    shapiro.test(residuals(mod7))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # Interaction with site
    mod8 <- lm(cat ~ nfplant*site, data = CDC25)
    summary(mod8)
    #No significance

    plot(density(resid(lm(cat ~ nfplant*site, data = CDC25))))

    # Check residuals of model
    simulationOutput1 <- simulateResiduals(mod8, plot = T)
    # Skewed data-deviation significant 

    shapiro.test(residuals(mod8))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # Site as a random effect
    mod9 <- lmer(cat ~ winter + (1|site), data = CDC25)
    summary(mod9)
    ##Intercept is  significant meaning the expected value of the dependent variable is 
    #significantly different from zero when all independent variables in the model is 0.
    #Winter is slightly significant.

    # Check residuals of model
    simulationOutput2 <- simulateResiduals(mod9, plot = T)
    # Skewed data-slightly significant

    shapiro.test(residuals(mod9))
    #Not normally distributed and the predicted estimate is significantly different from 0 

    # monthly precipitation as a random effect with cats and elevation
    mod10 <- lmer(cat ~ elev + (1|mprec), data = CDC25)
    summary(mod10)
    #No significance

    # Check residuals of model
    simulationOutput2 <- simulateResiduals(mod10, plot = T)
    # Skewed data-slightly significant

    shapiro.test(residuals(mod10))
    #Not normally distributed and the predicted estimate is significantly different from 0 


    df <- data.frame(x =cat, y = elev)
    df
    t.test(cat, elev)
    #Slight negative correlation that is marked as significant
    #Meaning, the lower the elevation, there are slightly more cats 

    catelev<- ggplot(CDC25, aes(x = cat, y = elev)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot() 
    #Does show that slight negative correlation 

    catelev

    df1 <- data.frame(x =cat, y = mprec)
    df1
    t.test(cat, mprec)
    #Not significant 

    catmprec<- ggplot(CDC25, aes(x = cat, y = mprec)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot() 
    #Does show that slight negative correlation 

    catmprec

    df2 <- data.frame(x =cat, y = nfplant)
    df2
    t.test(cat, nfplant)
    #Positive correlation that is significant 
    #Meaning, the higher the number of non-focal plants, the higher the number of cats 

    catnfplant <-ggplot(CDC25, aes(x = cat, y = nfplant)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot() 
    #Does show that  positive correlation 

    catnfplant 

    df3 <- data.frame(x =cat, y = winter)
    df3
    t.test(cat, winter)
    #Positive correlation that is significant 
    #Meaning, the higher the number of non-focal plants, the higher the number of cats 

    catwinter<- ggplot(CDC25, aes(x = cat, y = CDC25$`Winter Prior Prep (cm)`)) +
      geom_jitter(width = .2) + geom_smooth(method = "lm", formula = y~x) + theme_cowplot() 
    #Does show that slight positive correlation 

    catwinter

    #From this analysis, we can conclude so far: 
    ##Caterpillar count and monthly precipitation are not correlated or significant
    ##Caterpillar count and cone count are not correlated and significant
    ##Caterpillar count and dbh are not correlated and significant
    ##Caterpillar count and non-focal plant count are correlated and significant
    ##Caterpillar count and winter precep are correlated and significant

## Single-group Bayesian SEM

    #Probability of seeing a caterpillar is distributed binomially, 
    #with probability of having a non-focal plant present when we see that Caterpillar in our
    #80 plots 
      

    Pr(cat|nfplant,80)

    # Check and remove for NAs
    anyNA(CDC25)
    CDC25 <- na.omit(CDC25)  # Remove NAs

    CDC25[is.na(CDC25)] <- mean(CDC25, na.rm = TRUE)

    library(purrr)

    result <- map(CDC25, ~ {
      tryCatch({
        data = (CDC25)  # Replace .x with the actual item from flist if needed
      }, error = function(e) {
        message("Error: ", e$message)
        return(NULL)  # or some default value
      })
    })


    flist<-alist(
      cat~dnorm(mu, sigma),
      mu~dnorm(cat, nfplant),
      sigma~dunif(0,50) #flat prior
    )


    #Multivariate linear model 
    #Looking at the probability of seeing a caterpillar based on the number of non-focal 
    #plants based on site 

    # 0) Load needed pkgs
    library(dplyr)
    library(readr)
    library(janitor)

    # 1) Clean names (you already did this)
    CDC25 <- CDC25 |> janitor::clean_names()

    # 2) Create numeric nfplant and its z-score
    CDC25 <- CDC25 |>
      mutate(
        nfplant_raw = non_focal_tree_plants_present,
        nfplant = if (is.numeric(nfplant_raw)) nfplant_raw
        else readr::parse_number(as.character(nfplant_raw)),
        snfplant = as.numeric(scale(nfplant))
      )

    # 3) Quick checks
    str(CDC25$nfplant)
    summary(CDC25$nfplant)
    head(CDC25$snfplant)


    install.packages(rethinking)
    library(rethinking)
    data<-CDC25
    cat<-CDC25$caterpillar_count

    # z-score the column (mean 0, sd 1)
    CDC25$snfplant <- as.numeric(scale(CDC25$non_focal_tree_plants_present))
    #The new vector snfplant will have mean ≈ 0 and SD ≈ 1.
    #It’s now on the same scale as other standardized predictors—very useful for SEM or 
    #regression, because it: makes coefficients comparable across variables,improves numerical 
    #stability in estimation, and ensures that a 1-unit change always means “one standard deviation.

    #Don't standardize cat. Keep the original (unstandardized) continuous outcome in the dataset 
    #for interpretability, and rely on standardized = TRUE output for standardized comparisons.


    install.packages("lavaan")
    library(lavaan)


    model <- '
     caterpillar_count~ snfplant
    '


    fit <- sem(model, data = CDC25)
    summary(fit, fit.measures = TRUE, standardized = TRUE)
    #Standardized slope (β) ≈ −0.001 → no relationship. P = 0.99 → completely non-significant.
    #R² ≈ 0.00 → snfplant does not explain variability in cat. Model fit perfect but meaningless 
    #(single regression).Conclusion is:The standardized non-focal plant predictor (snfplant)
    #was not associated with cat.


    ##Posterior 
    install.packages(blavaan)
    library(blavaan)

    # Bayesian SEM equivalent to your lavaan model
    model <- '
     caterpillar_count ~ snfplant
    '

    # Fit Bayesian model
    bfit <- bsem(model, data = CDC25, burnin = 1000, sample = 5000)
    #This model fits well

    # Summarize posterior samples
    summary(bfit, standardized = TRUE)
    #The posterior mean for the slope (cat ~ snfplant) is −0.001,with a 95% credible interval from 
    #−0.214 to +0.211.Given the wide, symmetric interval centered on zero, the data provide no 
    #evidence for a meaningful relationship between snfplant and cat.The posterior predictive p-value 
    #(PPP = 0.505) indicates that the model fits the data well — the lack of effect is genuine, 
    #not due to model misfit.

    ###Visualize posterior
    # --- Get posterior draws from blavaan ---
    # Returns a list of matrices: one matrix per chain (iterations x parameters)
    post_list <- blavaan::blavInspect(bfit, "mcmc")

    # Bind chains into one matrix of draws
    post_mat  <- do.call(rbind, post_list)

    # See what the parameter column is called
    grep("cat.*~.*snfplant", colnames(post_mat), value = TRUE)
    par_name <- grep("caterpillar_count*~.*snfplant", colnames(post_mat), value = TRUE)[1]

    # --- Plot with bayesplot ---
    library(bayesplot)
    bayesplot::mcmc_areas(
      post_mat,
      pars = "caterpillar_count~snfplant",
      prob = 0.95
    ) +
      ggtitle("Posterior for cat ~ snfplant") +
      xlab("Standardized coefficient")


    ##SEM
    semPaths(fit, 
             what = "std",               # Standardized estimates
             residuals = TRUE,           # Show residuals
             edge.label.cex = 0.8,       # Font size for edge labels
             sizeMan = 8,                # Size of manifest variables
             sizeLat = 10,               # Size of latent variables
             layout = "tree",            # Layout of the graph
             intercepts = FALSE)         # Exclude intercepts from the plot


    --------------------------------------------------------------------------------

    #Probability of seeing a caterpillar is distributed binomially, 
    #with probability of having a burn present when we see that Caterpillar in our
    #80 plots 

      #Set working directory
      setwd("/Users/chloecollier/Desktop")

    #Read in the file we will be working with
    library(readr)
    CDC25<- read_csv("Combined-Data-Cat.csv")
      
    # Install and load fastDummies package
    install.packages("fastDummies")
    library(fastDummies)

    #We use dummy variables so that categorical groups (like seasons or treatment arms) can be 
    ##represented numerically in regression/SEM (0 or 1), allowing the model to estimate each group’s effect 
    ###relative to a baseline in a statistically valid way.

    # 1) Normalize weird spaces in all column names
    names(CDC25) <- gsub("\u00A0", " ", names(CDC25))    # replace NBSP with regular space
    names(CDC25) <- gsub("[[:space:]]+", " ", names(CDC25))  # collapse any runs of spaces/tabs
    names(CDC25) <- trimws(names(CDC25))                 # trim ends

    # 2) Verify the exact column name now
    grep("Year", names(CDC25), value = TRUE)
    # Expect to see: "Year Burned"

    # 3) Create dummies for that column
    library(fastDummies)
    # Make sure Year Burned is a factor so levels are stable
    CDC25$`Year Burned` <- factor(CDC25$`Year Burned`)
    cat <- factor(CDC25$`Caterpillar Count`)

    CDC25_dummies <- dummy_cols(
      CDC25,
      select_columns = "Year Burned",
      remove_first_dummy = TRUE,      # <- avoids multicollinearity
      remove_selected_columns = TRUE  # <- drops the original "Year Burned"
    )

    colnames(CDC25_dummies)           # look for "Year_Burned_2024", "Year_Burned_2025", etc.

    CDC25_dummies$z_dbh  <- scale(CDC25_dummies$`dbh(m)`)
    CDC25_dummies$z_ysb  <- scale(CDC25_dummies$`Monthly Prep (cm)`)  # or your chosen YSB variable
    CDC25_dummies$z_rich <- scale(CDC25_dummies$`Caterpillar Count`)  # example if you need a z version


    model1 <- '
      cat ~ CDC25$`Year Burned` 

    '

    # Fit the model
    fit <- sem(model, data = CDC25_dummies)

    # Summarize results
    summary(fit, fit.measures = TRUE, standardized = TRUE)

    install.packages("semPlot")
    library(semPlot)

    semPaths(fit, 
             what = "std",               # Standardized estimates
             residuals = TRUE,           # Show residuals
             edge.label.cex = 0.8,       # Font size for edge labels
             sizeMan = 8,                # Size of manifest variables
             sizeLat = 10,               # Size of latent variables
             layout = "tree",            # Layout of the graph
             intercepts = FALSE)         # Exclude intercepts from the plot

    library(dplyr)
    library(tidyr)
    library(ggplot2)

    post <- as.matrix(fit)              # posterior draws from brms
    colnames(post)[1:20]                # inspect names

    # Pick & relabel coefficients (adjust names to your model)
    burn_names <- c("b_burnedBurned","b_burned1","b_burnedTRUE")
    burn_present <- intersect(burn_names, colnames(post))

    df_long <- as.data.frame(post) |>
      dplyr::transmute(
        Intercept          = .data[["b_Intercept"]],
        `Precip Effect`    = .data[["b_mprec"]],
        `Non-focal Effect` = .data[["b_nfplant"]],
        `DBH Effect`       = .data[["b_dbh"]],
        `Burned Effect`    = if (length(burn_present)) .data[[burn_present]] else NA_real_
      ) |>
      tidyr::pivot_longer(everything(), names_to = "Parameter", values_to = "Posterior") |>
      dplyr::filter(!is.na(Posterior))

    # Thin points if needed (version-robust)
    set.seed(1)
    df_long <- df_long |>
      dplyr::group_by(Parameter) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::filter(row_id %in% sample(row_id, min(1500, dplyr::n()), replace = FALSE)) |>
      dplyr::ungroup() |>
      dplyr::select(-row_id)

    ggplot(df_long, aes(x = Parameter, y = Posterior)) +
      geom_boxplot(fill = "skyblue", alpha = 0.5, width = 0.55, outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.18, height = 0),
                 alpha = 0.45, size = 1.1, color = "blue") +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = "Posterior Distribution for Caterpillar Model",
           x = "Parameters", y = "Posterior Samples") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 28, hjust = 1))
    ------------------------------------------------------------------------------

      #Probability of seeing a caterpillar is distributed binomially, 
    #with probability of having non-focal species present when we see that Caterpillar in our
    #80 plots 


    # Load required packages
    library(lavaan)
    library(semPlot)

    #str(CDC25$`Caterpillar Count`)--> continuous do not make as a dummy variable 
    #str(CDC25$`Non-focal tree plants present`)--> continuous do not make as a dummy variable 

    set.seed(1)
    CDC25 <- data.frame(
      x1 = rnorm(100),
      x2 = rnorm(100),
      # treat these as continuous variables (NOT dummies)
      cat      = rnorm(100),   # e.g., a continuous score
      nfplant  = rnorm(100),   # e.g., a continuous index
      y1 = rnorm(100),
      y2 = rnorm(100)
    )

    model_cont <- '
      y1 ~ x1 + x2 + cat + nfplant
      y2 ~ x1 + x2 + cat + nfplant
    '

    fit_cont <- lavaan::sem(model_cont, data = CDC25, estimator = "MLR")
    summary(fit_cont, fit.measures = TRUE, standardized = TRUE)
    #df = 0 (saturated) → χ²=0, RMSEA=0, CFI/TLI=1 are trivial (don’t inform fit).
    #So the linear model explains little variance here.


    summary(fit_cont, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
    # or:
    inspect(fit_cont, "r2")


    # install.packages("car")
    car::vif(lm(y1 ~ x1 + x2 + cat + nfplant, data = CDC25))
    car::vif(lm(y2 ~ x1 + x2 + cat + nfplant, data = CDC25))

    # quick polynomial terms
    CDC25$cat2     <- CDC25$cat^2
    CDC25$nfplant2 <- CDC25$nfplant^2

    model_poly <- '
      y1 ~ x1 + x2 + cat + cat2 + nfplant + nfplant2
      y2 ~ x1 + x2 + cat + cat2 + nfplant + nfplant2
    '
    fit_poly <- sem(model_poly, data = CDC25, estimator = "MLR")
    summary(fit_poly, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)


    CDC25$x1_cat <- CDC25$x1 * CDC25$cat
    CDC25$x2_cat <- CDC25$x2 * CDC25$cat

    model_int <- '
      y1 ~ x1 + x2 + cat + nfplant + x1_cat + x2_cat
      y2 ~ x1 + x2 + cat + nfplant + x1_cat + x2_cat
    '
    fit_int <- sem(model_int, data = CDC25, estimator = "MLR")
    summary(fit_int, fit.measures = TRUE, standardized = TRUE)


    numvars <- c("x1","x2","cat","nfplant")
    CDC25[ paste0("z_", numvars) ] <- lapply(CDC25[numvars], scale)

    model_z <- '
      y1 ~ z_x1 + z_x2 + z_cat + z_nfplant
      y2 ~ z_x1 + z_x2 + z_cat + z_nfplant
    '
    fit_z <- sem(model_z, data = CDC25, estimator = "MLR")
    summary(fit_z, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

    #Check to make sure multicollinerity is ok (which it is)
    #BUT all three models are saturated and uninformative. 

    ###NEXT STEPS-check robustness and interpretability 


    #Center, rebuils terms, refit 
    library(lavaan)
    dat <- within(CDC25, {
      cx1      <- scale(x1, center = TRUE, scale = FALSE)
      cx2      <- scale(x2, center = TRUE, scale = FALSE)
      ccat     <- scale(cat, center = TRUE, scale = FALSE)
      cnfplant <- scale(nfplant, center = TRUE, scale = FALSE)
      
      ccat2     <- as.numeric(ccat)^2
      cnfplant2 <- as.numeric(cnfplant)^2
      cx1_ccat  <- as.numeric(cx1) * as.numeric(ccat)
      cx2_ccat  <- as.numeric(cx2) * as.numeric(ccat)
    })

    model_poly_c <- '
      y1 ~ cx1 + cx2 + ccat + ccat2 + cnfplant + cnfplant2
      y2 ~ cx1 + cx2 + ccat + ccat2 + cnfplant + cnfplant2
    '
    fit_poly_c <- sem(model_poly_c, data = dat, estimator = "MLR")

    model_int_c <- '
      y1 ~ cx1 + cx2 + ccat + cnfplant + cx1_ccat + cx2_ccat
      y2 ~ cx1 + cx2 + ccat + cnfplant + cx1_ccat + cx2_ccat
    '
    fit_int_c <- sem(model_int_c, data = dat, estimator = "MLR")

    summary(fit_poly_c, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
    summary(fit_int_c,  fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)

    # Compare by AIC/BIC
    AIC(fit_poly_c, fit_int_c); BIC(fit_poly_c, fit_int_c)

    #After centering, we detected a small but significant x2 × cat interaction predicting y₂ (β = 0.31,
    #p = 0.007; R² = 0.10) and a concave quadratic effect of nfplant on y₁ (β for nfplant² = −0.18,
    #p = 0.035; R² = 0.09); all other terms were non-significant, and information criteria favored the 
    #interaction specification.”

    #Probe the y2 interaction (simple slopes)
    # Get coefficients
    b <- parameterEstimates(fit_int_c)
    bx2     <- b$est[b$lhs=="y2" & b$rhs=="cx2"      & b$op=="~"]
    bx2_cc  <- b$est[b$lhs=="y2" & b$rhs=="cx2_ccat" & b$op=="~"]

    # SD of centered cat is original sd(cat)
    sd_cat <- sd(CDC25$cat, na.rm=TRUE)

    # Simple slopes
    slope_mean <- bx2 + bx2_cc * 0
    slope_low  <- bx2 + bx2_cc * (-sd_cat)
    slope_high <- bx2 + bx2_cc * ( sd_cat)

    c(slope_low=slope_low, slope_mean=slope_mean, slope_high=slope_high)

    #Slopes + SEs + 95% CIs programmatically
    # --- pull estimates and the vcov in the SAME naming system as V ---
    pe <- parameterEstimates(fit_int_c)              # for the point estimates
    V  <- lavInspect(fit_int_c, "vcov")              # var-cov of free params (robust with MLR)

    # grab the two coefficients we need: y2 ~ cx2  and  y2 ~ cx2_ccat
    b_x2    <- pe$est[pe$lhs=="y2" & pe$op=="~" & pe$rhs=="cx2"]
    b_int   <- pe$est[pe$lhs=="y2" & pe$op=="~" & pe$rhs=="cx2_ccat"]

    # find their positions inside V using rownames(V)
    ix_bx2  <- which(rownames(V) == "y2~cx2")
    ix_bint <- which(rownames(V) == "y2~cx2_ccat")

    # sanity checks (should be length 1 each)
    stopifnot(length(ix_bx2)==1, length(ix_bint)==1)

    var_bx2  <- V[ix_bx2, ix_bx2]
    var_bint <- V[ix_bint, ix_bint]
    cov_b    <- V[ix_bx2, ix_bint]

    # simple-slope function: slope(c) = b_x2 + c * b_int
    # Var[slope(c)] = Var(b_x2) + c^2 Var(b_int) + 2 c Cov(b_x2, b_int)
    ss_row <- function(c){
      slope <- b_x2 + c*b_int
      se    <- sqrt(var_bx2 + (c^2)*var_bint + 2*c*cov_b)
      z     <- slope / se
      p     <- 2*pnorm(-abs(z))
      ci    <- slope + c(-1,1)*1.96*se
      c(slope=slope, se=se, z=z, p=p, ci_lo=ci[1], ci_hi=ci[2])
    }

    sd_cat <- sd(CDC25$cat, na.rm=TRUE)
    out <- rbind(
      low  = ss_row(-sd_cat),
      mean = ss_row(0),
      high = ss_row(sd_cat)
    )
    round(out, 4)

    #Johnson–Neyman (where is it significant?)
    # --- get needed pieces from fit_int_c ---
    pe <- parameterEstimates(fit_int_c)
    V  <- lavInspect(fit_int_c, "vcov")

    b_x2    <- pe$est[pe$lhs=="y2" & pe$op=="~" & pe$rhs=="cx2"]
    b_int   <- pe$est[pe$lhs=="y2" & pe$op=="~" & pe$rhs=="cx2_ccat"]
    ix_bx2  <- which(rownames(V)=="y2~cx2")
    ix_bint <- which(rownames(V)=="y2~cx2_ccat")
    var_bx2 <- V[ix_bx2, ix_bx2]
    var_bint<- V[ix_bint, ix_bint]
    cov_b   <- V[ix_bx2, ix_bint]

    # z(c) = [b_x2 + c*b_int] / sqrt(var_bx2 + c^2*var_bint + 2*c*cov_b)
    zfun <- function(c) {
      num <- b_x2 + c*b_int
      den <- sqrt(var_bx2 + c^2*var_bint + 2*c*cov_b)
      num / den
    }

    # search over observed range of centered cat
    cc <- scale(CDC25$cat, center=TRUE, scale=FALSE)
    rng <- range(as.numeric(cc), na.rm=TRUE)

    # find roots for z = ±1.96 (two-sided 95%)
    f1 <- function(c) zfun(c) - 1.96
    f2 <- function(c) zfun(c) + 1.96

    cand <- c()
    for (f in list(f1, f2)) {
      # sample grid to find sign changes
      grid <- seq(rng[1], rng[2], length.out = 200)
      vals <- sapply(grid, f)
      sgn  <- sign(vals)
      flips <- which(diff(sgn)!=0 & is.finite(vals[-1]) & is.finite(vals[-length(vals)]))
      if (length(flips)) {
        for (i in flips) {
          r <- try(uniroot(f, c(grid[i], grid[i+1]))$root, silent=TRUE)
          if (!inherits(r, "try-error")) cand <- c(cand, r)
        }
      }
    }
    sort(unique(round(cand, 4)))

    library(ggplot2)

    slope_ci <- function(c){
      slope <- b_x2 + c*b_int
      se    <- sqrt(var_bx2 + c^2*var_bint + 2*c*cov_b)
      c(slope=slope, lo=slope-1.96*se, hi=slope+1.96*se)
    }

    grid <- data.frame(ccat = seq(rng[1], rng[2], length.out=200))
    tmp  <- t(apply(as.matrix(grid$ccat), 1, slope_ci))
    grid$slope <- tmp[, "slope"]
    grid$lo    <- tmp[, "lo"]
    grid$hi    <- tmp[, "hi"]

    ggplot(grid, aes(ccat, slope)) +
      geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.2) +
      geom_hline(yintercept=0, linetype="dashed") +
      geom_line(size=1) +
      labs(x="cat (centered)", y="Conditional slope of x2 → y2") +
      theme_minimal()

    ##Johnson–Neyman analysis showed that the conditional effect of x2 on y2 was significantly 
    #negative when cat < −0.54 (centered units), non-significant between −0.54 and 1.49, and 
    #significantly positive when cat > 1.49. Thus, higher cat values reverse and strengthen the 
    #x2→y2 association (hold variables for comparison)


    #Probability of seeing a caterpillar is distributed binomially, 
    #with probability of having a winter precipitation when we see that Caterpillar in our
    #80 plots 

## Multi-group Bayesian SEM (for 2023–2025)

    # --- Packages ---
    if (!requireNamespace("blavaan", quietly=TRUE)) install.packages("blavaan")
    if (!requireNamespace("bayesplot", quietly=TRUE)) install.packages("bayesplot")
    if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
    if (!requireNamespace("semPlot", quietly=TRUE)) install.packages("semPlot")

    library(blavaan)
    library(bayesplot)
    library(dplyr)
    library(semPlot)
      
      
    #Probability of seeing a caterpillar is distributed binomially, 
    #with probability of having non-focal species present with the additive of 
    # monthly precipitation when we see that Caterpillar in our
    #80 plots 

    # Load required packages
    library(lavaan)
    library(semPlot)

    # Create your dataset (adjust as necessary)
    CDC25<- data.frame(
      nfplant = rnorm(100),
      cat = rnorm(100),
      mprec=rnorm(100)
      # Exclude categorical variables
    )

    # Define the model without dummy variables
    model3 <- '
      cat ~ nfplant+mprec

    '

    # Fit the SEM model using the original data
    fit <- sem(model3, data = CDC25)

    # Plot the SEM model
    semPaths(fit, 
             what = "std", 
             residuals = TRUE, 
             edge.label.cex = 0.8, 
             sizeMan = 8, 
             sizeLat = 10, 
             layout = "tree", 
             intercepts = FALSE)

    sum(is.na(CDC25))

    # Summarize posterior samples
    summary(fit, standardized = TRUE)
    #There’s at best a weak tendency for higher nfplant to be associated with higher cat; 
    #mprec doesn’t add much once nfplant is in.

    m <- lm(cat ~ nfplant + mprec, data = CDC25)
    summary(m)          # R², p's
    car::Anova(m, type=2)  # partial sums of squares
    #Not significant


    # --- 1) Make sure season is numeric/factor with 3 years (2023–2025)
    table(df_mod$season)

    # --- 2) Multilevel Bayesian model ---
    # Level-1 (within years): plot-level predictors
    # Level-2 (between years): z_winter varies between years only
    # Priors: Normal(0,1) on standardized slopes; LKJ(2) on correlations

    model_bayes_ml <- '
      # =====================
      level: 1
      # =====================
        z_richness ~ prior("normal(0,1)") * a2 * z_ysb
        z_cones    ~ prior("normal(0,1)") * b1 * z_dbh
                    + prior("normal(0,1)") * b2 * z_ysb
        cat        ~ prior("normal(0,1)") * c1 * z_richness
                    + prior("normal(0,1)") * c2 * z_cones
                    + prior("normal(0,1)") * c3 * z_ysb
                    + prior("normal(0,1)") * c5 * z_dbh



      # =====================
      level: 2
      # =====================
        # between-year effects of winter
        z_richness ~ prior("normal(0,1)") * a1 * z_winter
        z_cones    ~ prior("normal(0,1)") * bw * z_winter
        cat        ~ prior("normal(0,1)") * c4 * z_winter

    '

    # --- 3) Fit Bayesian multilevel SEM ---
    bfit_ml <- blavaan::bsem(
      model_bayes_ml,
      data      = df_mod,
      cluster   = "season",     # year-level clustering
      n.chains  = 3,
      burnin    = 2000,
      sample    = 5000,
      adapt     = 1000,
      seed      = 42
    )

    # --- 4) Summaries ---
    summary(bfit_ml, fit.measures = TRUE, standardized = TRUE, priors = TRUE)

    ##Not great results, lets try again. Need tighter priors.-----------------------

    dp_ok <- dpriors(
      # paths/loadings/intercepts (z-scales)
      beta   = "normal(0,1)",
      lambda = "normal(0,1)",
      alpha  = "normal(0,1)",
      nu     = "normal(0,1)",
      # SD priors: keep gamma, but make them more informative than gamma(1,.5)
      # Center roughly near 1 SD with shorter tails
      theta  = "gamma(2,0.5)[sd]",  # observed residual SDs  (mean = 1)
      psi    = "gamma(2,0.5)[sd]"   # latent SDs            (mean = 1)
    )

    bfit_ml2 <- blavaan::bsem(
      model_bayes_ml,
      data    = df_mod,
      cluster = "season",
      n.chains = 3, burnin = 3000, sample = 8000, adapt = 1500,
      dp = dp_ok, seed = 42
    )

    # --- 4) Summaries ---
    summary(bfit_ml2, fit.measures = TRUE, standardized = TRUE, priors = TRUE)

    ##Diagnostics look fine: R̂ = 1.00 across the board; mixing is OK. The issue is
    #identifiability at Level-2, not convergence.

    ##Try again with season as a fixed effect---------------------------------------

    # Make 2023 the reference (baseline) level
    df_mod$season <- relevel(df_mod$season, ref = "2023")

    # Build design matrix with intercept
    X <- model.matrix(~ season, data = df_mod)
    colnames(X)
    # -> "(Intercept)", "season2024", "season2025"

    # Add the two dummy columns to df_mod
    df_mod$season2024 <- X[, "season2024"]
    df_mod$season2025 <- X[, "season2025"]

    # In your lavaan syntax, add these as covariates on each DV:
    model_bayes_fe <- "
      z_richness ~ z_ysb + season2024 + season2025
      z_cones    ~ z_dbh + z_ysb + season2024 + season2025
      cat        ~ z_richness + z_cones + z_ysb + z_dbh + season2024 + season2025
    "

    --------------------------------------------------------------------------------
    #Probability of seeing a caterpillar is distributed binomially, 
    #with probability of having non-focal species present with the additive of 
    # monthly precipitation when we see that Caterpillar in our
    #80 plots 

    library(blavaan)


    dp_fe <- dpriors(
      beta   = "normal(0,1)",
      lambda = "normal(0,1)",
      alpha  = "normal(0,1)",
      nu     = "normal(0,1)",
      theta  = "gamma(2,0.5)[sd]"   # (mean ~1 on SD scale; a bit tighter than default)
    )

    bfit_fe <- blavaan::bsem(
      model_bayes_fe,
      data     = df_mod,
      n.chains = 3,
      burnin   = 3000,
      sample   = 8000,
      adapt    = 1500,
      dp       = dp_fe,   # ← ensure this is here
      seed     = 42
    )

    summary(bfit_fe, fit.measures = TRUE, standardized = TRUE, priors = TRUE)
    #The model fit well (PPP = 0.59; all R̂ = 1.00). z_ysb positively predicted cones (β = 0.49, 95%
    ##CrI 0.19–0.79). Cones were lower in 2024 vs 2023 (β = −0.71, 95% CrI −1.36 to −0.05). Other 
    ###effects’ CrIs overlapped zero.





    library(dplyr)
    library(janitor)

    CDC25 <- CDC25 %>% janitor::clean_names()

    names(CDC25)
    # now you'll have:
    # "plot_num", "plot_diam", "site", "lat", "long", "elev_ft",
    # "date", "center_plant_gen", "dbh_m", "cone_count", "caterpillar_count",
    # "day_prep_cm", "monthly_prep_cm", "winter_precip_cm",
    # "year_burned", "non_focal_tree_plants_present"

    if (is.logical(CDC25$non_focal_tree_plants_present)) {
      CDC25$non_focal_tree_plants_present <- as.integer(CDC25$non_focal_tree_plants_present)
    } else if (is.factor(CDC25$non_focal_tree_plants_present)) {
      CDC25$non_focal_tree_plants_present <- as.integer(CDC25$non_focal_tree_plants_present) - 1
    }


    library(lavaan)

    model4 <- '
      caterpillar_count ~ non_focal_tree_plants_present + winter_precip_cm
    '

    fit <- sem(model4, data = CDC25, estimator = "MLR", missing = "fiml")
    summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)


    # Plot the SEM model
    semPaths(fit, 
             what = "std", 
             residuals = TRUE, 
             edge.label.cex = 0.8, 
             sizeMan = 8, 
             sizeLat = 10, 
             layout = "tree", 
             intercepts = FALSE)


    library(MASS)

    # Use your cleaned names (adjust if yours differ)
    m_pois <- glm(caterpillar_count ~ non_focal_tree_plants_present + winter_precip_cm,
                  data = CDC25, family = poisson())

    m_nb   <- MASS::glm.nb(caterpillar_count ~ non_focal_tree_plants_present + winter_precip_cm,
                           data = CDC25)

    AIC(m_pois, m_nb)    # m_nb is lower
    summary(m_nb)        # check coeff signs, p-values

    #Negative-binomial regression indicated that caterpillar counts increased with non-focal
    #tree plants (IRR = 1.004 per unit, 95% CI 1.002–1.006, p < .001) and decreased with winter 
    #precipitation (IRR = 0.94 per cm, 95% CI 0.88–1.00, p = .046). Overdispersion (θ ≈ 0.68) 
    #justified the NB link; model fit was adequate (Deviance 78.7/78; AIC 282).




    library(DHARMa)
    sim <- simulateResiduals(m_nb, n = 1000)
    plot(sim)
    testDispersion(sim); testZeroInflation(sim)

    newdat <- with(CDC25, expand.grid(
      non_focal_tree_plants_present = quantile(non_focal_tree_plants_present, c(.1,.5,.9), na.rm=TRUE),
      winter_precip_cm = seq(min(winter_precip_cm, na.rm=TRUE), max(winter_precip_cm, na.rm=TRUE), length.out=100)
    ))
    newdat$pred <- predict(m_nb, newdata=newdat, type="response")

    #DHARMa dispersion test: 1.05, p = 0.522 → no over/under-dispersion problem.
    #Zero-inflation test: ratio ≈ 1.007, p = 1 → no excess zeros.
    #Residual plots ran without complaints → diagnostics look fine.


    library(ggplot2)

    # label the non_focal levels
    nf_q <- quantile(CDC25$non_focal_tree_plants_present, c(.1,.5,.9), na.rm=TRUE)
    newdat$nf_level <- factor(newdat$non_focal_tree_plants_present,
                              levels = as.numeric(nf_q),
                              labels = c("Low (10th pct)", "Median (50th)", "High (90th)"))

    # get link-scale preds + SE for CIs, then transform
    pr <- predict(m_nb, newdata = newdat, type = "link", se.fit = TRUE)
    newdat$fit_link <- pr$fit
    newdat$se_link  <- pr$se.fit
    newdat$fit_resp <- exp(newdat$fit_link)
    newdat$lo_resp  <- exp(newdat$fit_link - 1.96*newdat$se_link)
    newdat$hi_resp  <- exp(newdat$fit_link + 1.96*newdat$se_link)

    ggplot(newdat, aes(winter_precip_cm, fit_resp, linetype = nf_level)) +
      geom_ribbon(aes(ymin = lo_resp, ymax = hi_resp), alpha = 0.15) +
      geom_line(linewidth = 1) +
      labs(x = "Winter precipitation (cm)",
           y = "Predicted caterpillar count",
           linetype = "Non-focal plants") +
      theme_minimal()



    library(ggplot2)

    # ensure nf_level is ordered properly
    newdat$nf_level <- factor(
      newdat$nf_level,
      levels = c("Low (10th pct)", "Median (50th)", "High (90th)")
    )

    ggplot(newdat, aes(x = winter_precip_cm, y = fit_resp, color = nf_level, fill = nf_level)) +
      # confidence ribbon with matching color but transparent
      geom_ribbon(aes(ymin = lo_resp, ymax = hi_resp, fill = nf_level),
                  alpha = 0.15, color = NA) +
      geom_line(linewidth = 1) +
      
      # clear labeling and legend title
      labs(
        x = "Winter precipitation (cm)",
        y = "Predicted caterpillar count",
        color = "Non-focal plants",
        fill = "Non-focal plants",
        title = "Effect of Winter Precipitation on Predicted Caterpillar Count",
        subtitle = "Lines show fitted values; shaded areas show 95% confidence intervals"
      ) +
      
      # nice clean minimal theme
      theme_minimal(base_size = 15) +
      theme(
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      
      # distinct, accessible color palette
      scale_color_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2")








    co <- coef(m_nb)
    irr_1cm  <- exp(co["winter_precip_cm"])                 # per +1 cm
    irr_5cm  <- exp(5*co["winter_precip_cm"])               # per +5 cm
    irr_nf1  <- exp(co["non_focal_tree_plants_present"])    # per +1 unit (plants)

    # Optional: across observed spans
    rng_precip <- range(CDC25$winter_precip_cm, na.rm=TRUE)
    delta_precip <- diff(rng_precip)
    irr_span <- exp(delta_precip * co["winter_precip_cm"])

    list(
      IRR_per_1cm   = irr_1cm,
      IRR_per_5cm   = irr_5cm,
      IRR_per_1plant= irr_nf1,
      IRR_over_precip_span = irr_span,
      precip_range_cm = rng_precip
    )

    #Winter precipitation: IRR ≈ exp(-0.0648) ≈ 0.94 per cm (≈ 6% fewer caterpillars for each additional cm).
    #Non-focal plants: IRR ≈ exp(0.0040) ≈ 1.004 per unit (≈ 0.4% more per additional plant).
    #Diagnostics (DHARMa) indicate no dispersion or zero-inflation issues; NB is appropriate.

    --------------------------------------------------------------------------------

      #Probability of seeing a caterpillar is distributed binomially, 
      #with probability of having non-focal species present with the additive of 
      # monthly precipitation, winter precip, DBH when we see that Caterpillar in our
      #80 plots 
      
      # Load required packages
      library(lavaan)
    library(semPlot)

    names(CDC25)[names(CDC25) == "dbh_m"] <- "dbh"
    names(CDC25)[names(CDC25) == "Winter Prior Prep (cm)"] <- "winter"

    # Create your dataset (adjust as necessary)
    CDC25<- data.frame(
      nfplant = rnorm(100),
      cat = rnorm(100),
      mprec=rnorm(100),
      dbh=rnorm(100),
      winter=rnorm(100)
      # Exclude categorical variables
    )

    # Define the model without dummy variables
    model3 <- '
      cat ~ nfplant+mprec+dbh+winter

    '

    # Fit the SEM model using the original data
    fit <- sem(model3, data = CDC25)

    # Plot the SEM model
    semPaths(fit, 
             what = "std", 
             residuals = TRUE, 
             edge.label.cex = 0.8, 
             sizeMan = 8, 
             sizeLat = 10, 
             layout = "tree", 
             intercepts = TRUE)



    # Summarize posterior samples
    summary(fit, standardized = TRUE)
    #There’s at best a weak tendency for higher nfplant to be associated with higher cat; 
    #mprec doesn’t add much once nfplant is in.

    m <- lm(cat ~ nfplant+mprec+dbh+winter, data = CDC25)
    summary(m)          # R², p's
    car::Anova(m, type=2)  # partial sums of squares
    #Only winter is significant 


    # install.packages("broom") # if needed
    library(ggplot2)
    library(broom)
    library(dplyr)

    # Tidy model with confidence intervals
    tcoefs <- broom::tidy(m, conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(term = recode(term,
                           nfplant = "Non-focal plants",
                           mprec   = "Monthly precip",
                           dbh     = "DBH",
                           winter  = "Winter precip"))

    ggplot(tcoefs, aes(x = estimate, y = reorder(term, estimate))) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
      labs(
        title = "Linear model coefficients (±95% CI)",
        x = "Estimate",
        y = NULL
      ) +
      theme_minimal(base_size = 14)




    # Build a prediction grid over winter; hold others at medians
    pred_grid <- data.frame(
      winter  = seq(min(CDC25$winter, na.rm = TRUE),
                    max(CDC25$winter, na.rm = TRUE),
                    length.out = 200),
      nfplant = median(CDC25$nfplant, na.rm = TRUE),
      mprec   = median(CDC25$mprec,   na.rm = TRUE),
      dbh     = median(CDC25$dbh,     na.rm = TRUE)
    )

    # Predict (identity link for gaussian lm)
    pr <- predict(m, newdata = pred_grid, se.fit = TRUE)
    pred_grid$fit <- pr$fit
    pred_grid$lo  <- pr$fit - 1.96 * pr$se.fit
    pred_grid$hi  <- pr$fit + 1.96 * pr$se.fit

    ggplot(pred_grid, aes(x = winter, y = fit)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
      geom_line(linewidth = 1) +
      labs(
        title = "Marginal effect of winter precipitation",
        subtitle = "Other covariates held at their medians",
        x = "Winter precipitation (standardized units, in your sim)",
        y = "Predicted cat"
      ) +
      theme_minimal(base_size = 14)





    ##Graphing posteriors

    library(dplyr)
    library(tidyr)
    library(ggplot2)

    # post can be a matrix or data.frame with columns for each parameter
    # Example: post <- as.matrix(fit)  # from brms/rstan

    # 1) Pick the columns that correspond to your coefficients.
    #    The regex handles names like b_mprec, b_dbh.m., etc.
    keep <- grep("^b?[_\\.]*?(Intercept|precip|mprec|burn|nfplant|non.?focal|dbh)",
                 colnames(CDC25), ignore.case = TRUE, value = TRUE)

    df_long <- as.data.frame(CDC25) |>
      dplyr::select(all_of(keep)) |>
      dplyr::rename(
        Intercept          = dplyr::matches("Intercept", ignore.case = TRUE),
        `Precip Effect`    = dplyr::matches("precip|mprec", ignore.case = TRUE),
        `Burned Effect`    = dplyr::matches("burn",   ignore.case = TRUE),
        `Non-focal Effect` = dplyr::matches("nfplant|non.?focal", ignore.case = TRUE),
        `DBH Effect`       = dplyr::matches("\\bdbh\\b|dbh", ignore.case = TRUE)
      ) |>
      tidyr::pivot_longer(everything(), names_to = "Parameter", values_to = "Posterior")


    df_long <- df_long |>
      dplyr::group_by(Parameter) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::filter(row_id %in% sample(row_id, min(1500, dplyr::n()), replace = FALSE)) |>
      dplyr::ungroup() |>
      dplyr::select(-row_id)   # <-- explicitly dplyr




    # (Optional) Subsample points so the scatter isn't too dense
    set.seed(1)
    df_long <- df_long |>
      group_by(Parameter) |>
      mutate(row_id = row_number()) |>
      filter(row_id %in% sample(row_id, min(1500, n()), replace = FALSE)) |>
      ungroup() |>
      select(-row_id)

    # 2) Plot: boxplot + jitter, with labels and theme
    ggplot(df_long, aes(x = Parameter, y = Posterior)) +
      geom_boxplot(fill = "skyblue", alpha = 0.5, width = 0.55, outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.18, height = 0),
                 alpha = 0.45, size = 1.1, color = "blue") +
      labs(
        title = "Posterior Distribution for Caterpillar Model",
        x = "Parameters",
        y = "Posterior Samples"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 28, hjust = 1))


    # install.packages("brms")  # if needed
    library(brms)

    # Example outcome (make one if you don’t have it yet)
    # Say "cat" is your response (standardize if appropriate)
    if (!"cat" %in% names(CDC25)) CDC25$cat <- rnorm(nrow(CDC25), 0.2*scale(CDC25$nfplant) + rnorm(nrow(CDC25),0,1))

    set.seed(123)

    CDC25$burned <- factor(
      rbinom(nrow(CDC25), 1, 0.45),
      labels = c("No burn", "Burned")
    )

    table(CDC25$burned)

    # Model with continuous and categorical predictors
    fit <- brm(
      cat ~ mprec + nfplant + dbh + burned,
      data = CDC25,
      family = gaussian(),
      cores = 4, chains = 4, iter = 2000, seed = 123
    )


    library(bayesplot)
    library(dplyr)
    library(tidyr)
    library(ggplot2)

    post <- as.matrix(fit)

    # Grab and relabel coefficients
    keep <- grep("^b_", colnames(post), value = TRUE)
    df_long <- as.data.frame(post) |>
      dplyr::select(all_of(keep)) |>
      dplyr::rename(
        Intercept          = tidyselect::all_of("b_Intercept"),
        `Precip Effect`    = tidyselect::all_of("b_mprec"),
        `Non-focal Effect` = tidyselect::all_of("b_nfplant"),
        `DBH Effect`       = tidyselect::all_of("b_dbh")
      )

    # Burned coefficient: try common names and keep if present
    burn_names <- c("b_burnedBurned","b_burned1","b_burnedTRUE")
    burn_present <- intersect(burn_names, colnames(post))
    if (length(burn_present) == 1) {
      df_long <- df_long |>
        dplyr::mutate(`Burned Effect` = .data[[burn_present]])
    }

    # Long format
    df_long <- df_long |>
      tidyr::pivot_longer(everything(), names_to = "Parameter", values_to = "Posterior")

    # (Optional) thin points so the scatter isn’t too dense
    set.seed(1)
    df_long <- df_long |>
      dplyr::group_by(Parameter) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::filter(row_id %in% sample(row_id, min(1500, dplyr::n()), replace = FALSE)) |>
      dplyr::ungroup() |>
      dplyr::select(-row_id)


    # Plot
    ggplot(df_long, aes(x = Parameter, y = Posterior)) +
      geom_boxplot(fill = "skyblue", alpha = 0.5, width = 0.55, outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.18, height = 0),
                 alpha = 0.45, size = 1.1, color = "blue") +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(
        title = "Posterior Distribution for Caterpillar Model",
        x = "Parameters",
        y = "Posterior Samples"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 28, hjust = 1))


    #Intercept → baseline predicted value when all predictors are 0 (or at their reference level).
    #b_burnedBurned / Burned Effect → how much the “Burned” category differs from “No burn.”
    #DBH Effect → effect of tree diameter on the outcome.
    #Non-focal Effect → effect of surrounding plant density.
    #Precip Effect → effect of precipitation.


    #Violin plot 
    ggplot(df_long, aes(x = Parameter, y = Posterior, fill = Parameter)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = "Posterior Distributions (Violin View)",
           x = "Parameter", y = "Posterior Samples") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")

    #For intercept: When environmental covariates (DBH, precipitation, and non-focal vegetation) 
    #are at their typical values and sites are unburned, the expected caterpillar abundance is 
    #modestly positive, suggesting a consistent baseline presence of caterpillars across study sites.




    ####Interactions 
    fit <- brm(
      cat ~ winter * nfplant * dbh * burned,
      data = CDC25,
      family = gaussian(),
      cores = 4, chains = 4, iter = 2000, seed = 123
    )


    library(bayesplot)
    library(dplyr)
    library(tidyr)
    library(ggplot2)

    post <- as.matrix(fit)

    # Grab and relabel coefficients
    keep <- grep("^b_", colnames(post), value = TRUE)
    df_long <- as.data.frame(post) |>
      dplyr::select(all_of(keep)) |>
      dplyr::rename(
        Intercept          = tidyselect::all_of("b_Intercept"),
        `Precip Effect`    = tidyselect::all_of("b_winter"),
        `Non-focal Effect` = tidyselect::all_of("b_nfplant"),
        `DBH Effect`       = tidyselect::all_of("b_dbh")
      )

    # Burned coefficient: try common names and keep if present
    burn_names <- c("b_burnedBurned","b_burned1","b_burnedTRUE")
    burn_present <- intersect(burn_names, colnames(post))
    if (length(burn_present) == 1) {
      df_long <- df_long |>
        dplyr::mutate(`Burned Effect` = .data[[burn_present]])
    }

    # Long format
    df_long <- df_long |>
      tidyr::pivot_longer(everything(), names_to = "Parameter", values_to = "Posterior")

    # (Optional) thin points so the scatter isn’t too dense
    set.seed(1)
    df_long <- df_long |>
      dplyr::group_by(Parameter) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::filter(row_id %in% sample(row_id, min(1500, dplyr::n()), replace = FALSE)) |>
      dplyr::ungroup() |>
      dplyr::select(-row_id)


    #Violin plot 
    ggplot(df_long, aes(x = Parameter, y = Posterior, fill = Parameter)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = "Posterior Distributions (Violin View)",
           x = "Parameter", y = "Posterior Samples") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")





    library(posterior)
    library(tidybayes)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(forcats)
    library(ggplot2)

    install.packages("tidybayes")  # run this once if not installed
    library(tidybayes)


    # 1) Get draws for all regression coefficients (b_* terms)
    draws <- as_draws_df(fit)

    coefs <- draws %>%
      dplyr::select(starts_with("b_")) %>%
      pivot_longer(everything(), names_to = "term", values_to = "value") %>%
      # 2) Classify by interaction order and make pretty labels
      mutate(
        order = stringr::str_count(term, ":"),                 # 0=main,1=2-way,2=3-way,3=4-way
        group = dplyr::recode(order,
                              `0` = "Main effects",
                              `1` = "Two-way interactions",
                              `2` = "Three-way interactions",
                              `3` = "Four-way interaction"),
        pretty = term %>%
          str_remove("^b_") %>%
          str_replace_all(":", " × ") %>%                      # nicer interaction glyph
          str_replace("burnedBurned", "Burned")
      )

    # 3) Summaries (50% & 95% CrIs) + "excludes 0?" flag
    sum_df <- coefs %>%
      group_by(term, pretty, group) %>%
      median_qi(value, .width = c(.95, .50)) %>%   # now available
      ungroup() %>%
      group_by(term, pretty, group) %>%
      mutate(excludes0 = any(.width == 0.95 & (.lower > 0 | .upper < 0))) %>%
      ungroup()


    # 4) Keep what you want to show (e.g., Main + Two-way); drop the intercept if desired
    sum_show <- sum_df %>%
      filter(group %in% c("Main effects", "Two-way interactions")) %>%
      filter(pretty != "Intercept")

    # 5) Order terms within facet by median


    library(dplyr)
    library(forcats)

    # Compute one ordering value per (group, pretty) from the 50% interval row,
    # then reorder factors within each facet using that value.
    sum_show <- sum_show %>%
      group_by(group, pretty) %>%
      mutate(order_val = median(value[.width == 0.50], na.rm = TRUE)) %>%
      ungroup() %>%
      group_by(group) %>%
      mutate(pretty = forcats::fct_reorder(pretty, order_val)) %>%
      ungroup()

    # 6) Plot: forest-style, faceted, with 95% and 50% intervals
    ggplot(sum_show, aes(x = value, y = pretty)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_errorbarh(data = ~ dplyr::filter(.x, .width == 0.95),
                     aes(xmin = .lower, xmax = .upper, color = excludes0),
                     height = 0.25, size = 0.8) +
      geom_errorbarh(data = ~ dplyr::filter(.x, .width == 0.50),
                     aes(xmin = .lower, xmax = .upper),
                     height = 0.25, size = 1.4) +
      geom_point(data = ~ dplyr::filter(.x, .width == 0.50), size = 1.6) +
      facet_wrap(~ group, scales = "free_y", ncol = 1) +
      scale_color_manual(values = c(`TRUE` = "steelblue4", `FALSE` = "grey55"),
                         name = "95% CrI excludes 0?") +
      labs(title = "Posterior coefficients by order of interaction",
           x = "Coefficient (median with 50% and 95% CrIs)",
           y = NULL) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")

    #The model suggests that caterpillar abundance is shaped more by subtle habitat structure 
    #(tree size, surrounding vegetation) than by burn history or precipitation alone. None of the 
    #predictors show high certainty, but DBH and vegetative density consistently lean positive. 
    #Interactions indicate that structural habitat features (DBH × vegetation) may modestly reinforce 
    #each other, while climatic and disturbance variables (winter × burned) show little synergy.
