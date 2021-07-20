Testing for generality of pollinator recognition in *Heliconia*:
Frequentist analysis of aviary data
================
D.G. Gannon, A.S. Hadley, U.G. Kormann, F.A. Jones, M.G. Betts

### R packages

``` r
knitr::opts_chunk$set(echo = TRUE)

    require(tidyverse)
    require(here)
    require(lme4)
    require(gridExtra)
```

### Load data

``` r
# Load data
  av <- read_csv(file = here("Data", "aviaries_data.csv"))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Flower_ID = col_double(),
    ##   Experiment = col_double(),
    ##   Species = col_character(),
    ##   Focal_Loc = col_character(),
    ##   Plant = col_character(),
    ##   Treatment = col_character(),
    ##   HP_time = col_double(),
    ##   PD_Loc = col_character(),
    ##   PD_plant = col_character(),
    ##   Bird = col_character(),
    ##   Date = col_double(),
    ##   Tube_Count = col_double()
    ## )

``` r
# Sort rows by species, then plant individual ID, then experiment ID
  av <- av[with(av, order(Species, Plant, Experiment)), ]
  
# Create factors for species, treatment, and plant individual
  av$f.species <- as.factor(av$Species)
  av$Treatment2 <- av$Treatment

# relable bird species as a treatment
  av$Treatment2[av$Treatment2 == "RTAH"] <- "SB"
  av$Treatment2[av$Treatment2 == "GREH"] <- "LB"
  
  av$f.trtmnt <- factor(av$Treatment2, levels = c("HP", "SB", "LB"))
  av$f.plant <- factor(av$Plant)
```

## Summary

To test whether pollen germination and tube growth is dependent on the
identity and morphology of a floral visitor, we conducted pollination
experiments with captive hummingbirds inside portable aviaries. In these
experiments, we randomly assigned flowers to one of three experimental
treatments: 1) hand-pollination only (HP treatment); 2) hand pollination
followed by a visit from a pollen-free rufous-tailed hummingbird (short,
straight bills; SB treatment); 3) hand pollination followed by a visit
from a pollen-free green hermit hummingbird (long, slightly decurved
bill; LB treatment).

This experiment follows a factorial design with factors for the plant
species and the treatment, but due to limited numbers of flowering
plants, we often conducted multiple experiments using flowers (the
experimental units) from the same plant. We therefore fitted the model
described below.

### Model

Let ![y\_{ijkl} \\in
\\mathbb{N}](https://latex.codecogs.com/png.latex?y_%7Bijkl%7D%20%5Cin%20%5Cmathbb%7BN%7D
"y_{ijkl} \\in \\mathbb{N}") be the number of pollen tubes scored in the
![l^\\text{th}](https://latex.codecogs.com/png.latex?l%5E%5Ctext%7Bth%7D
"l^\\text{th}") flower from the
![k^\\text{th}](https://latex.codecogs.com/png.latex?k%5E%5Ctext%7Bth%7D
"k^\\text{th}") plant of species
![j=1,..,4](https://latex.codecogs.com/png.latex?j%3D1%2C..%2C4
"j=1,..,4") and ![i](https://latex.codecogs.com/png.latex?i "i") index
the treatment
(![i=1,2,3](https://latex.codecogs.com/png.latex?i%3D1%2C2%2C3
"i=1,2,3")). We assume that ![y\_{ijkl} \\overset{iid}{\\sim}
\\text{Poisson}(\\lambda\_{ijk})](https://latex.codecogs.com/png.latex?y_%7Bijkl%7D%20%5Coverset%7Biid%7D%7B%5Csim%7D%20%5Ctext%7BPoisson%7D%28%5Clambda_%7Bijk%7D%29
"y_{ijkl} \\overset{iid}{\\sim} \\text{Poisson}(\\lambda_{ijk})") for
![l=1,...,n\_k](https://latex.codecogs.com/png.latex?l%3D1%2C...%2Cn_k
"l=1,...,n_k"). The model for the experiment can be written as

  
![&#10;\\log(\\lambda\_{ijk}) = \\mu + \\alpha\_i + \\beta\_j +
(\\alpha\\beta)\_{ij} +
\\gamma\_{k(j)}&#10;](https://latex.codecogs.com/png.latex?%0A%5Clog%28%5Clambda_%7Bijk%7D%29%20%3D%20%5Cmu%20%2B%20%5Calpha_i%20%2B%20%5Cbeta_j%20%2B%20%28%5Calpha%5Cbeta%29_%7Bij%7D%20%2B%20%5Cgamma_%7Bk%28j%29%7D%0A
"
\\log(\\lambda_{ijk}) = \\mu + \\alpha_i + \\beta_j + (\\alpha\\beta)_{ij} + \\gamma_{k(j)}
")  

where

  - ![\\mu](https://latex.codecogs.com/png.latex?%5Cmu "\\mu") is the
    overall mean log-pollen tube rate,

  - ![\\alpha\_i,\\ \\sum\_{i=1}^3
    \\alpha\_i=0](https://latex.codecogs.com/png.latex?%5Calpha_i%2C%5C%20%5Csum_%7Bi%3D1%7D%5E3%20%5Calpha_i%3D0
    "\\alpha_i,\\ \\sum_{i=1}^3 \\alpha_i=0"), is the average deviation
    from the mean for flowers that received treatment
    ![i](https://latex.codecogs.com/png.latex?i "i"),

  - ![\\beta\_j, \\sum\_{j=1}^4 \\beta\_j
    = 0](https://latex.codecogs.com/png.latex?%5Cbeta_j%2C%20%5Csum_%7Bj%3D1%7D%5E4%20%5Cbeta_j%20%3D%200
    "\\beta_j, \\sum_{j=1}^4 \\beta_j = 0"), is the average deviation
    from the mean for flowers of plant species
    ![j](https://latex.codecogs.com/png.latex?j "j"),

  - ![(\\alpha\\beta)\_{ij}, \\sum\_{i=1}^3 (\\alpha\\beta)\_{ij} =
    \\sum\_{j=1}^4 (\\alpha\\beta)\_{ij}
    = 0](https://latex.codecogs.com/png.latex?%28%5Calpha%5Cbeta%29_%7Bij%7D%2C%20%5Csum_%7Bi%3D1%7D%5E3%20%28%5Calpha%5Cbeta%29_%7Bij%7D%20%3D%20%5Csum_%7Bj%3D1%7D%5E4%20%28%5Calpha%5Cbeta%29_%7Bij%7D%20%3D%200
    "(\\alpha\\beta)_{ij}, \\sum_{i=1}^3 (\\alpha\\beta)_{ij} = \\sum_{j=1}^4 (\\alpha\\beta)_{ij} = 0"),
    is the species\*treatment interaction effect,

  - ![\\gamma\_{k(j)} \\overset{iid}{\\sim} \\mathcal{N}(0,
    \\sigma\_\\gamma^2)](https://latex.codecogs.com/png.latex?%5Cgamma_%7Bk%28j%29%7D%20%5Coverset%7Biid%7D%7B%5Csim%7D%20%5Cmathcal%7BN%7D%280%2C%20%5Csigma_%5Cgamma%5E2%29
    "\\gamma_{k(j)} \\overset{iid}{\\sim} \\mathcal{N}(0, \\sigma_\\gamma^2)"),
    for all plants
    ![k=1,2,...,K\_j](https://latex.codecogs.com/png.latex?k%3D1%2C2%2C...%2CK_j
    "k=1,2,...,K_j") of all species
    ![j=1,2,3,4](https://latex.codecogs.com/png.latex?j%3D1%2C2%2C3%2C4
    "j=1,2,3,4"), is a random effect for the plant that is nested within
    species.

For ease in defining the model in standard R packages (`lme4` (Bates et
al. 2015)), we reparameterized the model to follow a regression
parameterization with the reference level being hand pollinated flowers
of *Heliconia hirsuta*. The model therefore reads

  
![&#10;\\log{\\boldsymbol \\lambda} = {\\bf X}{\\boldsymbol \\beta} +
{\\bf
Z}{\\boldsymbol\\gamma},&#10;](https://latex.codecogs.com/png.latex?%0A%5Clog%7B%5Cboldsymbol%20%5Clambda%7D%20%3D%20%7B%5Cbf%20X%7D%7B%5Cboldsymbol%20%5Cbeta%7D%20%2B%20%7B%5Cbf%20Z%7D%7B%5Cboldsymbol%5Cgamma%7D%2C%0A
"
\\log{\\boldsymbol \\lambda} = {\\bf X}{\\boldsymbol \\beta} + {\\bf Z}{\\boldsymbol\\gamma},
")  

where

  - ![\\underset{\\scriptsize n\\times p}{\\bf
    X}](https://latex.codecogs.com/png.latex?%5Cunderset%7B%5Cscriptsize%20n%5Ctimes%20p%7D%7B%5Cbf%20X%7D
    "\\underset{\\scriptsize n\\times p}{\\bf X}") is a matrix of
    indicator variables indicating the species, treatment, and plant
    individual for the
    ![i^\\text{th}](https://latex.codecogs.com/png.latex?i%5E%5Ctext%7Bth%7D
    "i^\\text{th}") observation,
    ![i=1,...,n](https://latex.codecogs.com/png.latex?i%3D1%2C...%2Cn
    "i=1,...,n");

  - ![\\boldsymbol
    \\beta](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cbeta
    "\\boldsymbol \\beta") is a
    ![p](https://latex.codecogs.com/png.latex?p "p")-vector of
    regression coefficients;

  - ![\\underset{\\scriptsize n\\times q}{\\bf
    Z}](https://latex.codecogs.com/png.latex?%5Cunderset%7B%5Cscriptsize%20n%5Ctimes%20q%7D%7B%5Cbf%20Z%7D
    "\\underset{\\scriptsize n\\times q}{\\bf Z}") is an indicator
    matrix indicating which plant was used in experiment
    ![i=1,2,...,n](https://latex.codecogs.com/png.latex?i%3D1%2C2%2C...%2Cn
    "i=1,2,...,n");

  - ![\\boldsymbol
    \\gamma](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cgamma
    "\\boldsymbol \\gamma") is a
    ![q](https://latex.codecogs.com/png.latex?q "q")-vector of random
    intercepts for a given plant.

## Fit the model

``` r
 fit <- glmer(
    Tube_Count ~ f.trtmnt*f.species + (1|f.species:f.plant),
    data = av,
    family = "poisson",
    control = glmerControl(optimizer = "bobyqa")
 )
```

### Model predictions

``` r
# function used to bootstrap predicted means
  predict_means <- function(mod){
    
    return(exp(X%*%mod@beta))
    
  }

# Create model matrix for predictions
  df_new <- data.frame(
    species = factor(
      rep(unique(av$Species), each=length(unique(av$Treatment))),
      levels = levels(av$f.species)
    ),
    treatment = factor(
      rep(c("HP", "SB", "LB"), length(unique(av$Species))),
      levels = levels(av$f.trtmnt)
    )
  )
  
  X <- model.matrix(~treatment*species, data = df_new)
  
# bootstrap the predictions
  
  bootpreds <- bootMer(
    fit,
    FUN = predict_means,
    nsim = 1000,
    verbose = T,
    PBargs = list(style=3)
  )
  
# save bootstrapped means for potentual future use
  # saveRDS(bootpreds, file = here(
  #   "Data",
  #   "boot_av_means.rds"
  # ))
  
# combine into new dataframe for plotting
  df_new <- df_new %>% mutate(
    pred = as.double(bootpreds$t0),
    se = as.double(
      apply(bootpreds$t, 2, sd)
    ),
    low = as.double(
      apply(bootpreds$t, 2, quantile, probs=0.025)
    ),
    high = as.double(
      apply(bootpreds$t,2, quantile, probs=0.975)
    )
  )
```

## Plotting results

**Bring in data on open pollination pollen tube rates for reference**

``` r
# read in open pollination data

 hpop <- read_csv(file = here("Data", "HPvsOP_all.csv"))

 sp <- c("H.hirsuta", "H.rostrata",
        "H.tortuosa", "H.wagneriana")

 hpop_sub <- hpop[which(hpop$Treatment == "OP" & 
                         hpop$Species %in% sp), ]

# come up with a grand mean

hpop_sub_grp <- group_by(hpop_sub, Species)

# summarize the tube rates in each species
 op_means_sp_pl <- hpop_sub_grp %>% 
   summarise(.,
            nplants=length(unique(Plant)),
            tube_rate=mean(Tube_Count),
            sd=sd(Tube_Count)
            )

# add se
 op_means_sp_pl <- op_means_sp_pl %>%
   mutate(
     se = sd/sqrt(nplants)
   )
```

### Code used to create Fig 2 in the MS

``` r
  mean_theme <- theme(panel.background = element_blank(),
                      axis.line = element_line(colour = "darkgrey"),
                      legend.position = "none",
                      axis.text = element_text(size = 12, color="black"),
                      plot.title = element_text(face = "italic",
                                                size=16))
  
  hir <- ggplot(data = df_new[1:3,], aes(x=treatment, y=pred))+
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[1],
      color="darkgrey", 
      linetype="longdash"
    ) +
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[1] + 
        op_means_sp_pl$se[1],
      color="grey",
      linetype="dashed",
      size=0.5
    ) +
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[1] - 
        op_means_sp_pl$se[1],
      color="grey",
      linetype="dashed",
      size=0.5
    ) +    
    geom_errorbar(aes(
      ymin=pred-se,
      ymax=pred+se
    ), width=0)+
    geom_point(size=4, color="white")+
    geom_point(size=3, color="black")+
    ylim(c(0,2))+
    mean_theme+
    xlab("")+
    ylab("")+
    ggtitle("H. hirsuta")
  
  
  ros <- ggplot(data = df_new[4:6,], aes(x=treatment, y=pred))+
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[2],
      color="darkgrey", 
      linetype="longdash"
    ) +
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[2] + 
        op_means_sp_pl$se[2],
      color="grey",
      linetype="dashed",
      size=0.5
    ) +
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[2] - 
        op_means_sp_pl$se[2],
      color="grey",
      linetype="dashed",
      size=0.5
    ) +  
    geom_errorbar(aes(
      ymin=pred-se,
      ymax=pred+se
    ), width=0)+
    geom_point(size=4, color="white")+
    geom_point(size=3, color="black")+
    mean_theme+
    theme(axis.text.x=element_blank())+
    ylim(c(0,2))+
    xlab("")+
    ylab("")+
    ggtitle("H. rostrata")
  
 
  tor <- ggplot(data = df_new[7:9,], aes(x=treatment, y=pred))+
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[3],
      color="darkgrey", 
      linetype="longdash"
    ) +
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[3] + 
        op_means_sp_pl$se[3],
      color="grey",
      linetype="dashed",
      size=0.5
    ) +
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[3] - 
        op_means_sp_pl$se[3],
      color="grey",
      linetype="dashed",
      size=0.5
    ) + 
    geom_errorbar(aes(
      ymin=pred-se,
      ymax=pred+se
    ), width=0)+
    geom_point(size=4, color="white")+
    geom_point(size=3, color="black")+
    ylim(c(0,2))+
    mean_theme+
    theme(axis.text.x=element_blank())+
    xlab("")+
    ylab("")+
    ggtitle("H. tortuosa")
  
  wag <- ggplot(data = df_new[10:12,], aes(x=treatment, y=pred))+
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[4],
      color="darkgrey", 
      linetype="longdash"
    ) +
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[4] + 
        op_means_sp_pl$se[4],
      color="grey",
      linetype="dashed",
      size=0.5
    ) +
    geom_hline(
      yintercept = op_means_sp_pl$tube_rate[4] - 
        op_means_sp_pl$se[4],
      color="grey",
      linetype="dashed",
      size=0.5
    ) + 
    geom_errorbar(aes(
      ymin=pred-se,
      ymax=pred+se
    ), width=0)+
    geom_point(size=4, color="white")+
    geom_point(size=3, color="black")+
    ylim(c(0,2))+
    mean_theme+
    xlab("")+
    ylab("")+
    ggtitle("H. wagneriana")
  
 png(filename = here("Figures", "means_plot_freq.png"), units = "px",
      height = 1800, width = 1800, res = 300)
   grid.arrange(tor,ros,hir,wag, nrow=2,ncol=2)
 dev.off()
```

## Various statistics for the MS

``` r
# get covariance matrix of regression coefficients
  V <- matrix(
    data = vcov(fit)@x,
    nrow = length(fit@beta),
    ncol = length(fit@beta)
  )
    
# Difference b/w hand-pollinated controls and GREH for HETO
  c_hpgreh_tor <- c(0,0,1,0,0,0,0,0,0,1,0,0)
  # mean
  exp(c_hpgreh_tor%*%fit@beta)
  # CI
  exp(c(
    fit@beta[10] - 1.96*sqrt(c_hpgreh_tor%*%V%*%c_hpgreh_tor),
    fit@beta[10] + 1.96*sqrt(c_hpgreh_tor%*%V%*%c_hpgreh_tor)
  ))
  
# Difference b/w GREH and RTAH for HETO
  c_rtah_greh_tor <- c(0,-1,1,0,0,0,0,0,-1,1,0,0)
  exp(c_rtah_greh_tor%*%fit@beta)
  # CI
  exp(c(
    c_rtah_greh_tor%*%fit@beta - 1.96*sqrt(c_rtah_greh_tor%*%V%*%c_rtah_greh_tor),
    c_rtah_greh_tor%*%fit@beta + 1.96*sqrt(c_rtah_greh_tor%*%V%*%c_rtah_greh_tor)
  ))
```

### References

<div id="refs" class="references">

<div id="ref-bates2015">

Bates, Douglas, Martin Machler, Ben Bolker, and Steve Walker. 2015.
“Fitting Linear Mixed-Effects Models Using Lme4.” *Journal of
Statistical Software* 67 (1): 1–48.
<https://doi.org/10.18637/jss.v067.i01>.

</div>

</div>
