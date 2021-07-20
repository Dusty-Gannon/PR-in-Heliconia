Testing for generality of pollinator recognition in *Heliconia*:
Analysis of pollen export
================
D.G. Gannon, A.S. Hadley, U.G. Kormann, F.A. Jones, M.G. Betts

``` r
knitr::opts_chunk$set(echo = TRUE)

# Packages
  require(tidyverse)
```

    ## Loading required package: tidyverse

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.2     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.3     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
  require(here)
```

    ## Loading required package: here

    ## here() starts at /Users/dusty/Documents/Heliconia/PR-in-Heliconia

## Load data

``` r
# aviary experiment data
  av <- read_csv(file = here(
    "Data",
    "aviaries_data.csv"
  ))
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
# Load pollen slides data
  pol <- read_csv(file = here(
    "Data",
    "pollen_slides.csv"
  ))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Flower_ID = col_character(),
    ##   Count = col_double()
    ## )

``` r
# pollen counter included blank rows
  pol <- pol[complete.cases(pol), ]
  
# merge the two datasets
  
  pol_all <- merge(
    pol,
    av
  )
  
# now add back data from failed RNAseq experiments
## all these visits were by GREH
  
  pol_all_sub <- pol_all[,which(
    names(pol_all) %in% c("Flower_ID", "Count", "Species", "Treatment")
  )]
  
  rna <- pol[str_which(pol$Flower_ID, "R"), ]
  rna <- rna %>% mutate(
    Species = rep("H.tortuosa", nrow(rna)),
    Treatment = rep("GREH", nrow(rna))
  )
  
  pol_all_sub <- rbind(pol_all_sub, rna)
```

## Summarize data

``` r
# summarize for tortuosa
  group_by(pol_all_sub, Species, Treatment) %>%
    summarise(
      mean_count = mean(Count),
      sd_cound = sd(Count),
      n=n()
    )
```

    ## `summarise()` has grouped output by 'Species'. You can override using the `.groups` argument.

    ## # A tibble: 4 x 5
    ## # Groups:   Species [3]
    ##   Species      Treatment mean_count sd_cound     n
    ##   <chr>        <chr>          <dbl>    <dbl> <int>
    ## 1 H.rostrata   RTAH            202      163.     3
    ## 2 H.tortuosa   GREH           1222.     696.    17
    ## 3 H.tortuosa   RTAH            125      211.     3
    ## 4 H.wagneriana RTAH              0       NA      1

``` r
# Just out of curiosity

  heto_pol <- subset(pol_all_sub, Species=="H.tortuosa")
  
  t.test(
    x=heto_pol$Count[heto_pol$Treatment=="RTAH"],
    y=heto_pol$Count[heto_pol$Treatment=="GREH"]
  )
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  heto_pol$Count[heto_pol$Treatment == "RTAH"] and heto_pol$Count[heto_pol$Treatment == "GREH"]
    ## t = -5.2678, df = 11.643, p-value = 0.0002195
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -1552.0260  -641.6211
    ## sample estimates:
    ## mean of x mean of y 
    ##   125.000  1221.824
