Testing for generality of pollinator recognition in *Heliconia*:
Analysis of aviary data
================
D.G. Gannon, A.S. Hadley, U.G. Kormann, F.A. Jones, M.G. Betts

``` r
library(tidyverse)
library(here)

dist_mat <- function(X){
  G <- X%*%t(X)
  ones <- rep(1, nrow(X))
  d <- diag(G)
  D <- d%*%t(ones) - 2*G + ones%*%t(d)
  return(D)
}

euc_dist <- function(v1, v2){
  sq_diff <- (v1-v2)^2
  
  return(sqrt(sum(sq_diff)))
}
```

Here we detail the calculations for hummingbird bill and *Heliconia*
flower trait matching.

### Hummingbird bill morphology

For hummingbird bill shape, we photographed wild-caught hummingbirds’
bills against 0.25-inch gridded paper before releasing them. We then
used `ImageJ` to measure the length of the bill along the outer (top)
surface of the bill using a the line segment feature. To measure
curvature, we utilized the methods of Temeles et al. (2009) in which we
related the outer edge of the bill to a circle. The curvature of a
circle is measured as the inverse of its radius. To estimate the radius
of a circle with the same curvature as a bird’s bill, we connected two
points on the curved section of the bill, effectively drawing a chord on
the circle. We then measured the angle between the chord and the line
tangent to the bill at one end of the chord, known as the *angle of
declension*. Straight-forward trigonomentric identities then yield
![\\hat{r} =
\\frac{c/2}{\\sin{\\theta}}](https://latex.codecogs.com/png.latex?%5Chat%7Br%7D%20%3D%20%5Cfrac%7Bc%2F2%7D%7B%5Csin%7B%5Ctheta%7D%7D
"\\hat{r} = \\frac{c/2}{\\sin{\\theta}}"), where
![c](https://latex.codecogs.com/png.latex?c "c") is the length of the
chord (mm) and ![\\theta](https://latex.codecogs.com/png.latex?%5Ctheta
"\\theta") is the angle of declension (radians). Thus, we use
![\\hat{r}^{-1}](https://latex.codecogs.com/png.latex?%5Chat%7Br%7D%5E%7B-1%7D
"\\hat{r}^{-1}") as an estimate of the curvature of a bill or flower.

``` r
# Load bird data

  bills <- read_csv(
    file = here("Data", "bird_bill_morph.csv")
  )
```

    ## Parsed with column specification:
    ## cols(
    ##   Species = col_character(),
    ##   Bird = col_character(),
    ##   Date_caught = col_double(),
    ##   Location = col_character(),
    ##   bandORcolor = col_character(),
    ##   Sex = col_character(),
    ##   bill_length_mm = col_double(),
    ##   bill_chordL_mm = col_double(),
    ##   bill_chordtan_angle_deg = col_double()
    ## )

``` r
# compute curvature for bird bills
  bills <- bills %>% mutate(
    bill_chordtan_angle_rad = bill_chordtan_angle_deg*(pi/180)
  )
  
  bills <- bills %>% mutate(
    bill_curve = ((bill_chordL_mm/2)/(sin(bill_chordtan_angle_rad)))^(-1)
  )

# add column for the radius as well  
  bills <- bills %>% mutate(
    bill_r = (bill_chordL_mm/2)/(sin(bill_chordtan_angle_rad))
  )
```

### Plant morphology

``` r
  flowers <- read_csv(
    file = here("Data", "heliconia_morphology.csv")
  ) 
```

    ## Parsed with column specification:
    ## cols(
    ##   Species = col_character(),
    ##   Population = col_character(),
    ##   Plant = col_character(),
    ##   Ramet = col_double(),
    ##   Bract = col_double(),
    ##   OS_length = col_double(),
    ##   chord_length = col_double(),
    ##   total_length = col_double(),
    ##   angle = col_double(),
    ##   curvature = col_logical(),
    ##   date = col_double()
    ## )

``` r
# subset for focal species
  flowers <- flowers[
    which(flowers$Species %in%
            c("H.hirsuta", "H.rostrata",
              "H.tortuosa", "H.wagneriana")),
  ]

# compute curvature
  flowers <- flowers %>% mutate(
    angle_rad = angle*(pi/180)
  )
  
  flowers <- flowers %>% mutate(
    curvature = ((chord_length/2)/(sin(angle_rad)))^(-1)
  )
  
# compute radius as well
  flowers <- flowers %>% mutate(
    radius = (chord_length/2)/(sin(angle_rad))
  )
```

### Trait mismatch

We compute the average degree of mismatch as the distance between the
means for each species in the
![\\mathbb{R}^2](https://latex.codecogs.com/png.latex?%5Cmathbb%7BR%7D%5E2
"\\mathbb{R}^2") euclidean space with one axis as length (mm) and the
other the radius of the arc of the curve (mm).

``` r
  birds <- c("RTAH","GREH")
  plants <- c("H.hirsuta", "H.rostrata", "H.tortuosa", "H.wagneriana")
  
  mismatch <- as.data.frame(matrix(nrow = 8, ncol = 5))
  names(mismatch) <- c("bird", "plant", "mean_dist",
                       "sd_dist", "n")
  
  for(i in 1:length(birds)){
    for(j in 1:length(plants)){
      bird_i <- subset(bills, Species==birds[i])
      plant_i <- subset(flowers, Species==plants[j])
      nbirds <- nrow(bird_i)
      nflowers <- nrow(plant_i)
      N <- nbirds + nflowers
      
      # create matrix of bill measurements
      mat <- as.matrix(bird_i[,which(
        names(bird_i) %in% c("bill_length_mm", "bill_r")
      )])
      
      # bind matrix with flower measurements
      mat <- rbind(
        mat, as.matrix(
          plant_i[, which(
            names(flowers) %in% c("total_length", "radius")
          )]
        )
      )
      
      dist_ij <- dist_mat(mat)
      dists <- sqrt(as.double(dist_ij[(N-nflowers+1):N, 1:nbirds]))
      mismatch[(length(plants)*(i-1)+j), ] <-
        c(birds[i], plants[j], mean(dists), 
          sd(dists), length(dists))
    }
  }
```

### Some stats for the MS

``` r
  rtah <- subset(bills, Species=="RTAH")

  mean(rtah$bill_length_mm)
```

    ## [1] 21.60443

``` r
  sd(rtah$bill_length_mm)
```

    ## [1] 1.55278

``` r
  mean(rtah$bill_curve)
```

    ## [1] 0.01605901

``` r
  sd(rtah$bill_curve)
```

    ## [1] 0.002387717

``` r
 flwr_summary <- flowers[,
   which(names(flowers) %in% 
           c("Species", "total_length", "curvature"))
 ] %>% 
  group_by(Species) %>%
  summarise(., meanL=mean(total_length),
            s_L = sd(total_length),
            meanK=mean(curvature),
            s_K=sd(curvature))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
  flwr_summary
```

    ## # A tibble: 4 x 5
    ##   Species      meanL   s_L  meanK     s_K
    ##   <chr>        <dbl> <dbl>  <dbl>   <dbl>
    ## 1 H.hirsuta     51.6  1.31 0.0153 0.00412
    ## 2 H.rostrata    51.1  2.75 0.0184 0.00294
    ## 3 H.tortuosa    59.0  3.46 0.0265 0.00517
    ## 4 H.wagneriana  70.1  3.78 0.0169 0.00314

<div id="refs" class="references">

<div id="ref-temeles2009">

Temeles, Ethan J., Carolyn R. Koulouris, Sarah E. Sander, and W. John
Kress. 2009. “Effect of Flower Shape and Size on Foraging Performance
and Trade-Offs in a Tropical Hummingbird.” *Ecology* 90 (5): 1147–61.
<https://doi.org/10.1890/08-0695.1>.

</div>

</div>
