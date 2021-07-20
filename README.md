# Testing for pollinator recognition in *Heliconia*

Partner redundancy increases the robustness of plant-pollinator communities to extinctions of mutualistic partners. However, selection may reinforce floral traits that filter pollinator communities to promote pollination by efficient pollinators, reducing redundancy. This repository documents the analysis of data collected from experiments designed to test the generality of a recently described, cryptic pollinator filter termed ‘pollinator recognition’ (PR) by [Betts et al. (2015)](https://www.researchgate.net/publication/273156439_Pollinator_recognition_by_a_keystone_tropical_plant). Pollinator recognition was first documented in *Heliconia tortuosa* on the grounds that with pollen quality experimentally held constant, pollen tube germination occurred following visits from some hummingbird species but not others. Such a cryptic pollinator filter  could drastically reduce the realized number of pollinators compared to the number of floral visitors. We assessed the prevalence of PR in four taxa spread widely across the Heliconiaceae. 

With aviary experiments that standardized pollen quality and minimized variation in pollen quantity, we corroborated previous results that visits to *H. tortuosa* flowers by pollen-free hummingbirds with bill morphologies that match the flower morphology increased pollen tube rates compared to visits by morphologically mismatched hummingbirds as well as hand pollination alone. In *H. rostrata* we observed an increase in pollen tube rates in flowers visited by hummingbirds compared to hand pollination alone, regardless of the bill morphology. In two other species (*H. hirsuta* and *H. wagneriana*), we observed decreased pollen tube rates in flowers visited by hummingbirds, potentially due to pollen desiccation and removal by visiting hummingbirds. Finally, we were unable to corroborate previous results that nectar removal by morphologically matched hummingbird provides the mechanism of recognition. Our results highlight the complexities and variability in plant mating strategies and the need for further study of pollination of *Heliconia* plants.


### Repository organization

#### Data folder

All necessary data files can be found in the [Data](https://github.com/Dusty-Gannon/PR-in-Heliconia/tree/main/Data) folder.

* The file `aviaries_data.csv` includes all data from the single-vist aviary experiments necessary to reproduce the results in the manuscript.

  - Flower_ID: Individual identifier for the observation
  - Experiment: Experiment identifier
  - Species: Plant species (*H. tortuosa*, *H. rostrata*, *H. hirsuta*, or *H. wagneriana*)
  - Focal_Loc: General location of plant
  - Plant: Unique plant ID
  - Treatment: Experimental treatment. **HP** - hand pollination only; **RTAH** - hand pollination followed by a visit from a pollen-free Rufous-Tailed Hummingbird (short bill); **GREH** - hand pollination followed by a visit from a pollen-free Green Hermit Hummingbird (long, curved bill)
  - HP_time: Time of hand pollination
  - PD_Loc: General location of pollen donor
  - PD_plant: Unique identifier for pollen donor plant
  - Date: Date of experiment
  - Tube_Count: Number of pollen tubes found in the style

* The file `nectar_experiments.csv` includes all data from experiments designed to test for the mechanism of pollinator recognition. All field names are the same as above, but the treatments are as follows: **HP** - hand pollination only; **HPNE** - hand pollination followed by manual nectar removal using a pipette tip and syringe; **PM** - hand pollination followed by the insertion of the pipette tip without any attempt to remove nectar from the flower.

* The file `bird_bill_morph.csv` along with `heliconia_morphology.csv` contain data necessary to reproduce Fig 1 of the manuscript which compares the morphology the bird bills to flower shapes (see the [trait matching analysis]())

* The file `pollen_slides.csv` contains data on pollen removal after single visits to *H. tortuosa* flowers by different birds.

* The file `HPvsOP_all.csv` includes all the same fields at the aviary data (see above), but the treatments are either hand pollination only (**HP**) or open pollination (**OP**).


#### Analyses folder

RMarkdown documents with annotated analysis code can be found in the [Analyses](https://github.com/Dusty-Gannon/PR-in-Heliconia/tree/main/Analyses) folder. The `.md` are knitted versions of the `.Rmd` files that should render on GitHub. The results from the manuscript should be fully reproducible. 

Any file ending `_bayes.md` includes a Bayesian analysis with a more complicated variance structure than that used in the manuscript. For ease of comparison to [Betts et al. (2015)](https://www.researchgate.net/publication/273156439_Pollinator_recognition_by_a_keystone_tropical_plant), we fitted frequentist GLMMs to the data for the manuscript (files ending in `_freq.md`).

#### MS_files folder

A detailed description of the experimental, laboratory, and statistical methods can be found in the [supplementary methods](https://github.com/Dusty-Gannon/PR-in-Heliconia/blob/main/MS_files/Gannon_et_al_Appendix_S1.pdf) and in the [draft manuscript](https://github.com/Dusty-Gannon/PR-in-Heliconia/tree/main/MS_files).






