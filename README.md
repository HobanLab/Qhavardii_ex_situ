<h1>Quercus havardii- genetic diversity in ex situ collections</h1>
Raw data and R scripts for analyzing Quercus havardii genetic diversity in both in situ populations and ex situ collections, and for producing Figures for papers in preparation.
<h2>Overview</h2>
It is important to conserve genetic diversity in botanic gardens, to ensure their resilience into the future.  To ensure this we can measure genetic diversity directly or measure the breadth of geography and ecological diversity conserved.  The number of samples and how to collect them in the field (e.g. number of populations to visit) in order to ensure genetic diversity conservation, is not known.  Recent work has looked at minimum sample sizes to conserve rare species, but there are few if any guidelines for seed sampling for widespread but locally rare and threatened species.  We investigate this in the IUCN Red List Threatened species Quercus havardii, a desert adapted oak found in the Western United States in two disjunct regions.  We quantify how much genetic diversity is conserved in each of numerous botanic gardens, and we resample the in situ genetic data to determine minimum sample sizes with 'ideal' sampling according to best practice guidelines. 

<h1>Development</h1>

We welcome contributions from any individual, whether code, documentation, or issue tracking. All participants are expected to follow the code of conduct for this project.

Authors and Contributions:
*Bethany A. Zumwalde (BAZ)- co-led the writing, provided mentoring and advising, performed several analyses, led revisions after peer review
*Bailie Fredlock Munoz (BFM)- performed lab work, wrote first draft of manuscript
*Ross A. McCauley (RM)- field sampling and discussions
*Drew Duckett (DD)- field sampling, contributed to lab work, discussions, assisted with writing
*Emma Suzuki Spence (ESS)- performed some lab work, provided mentoring and advising
*Emily Beckman Bruns (EBB)- performed the geospatial analysis, contributed to writing and revisions
*Sean Hoban (SH)- project lead, field work, performed analysis and created figures, co-led the writing and revisions

<h1> Data Collection</h1>
Tissue for DNA analysis was sampled from Q havardii in 2016 in two separate trips referred to as East and West. The East sampling was performed by Sean and Drew in Texas, New Mexico, and Oklahoma. The West sampling was performed by Sean and Ross McCauley (Fort Lewis College) in Utah and Arizona with some additional samples collected in Colorado and New Mexico. Samples were also provided by collaborators. A total of 667 samples from 26 primary populations and 10 auxiliary populations of georeferenced locations were used in this study. Populations were chosen by contacting land managers of private and public land in the region, by consulting GBIF and SEINet, and via suggestions from the International Oak Society. The objective was to ensure populations were sampled throughout the geographic range. In this same sampling trip, we collected seeds for ex situ conservation. In the following years, 2017 and 2018, we sampled leaf tissue from these seedlings for DNA analysis. Leaves from 290 seedlings from 66 maternal trees from 26 populations were sampled.

<h1> File explanations</h1>
All .gen files are genpop files that contain the genetic data.  The differences among the .gen files are just different subsets of the data or different ways of combining populations, described as follows:
*QH_total_garden_by_pop.gen contains all the genotypes ex situ and in situ. Each population is separated by POP. The East populations come first, and the West populations come next, and then the garden populations
*QH_total_garden_by_pop_E.gen contains all the genotypes from the in situ populations and all the ex situ samples from East populations. Each population is separated by POP. The East populations come first, and the West populations come next, and then the garden populations
*QH_total_garden_by_pop_W.gen contains all the genotypes from the in situ populations and all the ex situ samples from West populations. Each population is separated by POP. The East populations come first, and the West populations come next, and then the garden populations
*QH_total_garden_for_FST.gen contains all the genotypes ex situ and in situ. The only difference from "garden_by_pop" is that all the East populations are merged under one POP and all the West populations are merged under one POP. Each garden population is separated by POP. (The purpose of this was to calcualte FST between each garden and the East and West.) As with the others, the East populations come first, and the West populations come next, and then the garden populations
There are also several .txt files
*naming_samples_by_gard.txt lists all the garden samples, and which populations, regions and maternal trees they come from
*reduced_prop_capture_lm.txt is the proportion of alleles captured in each allele category by each botanic garden (reduced means only considering alleles with >2 copies in the dataset)
*all_prop_capture_lm.txt is the proportion of alleles captured in each allele category by each botanic garden (all means considering all alleles no matter the number of copies)

<h1> Scripts for analysis</h1>
There is one main long analysis script- "Qhav_ex_situ_code.R"
The first code chunk calculates how much genetic diversity exists in botanic gardens, by comparing the ex situ dataset to the in situ dataset
The second code chunk subsamples the in situ populations to represent a seed sampler using ideal sampling (random from all populations) and thus to calculate a minimum sample size
The remainder of the code creates plots for the manuscript, and performs the linear models 
