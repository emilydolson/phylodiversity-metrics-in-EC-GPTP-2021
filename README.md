# What can phylogenetic metrics tell us about useful diversity in evolutionary algorithms?

This repository contains all code and supplemental material for "What can phylogenetic metrics tell us about useful diversity in evolutionary algorithms?", originally presented at [Genetic Programming in Theory and Practice](http://gptp-workshop.com/schedule.html), 2021. This paper is available as a [preprint](https://arxiv.org/abs/2108.12586) and in the [published conference proceedings](https://link.springer.com/chapter/10.1007/978-981-16-8113-4_4).

[![supplemental](https://img.shields.io/badge/go%20to-supplemental%20material-ff69b4)](https://emilydolson.github.io/phylodiversity-metrics-in-EC-GPTP-2021)
[![preprint](https://img.shields.io/badge/preprint-arXiv:2108.12586-brightgreen)](https://arxiv.org/abs/2108.12586)
[![paper](https://img.shields.io/badge/published%20in-GPTP%202021-yellow)](https://link.springer.com/chapter/10.1007/978-981-16-8113-4_4)
[![OSF](https://img.shields.io/badge/data%20%40%20OSF-https%3A%2F%2Fosf.io%2F6rndg%2F-blue)](https://osf.io/6rndg/) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5456710.svg)](https://doi.org/10.5281/zenodo.5456710)

  **Contents:**
  - [Abstract](https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021#abstract)
  - [Supplemental information](https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021#supplemental-information)
  - [Dependencies](https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021#dependencies)
  - [Compilation](https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021#compilation)
  - [Authors](https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021#authors)
  - [Research overview](https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021#research-overview)
  - [Repository contents](https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021#repository-contents)

## Abstract

> It is generally accepted that "diversity" is associated with success in evolutionary algorithms. However, diversity is a broad concept that can be measured and defined in a multitude of ways. To date, most evolutionary computation research has measured diversity using the richness and/or evenness of a particular genotypic or phenotypic property. While these metrics are informative, we hypothesize that other diversity metrics are more strongly predictive of success. Phylogenetic diversity metrics are a class of metrics popularly used in biology, which take into account the evolutionary history of a population. Here, we investigate the extent to which 1) these metrics provide different information than those traditionally used in evolutionary computation, and 2) these metrics better predict the long-term success of a run of evolutionary computation. We find that, in most cases, phylogenetic metrics behave meaningfully differently from other diversity metrics. Moreover, our results suggest that phylogenetic diversity is indeed a better predictor of success.

## Supplemental information

The supplemental information for this paper is [here](https://emilydolson.github.io/phylodiversity-metrics-in-EC-GPTP-2021).

## Dependencies

The C++ code to run these experiments requires:
- [Empirical](https://github.com/devosoft/Empirical)
- The [EC Ecology toolbox](https://github.com/emilydolson/ec_ecology_toolbox)

## Compilation

You can compile and run the code used in the paper as follows:

```{bash, shell_installation}
# Clone Empirical
git clone --recursive https://github.com/devosoft/Empirical.git 

# Clone EC-ecology-toolbox
git clone https://github.com/emilydolson/ec_ecology_toolbox.git

# Clone the repo for this project
git clone --recursive https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021.git

### Complex fitness landscapes

# Compile the executable to run experiments for this project
make

# Run an experiment. To set parameters, use command line flags
# e.g. to set the selection scheme, run ./ecology_parameter_sweep -SELECTION 2
# To see all options, run ./ecology_parameter_sweep --help
./ecology_parameter_sweep

### Exploration diagnostic

# all of the code for the exploration diagnostic lives in the exploration_diagnostic submodule
cd exploration_diagnostic
make
./dia_world

# the dia_world executable can be configured in the same way as the ecology_parameter_sweep executable

```

## Authors

- [Jose Guadalupe Hernandez](https://jgh9094.github.io/)
- [Alexander Lalejini](https://lalejini.com/)
- [Emily Dolson](http://emilyldolson.com/)

## Research overview

### Phenotypic diversity vs phylogenetic diversity.

In short, phenotypic diversity measures the diversity of [phenotypes](https://stackoverflow.com/questions/30002900/definitions-of-phenotype-and-genotype/30005949#30005949) in the population at any one point in time. Phylogenetic diveristy measures the diversity of evolutionary history represented in a population. We wrote a lot more about building phylogenies in the context of computational evolution [in this paper](https://github.com/emilydolson/interpreting_the_tape_of_life#metricvisualization-implementations).

As an example, the following figure shows two different phylogenies (ancestry trees). Arrows show parent-child relationships. Each node is a taxonomically unique phenotype (i.e., a phenotype with a unique evolutionary origin). For simplicity, leaf nodes in these diagrams are assumed to be the current set of taxa in the population; in reality, there could be non-leaf nodes corresponding to extant taxa. A) A population with high phenotypic diversity (phenotypic richness = 5) and low phylogenetic diversity (mean pairwise distance = 2). B) A population with low phenotypic diversity (phenotypic richness = 2) and high phylogenetic diversity (mean pairwise distance = 6).

![Example of populations with different levels of phenotypic and phylogenetic diveristy](
conceptual_fig.png)

### Research questions

1. Is phylogenetic diversity meaningfully different from phenotypic diversity in the context of evolutionary computation?

The answer to this question is important. Intuitively, we might think that since these are both types of diversity, they should correlate pretty closely. Given that phylogenetic diversity is more computationally intensive to measure, if we're going to argue that it's something evolutionary computation researchers should pay attention to (spoilers: we are!), we need to show that it is meaningfully different.

2. Is phylogenetic diversity more informative about outcomes in evolutionary computation than phenotypic diversity?

The importance of this question is more obvious. We know that diversity is centrally linked to the success of evolutionary algorithms. There are hints scattered across the literature that certain types of diversity are more "useful" to solving problems than others. So our goal is for this work to move us towards a better understanding of which types of diversity we should be promoting in evolutionary algorithms.

### Study design

We ran 5 selection schemes (random, tournament, fitness sharing, lexicase selection, and eco-ea) on 5 different problems (one designed to be a clean test environment, and 4 chosen to evoke the messy realities of real problems) and gathered a ludicrous amount of data. Here and in the paper, we attempt to focus very closely on getting answers to the two specific questions that we asked above (to avoid overwhelming ourselves or the reader with a firehose of data). There are many intriguing aspects of this data set that raise further questions, which we look forward to addressing in the future.

### Results

1. Phylogenetic diversity and phenotypic diversity behave differently to an extent that was even surprising to us. 

2. Phylogenetic diversity is more predictive of success than phenotypic diversity in the vast majority of cases. The differences are often substantial (check out our effect sizes!).

### Caveats/areas for future research

- Phylogenetic diversity and phenotypic diversity are both broad classes of metrics, and there is substantial variation in how different phylogenetic diversity metrics behave in different contexts.

- There is clearly variation in all of this over time and by fitness landscape.

## Measuring phylogenies

If you want to start measuring phylogenetic diversity in your own research, check out the following phylogeny-tracking libraries:

- [phylotrackpy](https://github.com/emilydolson/phylotrackpy/) (python)
- [Empirical phylogeny tracker](https://empirical--466.org.readthedocs.build/en/466/library/Evolve/systematics.html) (C++) - part of a larger library, but the library is header-only so you can just use the systematics manager module. Phylotrackpy is just a wrapper around this code and is currently better documented, so we currently reccomend reading the phylotrackpy documentation to even if you're using the C++ implementation

## Repository contents

- **.github/workflows**: Github actions to automatically build supplemental material
- **analysis**: Contains Rmarkdown (and compiled html) for all code used to analyze data for this project. This is broken into analysis for the exploration diagnostic and analysis for the other four problems. These files are compiled into the supplemental material
- **book**: The supplemental material (this is a pre-compiled backup - the full supplemental material is auto-generated via github actions)
- **exploration_diagnostic**: A submodule containing all code for the exploration diagnostic experiments.
- **scripts**: Contains all scripts used to wrangle experiments and data
  - count_convergence.py: A script to count the number of times a phylogeny rediscovers the same phenotype
  - de_dup_phylogeny.py: A script to remove duplicate lines from phylogenies (can happen when phylogeny files from multiple time points are concatenated)
  - generate_hpcc_jobs.py: The script used to submit experiments for this paper to our computing cluster
  - munge_data.py: The script used to pull all of the data into a single file.
- **source**: Contains the code used to generate the executable that was run in all complex fitness landscape experiments:
  - ecology_world.hpp: Contains the bulk of the code for the experiment
  - BitSorterMutators.hpp: Handles mutations for the sorting networks problem
  - TaskSet.hpp: Used to build the Logic-9 problem
  - TestcaseSet.hpp: Used for the program synthesis problem
  - org_types.hpp: Used to help wrangle the different representations used for each fitness landscape
  - **native**: Contains the .cc file that is ultimately compiled.
- **testcases**: Contains test cases (from the [program synthesis benchmark suite](https://thelmuth.github.io/GECCO_2015_Benchmarks_Materials/)) used for the program synthesis problems.
- DESCRIPION: For keeping tracking of R dependencies
- LICENSE: The MIT License
- Makefile: contains instructions for compiling the executable
- \_bookdown.yml: For the supplemental material
- \output.yml: For the supplemental material
- book.bib: For the supplemental material
- build_book.sh: Builds the supplemental material
- index.Rmd: Contains the beginning of the supplemental material
- packages.bib: Citations for all the R packages we used
- style.css: Styling for the supplemental material
- tail.Rmd: For the supplemental material
