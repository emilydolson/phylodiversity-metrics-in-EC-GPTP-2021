# What can phylogenetic metrics tell us about useful diversity in evolutionary algorithms?

This repository contains all code and supplemental material for "What can phylogenetic metrics tell us about useful diversity in evolutionary algorithms?", originally presented at [Genetic Programming in Theory and Practice](http://gptp-workshop.com/schedule.html), 2021. This paper is currently available as a preprint, and will eventually appear in the published conference proceedings.

[![preprint](https://img.shields.io/badge/preprint-arXiv:2108.12586-brightgreen)](https://arxiv.org/abs/2108.12586)
[![OSF](https://img.shields.io/badge/data%20%40%20OSF-https%3A%2F%2Fosf.io%2F6rndg%2F-blue)](https://osf.io/6rndg/) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4733407.svg)](https://doi.org/10.5281/zenodo.4733407)


  **Contents:**
  - [Abstract](https://github.com/emilydolson/phylodiversity-metrics-in-EC-GPTP-2021#abstract)
  - Dependencies

## Abstract

> It is generally accepted that "diversity" is associated with success in evolutionary algorithms. However, diversity is a broad concept that can be measured and defined in a multitude of ways. To date, most evolutionary computation research has measured diversity using the richness and/or evenness of a particular genotypic or phenotypic property. While these metrics are informative, we hypothesize that other diversity metrics are more strongly predictive of success. Phylogenetic diversity metrics are a class of metrics popularly used in biology, which take into account the evolutionary history of a population. Here, we investigate the extent to which 1) these metrics provide different information than those traditionally used in evolutionary computation, and 2) these metrics better predict the long-term success of a run of evolutionary computation. We find that, in most cases, phylogenetic metrics behave meaningfully differently from other diversity metrics. Moreover, our results suggest that phylogenetic diversity is indeed a better predictor of success.

## Supplemental information

The supplemental information for this paper is still under construction, but the key analyses can be found here:

- [Exploration diagnostic](https://emilydolson.github.io/phylodiversity-metrics-in-EC-GPTP-2021/analysis/Exploration_diagnostic_analysis.html)
- [Other fitness landscapes](https://emilydolson.github.io/phylodiversity-metrics-in-EC-GPTP-2021/analysis/other_problems_analysis.html)

## Dependencies

The C++ code to run these experiments requires:
- [Empirical](https://github.com/devosoft/Empirical)
- The [EC Ecology toolbox](https://github.com/emilydolson/ec_ecology_toolbox)
