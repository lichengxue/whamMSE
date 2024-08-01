# Management Strategy Evaluation Package for WHAM

The whamMSE package is a comprehensive toolbox designed to conduct management strategy evaluation (MSE) for supporting spatial management for complex, heterogeneous populations. The package is built upon a state-space age-structured stock assessment model - Woods Hole Assessment Model (WHAM). WHAM can incorporate time- and/or age-varying process errors (treated as random effects) on population processes, and can also incorporate environmental effects on population processes. This package provides tools to simulate and evaluate management strategies (combination of data collection methods, stock assessment models, and harvest control rules) by allowing users to explore various scenarios and their impacts on estimation and management performance for fish stocks with different life-history traits across different regions.

Capabilities of the whamMSE package:

1. Different fish life histories, spatial structures, and movement dynamics
2. Time- and/or age-varying population and fishery processes
3. Different management strategies
4. Comprehensive output analysis tools

We suggest walking through the vignettes to familiarize yourself with whamMSE: https://lichengxue.github.io/whamMSE.

# WHAM: a state-space age-structured assessment model

The Woods Hole Assessment Model (WHAM) is a general state-space age-structured stock assessment framework designed to include environmental effects on population processes. The state-space framework is attractive because it can estimate observation and process error, as well as naturally propagate random effect parameters in stock projections. WHAM can be configured to estimate a range of assessment models (see [Ex 1](https://timjmiller.github.io/wham/articles/ex1_basics.html) and [Ex 6](https://timjmiller.github.io/wham/articles/ex6_NAA.html)):

- statistical catch-at-age (SCAA) model with recruitments as fixed effects, 
- SCAA with recruitments as random effects
- "full state-space model", abundance at all ages are random effects

WHAM advances fisheries assessment because it can estimate constrained random deviations, i.e. random effects, on parameters such as:

- recruitment / numbers-at-age ([Ex 2](https://timjmiller.github.io/wham/articles/ex2_CPI_recruitment.html) and [Ex 6](https://timjmiller.github.io/wham/articles/ex6_NAA.html)),
- selectivity ([Ex 4](https://timjmiller.github.io/wham/articles/ex4_selectivity.html)),
- natural mortality ([Ex 5](https://timjmiller.github.io/wham/articles/ex5_GSI_M.html)), and
- environmental effects on the above ([Ex 2](https://timjmiller.github.io/wham/articles/ex2_CPI_recruitment.html) and [Ex 5](https://timjmiller.github.io/wham/articles/ex5_GSI_M.html))

A nice property of treating population and environmental processes as random effects is that their uncertainty is naturally propagated in projections/forecasts ([Ex 3](https://timjmiller.github.io/wham/articles/ex3_projections.html)).

Overview of WHAM presentation (Jan 8 2021): https://www.youtube.com/watch?v=o8vJvbIaOdE

Vignettes of WHAM: https://timjmiller.github.io/wham/articles.

## Installation

For the latest stable, tested release:

```
devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
```

## NOAA Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or
favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a
DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by
DOC or the United States Government.

****************************

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) | [National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [NOAA Fisheries](https://www.fisheries.noaa.gov/)


