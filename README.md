This is a full end to end analysis and prediction pipeline in the form of a Shiny App built off of the standalone, BreedStream R package (https://github.com/JohnSearl007/BreedStream).

This pipeline:
1) Allows the user to visually clean raw phenotypic files for erronous plot values.
2) Formats the cleaned phenotypic data in a standardized mannor. Removes border/filler plots. Calculates yield from plot combine weight and moisture readings (adjustable for different crops and plot sizes).
3) Generates diagnostic plots for quality control measures on the first stage model before the user proceeds further in analysis.
4) Allows the user to generate in silico genotypes.
5) Fits a Two-Stage model with many options including allowing the user to utilize an optimized H Matrix based on AIC values.
6) Supports the usage of mutli-trait and restircted interval selection indices.
7) Estimates additive and optionally dominant marker effects.
8) Estimates a usefulness criterion for use in the development of next cycle breeding crosses.
9) Enables deployment of Optimal Mate Allocation for the selection of next cycle breeding crosses.
10) Performs RR-BLUP predictions on tested and untested hyrbids.
