## Hypothesis testing for partial tail correlation in multivarate extremes

### 1. Script and its source files for simulation studies.

* Script_simulation.R: Main script outlineing the steps to reproduce the results in the manuscript.
* powertest.R: Functions for conducitng power tests.
* TransformedOperations.R: Functions that define transformed linear operations.
* functions.R: Functions that estimate parameters.

### 2. Script and its source files for applications.

* Script_Danube.R: Script for performing hypothesis tests and generating graphs for the Danube application.
* functions_graph.R: Functions for graph construction and visualization.

### Rdata files for simulation study:

* Output_unif.Rdata: Simulation outputs where the matrix C is drawn from a uniform distribution.
* Output_ar1.RData: Simulation outputs where the matrix C is created from an AR(1) stucture.

### Rdata files for applications:

* danube.rda: Danube dataset obtained from the graphicalExtremes pacakges.
* Data.Zip: Datasets obtained from Asadi.et al (2015), provided in the Data.zip file within the Data folder.
