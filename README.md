# Chemical Compound Time Trend Analysis

This document provides guidelines for preparing the dataset and performing time trend analysis for chemical compounds.


## Dataset Preparation

To ensure proper analysis, the dataset should be formatted as follows:

| Name  | Condition | Concentration | Group |
|-------|----------|--------------|-------|
| Compound A | Condition 1 | 10.5 | Group 1 |
| Compound B | Condition 2 | 15.2 | Group 2 |

Each column represents:
- **Name**: The name of the chemical compound.
- **Condition**: The experimental condition under which the data was collected.
- **Concentration**: The measured concentration of the compound.
- **Group**: The experimental group classification.

## Time Trend Analysis

To perform time trend analysis, follow these steps:

1. **Set the data path**  
   Ensure the dataset is stored as an `.xlsx` file. Specify the file path in the script.

2. **Run the analysis script**  
   Source the provided R script. Upon execution, the following outputs are expected:
   - The **time constant** will be displayed in the terminal.
   - A **trend plot** for the specified compound will be generated and saved in the same directory as the dataset.

By following these steps, you will be able to analyze the temporal trends of chemical compounds efficiently.

## Publication and Citation

The time trend analysis workflow described here has been applied in our published research on the chemical transformation of aging tire and artificial turf crumb rubber. The study has been published in [*Environmental Science & Technology*](https://pubs.acs.org/journal/esthag).

If you use this analysis workflow or the associated scripts, please cite the following article:

**Recommended citation:**

McMinn, M. H.; Tang, Y.; Berger, P.; Poisson, K.; Lima, A. O. T.; Stubbins, A.; Güler, A. T.; Tian, Z. *From the Road to the Field: Decoding Chemical Transformation in Aging Tire and Artificial Turf Crumb Rubber. Environmental Science & Technology **2026**, *60* (1), 1051–1062. https://doi.org/10.1021/acs.est.5c08260

**DOI:** [10.1021/acs.est.5c08260](https://doi.org/10.1021/acs.est.5c08260)

