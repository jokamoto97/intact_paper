This directory contains all data and code necessary to reproduce the analysis from Okamoto et al. 2022 https://doi.org/10.1101/2022.07.19.500651.
To run all analyses, first, be sure that you have the following software packages installed:

bedtools (v2.22.0),
R packages:
    stringr (1.4.0),
    fgsea (1.12.0),
    tidyr (1.1.3),
    ggplot2 (3.3.6),
    dplyr (1.0.8),
    data.table (1.14.0),
    qvalue (2.18.0),
    SQUAREM (>= 2021.1),
    bdsmatrix (>= 1.3-4),
    numDeriv (>= 2016.8-1.1)

Then, execute

```snakemake --cores 8 all```


The above command reproduces each figure, table, and supplemental file from the manuscript. Results are generated in the output/ directory.






If you are interested in a certain part of our work, you can save time by only running parts of the analysis:


# Simulation analysis:
*We do not include the code to generate the simulated data since it relies on individual-level GTEx genotypes.

```snakemake simulation_study```



# Hukku et al. data analysis (everything but the GO term enrichment analysis)

```snakemake real_data_hukku```



# Hukku et al. data analysis (GO term enrichment analysis)
*This step may take several hours to run.

```snakemake go_gse```

```snakemake summarize_go_gse```



# Sinnott-Armstrong data analysis

```snakemake real_data_sinnott_armstrong_bedtools_intersect```

```snakemake real_data_sinnott_armstrong```


# Support


Please contact xwen@umich.edu or jokamoto@umich.edu if you have any questions.
