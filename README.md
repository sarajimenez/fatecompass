# FateCompass

FateCompass is a flexible pipeline to estimate dynamic **stochastic trajectories** of gene expression and **transcription factor activities** during **cell-fate decision** using **single-cell RNA seq data**. We combined a probabilistic framework with RNA velocities ([La Manno et al. 2018](https://doi.org/10.1038/s41586-018-0414-6) and [Bergen et al. 2020](https://doi.org/10.1038/s41586-018-0414-6)) or a differentiation potential to estimate transition probabilities and perform stochastic simulations. Also, we implemented a linear model of gene regulation ([Balwierz et al.2014](http://www.genome.org/cgi/doi/10.1101/gr.169508.113)) to learn transcription factor activities. Taking into account dynamic changes and correlations, we identified lineage-specific regulators. 

![](images/fatecompass.png)


