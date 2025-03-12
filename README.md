[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6603909.svg)](https://doi.org/10.5281/zenodo.6603909) [![DOI](https://zenodo.org/badge/480557469.svg)](https://zenodo.org/badge/latestdoi/480557469)


# RUFF_CHOI_unfolded_protein_phase_separation

This repository contains the analysis code associated with the **Unfolded Protein Phase Separation** project, led by **Kiersten M. Ruff and Yoon Hee Choi**. This manuscript has been published under the title [**"Sequence grammar underlying unfolding and phase separation of globular proteins"**](https://www.sciencedirect.com/science/article/pii/S1097276522006074?via%3Dihub).

## Prerequisites

This analysis assumes a standard installation of MATLAB 9. This analysis also assumes a standard installation of Python 3 (=> **3.8.8**). For specific package requirements, see the environment.yml file, or  create a new conda environment containing all packages by running ```conda create -f environment.yml```. In addition to the analysis contained here, some simple statistical tests were performed using [GraphPad Prism v **8.0**](https://www.graphpad.com/scientific-software/prism/).

## Raw data

The .RAW files have been deposited via the [PRIDE][1]<sup>[1]</sup> partner repository to the [ProteomeXchange Consortium][2]<sup>[2]</sup> under the dataset identifier PXD033716. For convenience, the preprocessed identification and quantitation data (hereon termed raw data) have also been uploaded as an open-access [Zenodo dataset](https://doi.org/10.5281/zenodo.6603909). 

Various public databases were also queried as cited in the accompanying manuscript.

## Workflow

Initial processing of the novel mass spectrometry spectra files was completed using either [Proteome Discoverer] or [MaxQuant][3]<sup>[3]</sup>.

Individual analyses are presented within the ```src``` folder. Where processing order is important for individual analyses, Python scripts have been numbered and MATLAB scripts have been alphabetized and should be run in order before unnumbered / unalphabetized counterparts, respectively.

1. [published_Tm_disease](https://github.com/kierstenruff/RUFF_CHOI_unfolded_protein_phase_separation/tree/master/src/published_Tm_disease) - analyses related to Figure 1

2. [optoDroplets_csat](https://github.com/kierstenruff/RUFF_CHOI_unfolded_protein_phase_separation/tree/master/src/optoDroplets_csat) - analyses related to Figures 3, 4, 5, S2, and S3

3. [atomistic_simulations](https://github.com/kierstenruff/RUFF_CHOI_unfolded_protein_phase_separation/tree/master/src/atomistic_simulations) - analyses related to Figures 4A and 7

4. [UPOD_proteomics](https://github.com/kierstenruff/RUFF_CHOI_unfolded_protein_phase_separation/tree/master/src/UPOD_proteomics) - analyses related to Figures 6 and S5

## References

[1]: https://www.ebi.ac.uk/pride/archive/

1. Perez-Riverol, Yasset, Attila Csordas, Jingwen Bai, Manuel Bernal-Llinares, Suresh Hewapathirana, Deepti J Kundu, Avinash Inuganti, et al. “The PRIDE Database and Related Tools and Resources in 2019: Improving Support for Quantification Data.” Nucleic Acids Research 47, no. D1 (January 8, 2019): D442–50. https://doi.org/10.1093/nar/gky1106.

[2]: http://proteomecentral.proteomexchange.org

2. Deutsch, Eric W., Attila Csordas, Zhi Sun, Andrew Jarnuczak, Yasset Perez-Riverol, Tobias Ternent, David S. Campbell, et al. “The ProteomeXchange Consortium in 2017: Supporting the Cultural Change in Proteomics Public Data Deposition.” Nucleic Acids Research 45, no. Database issue (January 4, 2017): D1100–1106. https://doi.org/10.1093/nar/gkw936.

[Proteome Discoverer]: https://www.thermofisher.com/au/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html

[3]: https://www.maxquant.org/

3. Tyanova, Stefka, Tikira Temu, and Juergen Cox. “The MaxQuant Computational Platform for Mass Spectrometry-Based Shotgun Proteomics.” Nature Protocols 11, no. 12 (December 2016): 2301–19. https://doi.org/10.1038/nprot.2016.136.
