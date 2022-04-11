# Public repository template for scientific manuscripts

This repository contains a collection of templates designed to ease the creation of public-facing code and data repositories for scientific manuscripts. The templates provided are geared toward python-directed computational analyses of biological datasets, but could be adapted to suit other purposes at the users discretion.

## Getting started
---

### 1. Clone the template repository

To use this repository, first clone the repository using the ```Use this template``` button, and check the box to include all branches. It is recommended to title the repository using the *AUTHOR_Running-title syntax* for ease of use and access. Once cloned, make sure to locally checkout the branch of interest before returning to the master branch to merge: 

### 2. Select your template branch

Once cloned, find your template of interest (currently provided are *general*, *image-analysis* and *mass-spectrometry* examples) by inspecting each of the relevant branches. Once you have a local copy of the branch of interest, return to the ```master``` branch and overwrite it with your chosen template branch: ```git merge -Xtheirs <BranchName> --allow-unrelated-histories```. Optionally, you may also remove the remaining template branches using ```git branch -d <BranchName>``` (delete local branch) and ```git push origin --delete <BranchName>``` (delete the remote version) for each of the unused branches.

### 3. Add components

Add your codebase under the ```src``` folder and example data to the ```data``` folder. An optional ```utilities``` folder is also provided to house scripts that are accessed via relative imports in the ```src``` files. 

Add your environment file (an [environment.yml](environment.yml) example is provided), including relevant version constraints. For example, an environment file containing only the components explicitly installed in a conda environment can be generated using the command ```conda env export --from-history```.

Check (or replace) the [license](LICENSE) file, and ensure that it provides appropriate permissions for anyone wishing to repurpose your codebase/dataset. 

Finally, check (or replace) the [.gitignore](.gitignore) file. A standard python version is provided.

### 4. Edit the README

In the newly-updated main branch (after merging your template branch of interest), the [README](README.md) will update to the relevant front matter. Edit this file to reflect your manuscript details, and include a summary of the workflow(s) to be contained within the repository.

Some commonly-used references are provided in each template, as well as an example table for providing detailed workflow descriptions. Alternatively, users familiar with python may wish to try a specific workflow management method, such as [SnakeMake](https://snakemake.readthedocs.io/en/stable/) or [YAWL](https://yawlfoundation.github.io/).

### 5. *RECOMMENDED*: Zenodo repository for code archive and DOI creation

In addition to providing publically-accessible code via this GitHub repository, services such as [FigShare](https://figshare.com/) and [Zenodo](https://zenodo.org/) can connect with GitHub to generate archived and versioned copies of individual repositories. Using this integration, it is then possible to generate a DOI which can be included in your manuscript to point to the accompanying codebase here.

 It is recommended to do this once editing your cloned version of this repository is complete, and space is provided at the top of the README in each template to house the associated repository badges.


### 6. *OPTIONAL*: External repositories for storing raw data collections

Many disciplines maintain technique-specific repositories for all raw data associated with published manuscripts, for example the PRIDE repository housing raw mass spectrometry data. It is recommended to use these repositories in the first instance, and usage instructions given by the provider should be followed. In most cases, these databases will generate a unique identifier or DOI which can then be linked in the accompanying repository here. For example, space is provided for PRIDE identifiers for mass spectrometry data etc.

In the case of no specific repositories being available, general-use examples include [FigShare](https://figshare.com/) and [Zenodo](https://zenodo.org/). Both accept raw dataset submissions, generate a DOI and provide long-term public managed access to your dataset.

Once any external data respository DOI's have been generated, add the associated repository badges to the top of the README.

## Disclaimer
---

*This template repository was designed for personal use and is provided as-is. Whilst I endeavour to keep it up-to-date and respond to issues raised here, I can provide no guarantee of the completeness, accuracy, reliability, suitability or availability of the information, services and software contained here for your use case.*