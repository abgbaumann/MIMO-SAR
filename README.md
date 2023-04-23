# Data Processing of Multiple Input Multiple Output - Synthetic Aperture Radar (MIMO-SAR) acquistions

This repository provides the functions and scripts for processing radar acquisitions, especially the ones acquired with the Texas Instruments TIDEP-01012 device. The codes provided in this repository have been used in the following publications:

### [MIMO-SAR Interferometric Measurements for Structural Monitoring: Accuracy and Limitations](https://doi.org/10.3390/rs13214290)

[Andreas Baumann-Ouyang](https://orcid.org/0000-0002-8306-3642), [Jemil Avers Butt](https://orcid.org/0000-0002-0858-2813), [David Salido Monzú](https://orcid.org/0000-0003-4274-6874), and [Andreas Wieser](https://orcid.org/0000-0001-5804-2164)

### [Bridge deformations during train passage: monitoring multiple profiles using concurrently operating MIMO-SAR sensors](http://doi.org/10.4995/JISDM2022.2022.13620)

[Andreas Baumann-Ouyang](https://orcid.org/0000-0002-8306-3642), [Jemil Avers Butt](https://orcid.org/0000-0002-0858-2813), and [Andreas Wieser](https://orcid.org/0000-0001-5804-2164)

### [Estimating 3D displacement vectors from line-of-sight observations with application to MIMO-SAR](https://doi.org/10.1515/jag-2022-0035)

[Andreas Baumann-Ouyang](https://orcid.org/0000-0002-8306-3642), [Jemil Avers Butt](https://orcid.org/0000-0002-0858-2813), and [Andreas Wieser](https://orcid.org/0000-0001-5804-2164)

### [MIMO-SAR Interferometric Measurements for Wind Turbine Tower Deformation Monitoring](https://doi.org/10.3390/en16031518)

[Andreas Baumann-Ouyang](https://orcid.org/0000-0002-8306-3642), [Jemil Avers Butt](https://orcid.org/0000-0002-0858-2813), [Matej Varga](https://orcid.org/0000-0002-3453-169X), and [Andreas Wieser](https://orcid.org/0000-0001-5804-2164)


### Environment Setup

The code is provided in m-files, which is native to the MATLAB software. All the scripts have been tested on the MATLAB R2022a version but should also run for R2021a and higher. 

The code can be used directly after cloning this repository and opening in MATLAB. If the video tracking function is used, then is recommend to download and install davinci_draw_R2017a from http://davinci-draw.com/documentation/installation_and_quick_start/ to visualise the north arrow and some axes arrows.

### Citation

If you found this code useful, please consider citing one of the following papers:

```diff
@article{baumann2022,
  Author = {Baumann-Ouyang, Andreas and Butt, Jemil Avers and Wieser, Andreas},
  Title = {Estimating 3D displacement vectors from line-of-sight observations with application to MIMO-SAR},
  Journal  = {Journal of Applied Geodesy},
  Year = {2023},
  Publisher = {de Gruyter},
  DOI = {10.1515/jag-2022-0035}
}
```

```diff
@article{baumann2023windturbine,
  Author = {Baumann-Ouyang, Andreas and Butt, Jemil Avers and Varga, Matej and Wieser, Andreas},
  Title = {MIMO-SAR Interferometric Measurements for Wind Turbine Tower Deformation Monitoring},
  Journal = {Energies},
  Volume = {16},
  Year = {2023},
  Number = {3},
  Article-Number = {1518},
  ISSN = {1996-1073},
  DOI = {10.3390/en16031518}
}
```
