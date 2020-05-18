# Multimodal AFM on PMMA & P3HT 
Repository containing data, analysis workflow, and I/O and analysis scripts for processing multimodal atomic force 
microscopy data in the paper titled ["Identifying Nanoscale Structure–Function Relationships Using Multimodal Atomic 
Force Microscopy, Dimensionality Reduction, and Regression Techniques."](https://doi.org/10.1021/acs.jpclett.8b01003)
If you find the notebook and workflow here useful, we'd appreciate if you'd cite our paper. 

### About the Notebook
The notebook contains extensive annotation via comments and additional information regarding the analysis can also be 
found in the SI of the paper. Broadly, the analysis steps are as follow:

1. Register PiFM and cAFM image data.
2. Apply PCA to the hyperspectral data.
3. Perform PCR using the registered image data and first 10 PCs. 
4. Calculate error metrics (F-statistic, RSS, R<sup>2</sup>)
5. Error analysis showing utility of combining PiFM and cAFM image data. 
6. Hyperspectral unmixing results using standard inputs with a variety of techniques (NMF, ATGP, NFINDR, PPI, and VCA).
7. How each SI figure was generated. 

### Other Information
The notebook was written in Python 3.5. Aside from the typical Python computing packages (like scipy, numpy, and
sklearn), [dipy](https://dipy.org/) is required. 

If you are interested in registering multimodal image data, check out [this notebook](https://github.com/kongjy/DipyImageAlignentTutorial/blob/master/Image%20Registration%20Tutorial%20with%20Dipy.ipynb) where I provide a step-by-step
tutorial (beginning with installing Python with Conda) on how to do so with dipy. It also contains commentary on what's 
happening under the hood. If you use dipy for image registration, please cite the package according to the guidelines as
indicated by the developers of the package [here](https://dipy.org/documentation/1.1.1./cite/#a-note-on-citing-our-work).
### Contact:
The corresponding author on the paper is:
```
David S. Ginger, Ph.D.
University of Washington
Alvin L. and Verla R. Kwiram Endowed Professor of Chemistry
Adjunct Professor of Physics
Chief Scientist, UW Clean Energy Institute
Founding Co-Director, Northwest Institute for Materials Physics, Chemistry, and Technology (NW IMPACT)
Washington Research Foundation Distinguished Scholar
Associate Editor, Chemical Reviews
E-mail: dginger at uw dot edu
```
Questions regarding the code, analysis, and workflow can also be directed to:
```
Jessica Kong
University of Washington
Department of Chemistry
E-mail: kongjy at uw dot edu
```
