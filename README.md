/* This is the main README file for the project (it is the page that is displayed when you first go to the GitHub page)

Summary and important points: 

- use a combination of perturbations and measurements to help decode the system's dynamics, this package was designed to assist those attempting to understand biological dynamics by designing, performing, and analyzing perturbation scRNA-seq experiments.


-python dependencies and short general description:

    - Pandas:
        -> this library is used for data maniputation and analysis
        -> the main data structure in pandas is the DataFrame, which is a 2D table-like data structure

    - NumPy:
        -> Numerical Python
        -> used for scientific computing and data analysis 
        -> the main data structure in numpy is an array (multi-dimensional)    data structure that is efficient for storing and manipulating multi-dimensional homogeneous data
        -> has many built-in mathematical functions for operations such as logarithms, trigonometry and exponentials 
        
    - Matplotlib:
        -> used for creating high-quality 2D plots and graphs 
        -> for example, plot() function can be used to create a line plot, while scatter() can be used to create a scatter plot 
        -> provides functions for saving plots to various file formats 
        -> integrates well with other scientific computing libraries such as NumPy and Pandas 

    - seaborn:
        -> built on top of Matplotlib, and provides a high-level interface for creating statistical graphics 
        -> this library is used for data visualization and statistical graphics 
        -> provides functions for performing statistical analysis, including regession analysis, hypothesis testing, and multivariate analysis 

    - scipy:
        -> scientific python
        -> used for scientific computing and technical computing 
        -> includes algorithms for numerical computation, optimization, intergration, linear algebra, signal processing, etc. 

    - sklearn:
        -> this library is for machine learning built on top of NumPy and Matplotlib
        -> provides tools for data mining and data analysis including classifications, regression, and clustering
        -> includes tools for data preprocessing such as scaling and normalization
        -> functionalities for supervised and unsupervised machine learning
        -> functionalities for model selection and evaluation
        
    - statsmodels:
        -> provides classes and functions for statistical modeling and analysis 
        -> functionalities for regression analysis (linear, generalized linear, mixed-effects models, and time series regressions)

    - Facebook PCA (though sklearn SVD or PCA should also work):
        -> Facebook PCA is a technique to compress large datasets and reduce dimensionality of data while preserving its underlying structure 
        -> various python libraries can be used to perform PCA such as NumPy

    - networkx:
        -> used for creating and manipulating various types of graphs including directed and undirected graphs, and weighted and unweighted graphs
        -> you can search for nodes and edges, calculate various graph metrics, and apply various graph algorithms 
        -> provides tools for visualizing graphs 

    - Infomap (this can be removed):
        -> this package is a python implementation of the Infomap algorithm, which is a community detection algorthim for complex networks 
        -> can be used to analyze the community structure of various types of networks

    - goatools 
        -> tool set for working with Gene Ontology (GO) annotations and performing gene ontology analyses 
 


- the Power_Analysis_DOE folder in this repository contains an iPython notebook that shows a rough comparison of the pilot scRNA-seq to population RNA-seq of the same perturbation
- the original experiment used a high multiplicity of infection, meaning a lot of sgRNAs compared to cells to ensure a large fraction of the cells get delivered a perturbation

Computational Framework:
    Inputs:
        -> An expression matrix 
        -> Guide barcode (GBC) PCR data to pair perturbations with cell barcodes 
        -> A database of preassociated sgRNA/GBC pairs 
    
    Intermediate calculations: 
        -> Guide barcodes and cell barcodes have to be paired accurately, accounting for chimeric PCR products (when the PCR amplifies multiple templates simultaneoulsy, resulting in a PCR product that can contain a mixture of sequences from the different templates)
        -> A simple fitness calculation is possible by determining the difference between the initial abundances of a GBC and how many cells it appeared in

        -> MOI and detection probability are evaluated by comparing the observed number of GBCs/cell to a poission distribution that is first zero truncated, due to FACS selection for transduced cells, and then zero inflated, due to detection dropout (an iPython notebook showing an example of this can be found under the Cell_States folder) 
        -> A Cell state classifier (computational tool that uses machine learning algorithms to predict the state or identity of a single cell based on its gene expression profile) is defined on wildtype or control cells and then applied to all cells in an experiment (an iPython notebook showing an example of this can be found under the Cell_States folder). These classifications can used as outputs to be predicted (instead of gene expression) or as covariates in the model.
        -> The linear model integrating all covariates (and interactions terms as desired) is fit. An EM-like approach filters cells that look much more like control cells than perturbed cells (an iPython notebook showing an example of this can be found in the contrived-em_example.ipynb) 

    Outputs: 
    -> regulatory coefficients
    -> cell state effects are obtained by predicting the cell states based on the linear model instead of predicting gene expression
    -> cell size effects (genes detected or transcripts detected) can be predicted as well

*/

##########################################################################################################################################


<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/mimosca_logo.png" title="MIMOSCA" alt="MIMOSCA" height=99 width=372>

# Multiple Input Multiple Output Single Cell Analysis

Named after <a href="https://en.wikibooks.org/wiki/Control_Systems/MIMO_Systems">MIMO control systems</a>, and the associated methods in Systems Identification to use a combination of perturbations and measurements to help decode the system's dynamics, this package was designed to assist those attempting to understand biological dynamics by designing, performing, and analyzing perturbation scRNA-seq experiments.

Python Dependencies: pandas, numpy, matplotlib, seaborn, scipy, sklearn, statsmodels, <a href="https://github.com/facebook/fbpca">Facebook PCA</a> (though sklearn SVD or PCA should also work), networkx , and <a href="http://www.mapequation.org/code.html">Infomap</a> (this can be removed), and goatools 

## Related Resources
* <a href="http://www.sciencedirect.com/science/article/pii/S0092867416316105">Our paper</a>
* <a href="https://groups.google.com/forum/#!forum/perturb-seq">Perturb-seq Google Forum</a>
* <a href="http://biorxiv.org/content/early/2016/12/12/093237">Chimera Correction</a>  and <a href="https://github.com/asncd/schimera">Code</a>
* <a href="http://www.clontech.com/US/Products/Genome_Editing/CRISPR_Cas9/Resources/Online_tools_for_guide_RNA_design">sgRNA Design Tools</a>
* Addgene Plasmids:  <a href="https://www.addgene.org/85801/">pPS</a>, <a href="https://www.addgene.org/85967/">pBA439</a>
* <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90063">GEO Database link</a>



## Contents


* [FAQ](https://github.com/asncd/MIMOSCA/blob/master/README.md#faq)
* [Design of Experiments](https://github.com/asncd/MIMOSCA/blob/master/README.md#design-of-experiments--power-calculations)
* [Guide barcode, Cell barcode pairing](https://github.com/asncd/MIMOSCA/blob/master/README.md#guide-barcode-cell-barcode-pairing)
* [Computational Workflow](https://github.com/asncd/MIMOSCA/blob/master/README.md#computational-workflow)

Please let me know (here or in the Google Forum) if there are any areas that you'd like to see improved or new items to be added! 

## FAQ

* Q: How many perturbations can I do?
* A: It costs between $0.1/cell and $0.2/cell for commercial droplet scRNA-seq methods, it takes ~10 cells/perturbation to observe signature effects and ~100 cells/perturbation to see individual gene level effects robustly. Based on your budget, you can crunch the numbers (also see Design of Experiments section below)


## Design of Experiments & Power Calculations

<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/comp_knob.png" title="Experimental Design" alt="Experimental Design" height=360 width=432>

For a rough comparison of our pilot scRNA-seq to population RNA-seq of the same perturbation, see this <a href="https://github.com/asncd/MIMOSCA/blob/master/Power_Analysis_DOE/ost_ko_comparison.ipynb">iPython notebook</a>.

In designing Perturb-seq like experiments, there are a few key factors to keep in mind:

### Signatures vs. individual transcript-level phenotypes
Are you interested in broad transcriptional signatures or individual gene level differential expression? If the former, a rough approximation may be around 10 cells/perturbation. If the later, 100 or more cells may be required based on the effect size. 

A similar approximation for reads/cell would be a couple thousand for signatures and tens of thousands for gene-level.

### Library Size and Representation

As in any pooled screen, the representation of each perturbation in the library will vary. With genome wide CRISPR libraries the difference between the 10th and 90th percentile of a library is roughly 6-fold (Wang, 2013). Depending on how much a user wants to ensure every member of the library is represented, the cells/perturbation factor should be multiplied by an additional factor to reflect this variance.

### Using High MOI to infer genetic interactions

Our approach to use high MOI instead of either a single vector with multiple sgRNAs or vectors with different selection methods benefits from ease of implementation and the ability to represent a large diversity of combinations (only limited by the number of cells). 

However, challenges include a Poisson-like variance in the number of sgRNA/cell, sgRNA detection sensitivity, and the formation of PCR chimeras during the enrichment PCR procedure that can create misassignments. 

All three of these factors should be assessed in pilot experiments to troubleshoot. An example of such a pilot is shown below (modified from the Drop-seq style species mixing experiments): 

<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/species_mix.png" title="Species Mix" alt="SMIX" height=431 width=365>

## Guide Barcode, Cell Barcode Pairing

<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/perturbseq_readdistribution.png" title="pseq_plasmid" alt="pseq_plasmid">

The distribution of reads going to the Perturb-seq vector (antiparallel) from 10X RNA-seq is shown above. Note that while the expression of the construct is comprable to that of a housekeeping gene, only a fraction of the reads overlap with the 18bp barcode (colored section in the coverage track). As such, it is advisable in most cases when you have a short barcode to perform enrichment PCR to obtain sensitive GBC/CBC pairing.

## Computational Workflow

<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/comp_flow2.png" title="Overview" alt="Overview" height=662 width=569>

### Inputs
* An expression matrix output by a high throughput scRNA-seq protocol (such as <a href="http://mccarrolllab.com/dropseq/">Drop-seq</a> or <a href="https://support.10xgenomics.com/single-cell/software/pipelines/latest/what-is-cell-ranger">10X cellranger</a>)
* Guide barcode (GBC) PCR data to pair perturbations with cell barcodes (for certain applications this may be able to be directly obtained from the RNA-seq data
* A database of preassociated sgRNA/GBC pairs (either by Sanger sequencing or NGS)

### Intermediate Computation
* Guide barcodes and cell barcodes have to be paired accurately, accounting for chimeric PCR products
* A simple fitness calculation is possible by determining the difference between the initial abundances of a GBC and how many cells it appeared in
* MOI and detection probability are evaluated by comparing the observed number of GBCs/cell to a poission distribution that is first zero truncated, due to FACS selection for transduced cells, and then zero inflated, due to detection dropout (see this <a href="https://github.com/asncd/MIMOSCA/blob/master/Cell_States/DC_celltypes-wt-github.ipynb">iPython notebook</a>)
* A Cell state classifier is defined on wildtype or control cells and then applied to all cells in an experiment (see this <a href="https://github.com/asncd/MIMOSCA/blob/master/Cell_States/DC_celltypes-wt-github.ipynb">iPython notebook</a>). These classifications can used as outputs to be predicted (instead of gene expression) or as covariates in the model.
* The linear model integrating all covariates (and interactions terms as desired) is fit. An EM-like approach filters cells that look much more like control cells than perturbed cells (see this <a href="https://github.com/asncd/MIMOSCA/blob/master/contrived-em_example.ipynb">iPython notebook</a> for an example)

### Outputs

* The regulatory coefficient obtained from the model are the most informative output giving an estimate of what extent each covariate (perturbation, cell state, pairwise interaction between perturbations, etc) impacted a given gene.
* Cell state effects are obtained by predicting the cell states based on the linear model instead of predicting gene expression
* Cell size effects (genes detected or transcripts detected) can be predicted as well



