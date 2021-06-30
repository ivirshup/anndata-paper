---
title: 'anndata: Annotated data'
authors:
  - name: Isaac Virshup
    orcid: 0000-0002-1710-8945
    affiliation: "1,†"    
  - name: Sergei Rybakov
    affiliation: "2"
  - name: Philipp Angerer
    affiliation: "2,3,†"    
  - name: F. Alexander Wolf
    orcid: 0000-0002-8760-7838
    affiliation: "2,3,†"
affiliations:
 - name: University of Melbourne.
   index: 1
 - name: Helmholtz Munich, Institute of Computational Biology, Munich, Germany.
   index: 2
 - name: Cellarity, Cambridge, MA.
   index: 3
 - name: corresponding authors
   index: †
date: May 1st, 2021
bibliography: paper.bib
---

# Abstract

anndata is a python software package for handling annotated datasets both in memory and on disk. It focuses on enabling intuitive iterative data science workflows.


# Statement of need

In exploratory data analysis, generating insight from high-dimensional data is typically achieved through learning patterns that allow (i) to condense data into meaningful lower-dimensional representations and (ii) to assign semantic meaning to observations and variables.
Learning these patterns almost always involves workflows of iteratively training models on pre- and post-learned annotations of data, requiring to book-keep their representations and scalar annotations, such as labels and numerical scores.
anndata's purpose is to make such workflows as efficient as possible through a data structure that naturally integrates book-keeping with model training and analysis, well-integrated into the pydata ecosystem.

A particularly relevant use case with high degrees of iterations and many annotations involved concerns computational biology, where 
advances in single-cell high throughput sequencing (scRNA-seq) have given rise to new classes of analysis problems.
While previous bulk studies had few observations with known labels, current datasets have many observations with little sample level information, in high dimensions, and with a high degree of sparsity. Neither xarray nor pandas meet these needs, anndata offers the sparse data support and efficiency needed.

While analysis of bioinformatic data has been dominated by the R ecosystem, the recent explosion in popularity and availability of machine learnings tools in the Python ecosystem is important to take advantage of. anndata offers an in-memory representation that seamlessly integrates with the Python ecosystem, while offering a cross-ecosystem on-disk format that allows interfacing with the R ecosystem.


# Defining the AnnData object

AnnData was inspired by similar data structures within the R ecosystem, in particular, ExpressionSet, and single-cell related more recent alternatives, like SingleCellExperiment and the Seurat on-disk format.

Within the pydata ecosystem the closest package that would be amenable to serve this paradigm is xarray [@Hoyer2017], which enables to deal with highly complex labelled data tensors of arbitrary dimensions - keeping labels on data is useful [@Hoyer2017].
By contrast, the highly popular package pandas [@McKinney2010] operates only on `DataFrames`, that is, single tables of data.
anndata is positioned in between by providing the minimal additional structure to enable storing compact annotations and representations of high-dimensional data, making the book keeping during learning from it much simpler. Current learning practices in data analysis libraries such as scikit-learn [@Buitinck2013](Sec. 2.2) model input and output for each computation as set of arrays. To organize process, AnnData first defines a particular data semantics for it.


## AnnData's data semantics is designed for an exploratory data analysis workflow

AnnData models a dataset as a set of observations and variables with measured, annotated, and derived values \autoref{fig:overview}.
Observations (samples) and variables (features) here take much the same meaning as they do in the tidy data [@Wickham2014] framework.
Measured values are recorded on the cross product of observations and variables, e.g. on the `X` or in the `layers` attributes.
In this way, the observations and variables can be thought of as discretely valued principal dimensions of the dataset.
Each value on these dimensions is given a label, stored in `obs_names` and `var_names` respectivley.

Annotations and derived values can then be stored on the dimension specific axes.
Simple annotations and derived values which can be stored in a single vector are added to the main annotation dataframes for each axis, `obs` and `var`.
Learned representations are added to `obsm` and low-dimensional manifold structure to `obsp`. 
Annotations added here include values like alternative names (e.g. different identifier mappings or categories for each variable).
Derived values added here can be descriptive statistics (e.g. mean and variance) or categories from clustering.

In contrast to the R ecosystem, observations be in the rows and variables be in the columns is the convention of the modern classics of statistics [@Hastie2009] and machine learning [@Murphy2012], tidy data [@Wickham2014], the convention of dataframes both in R and Python, and the established statistics and machine learning packages in Python (statsmodels, scikit-learn).

![**Structure of the AnnData object.**
*(a)* The AnnData object is a collection of arrays aligned to the common dimensions of observations (`obs`) and variables (`var`).
This was designed to organize analysis results, fitting in with the common conventions of statistical/ machine learning.
Here, color is used to denote elements of the object, with "warm" colors selected for elements aligned to the observations and "cool" colors for elements aligned to variables.
The object is centered around the main data matrix `X`, whose two dimensions correspond to observations and variables respectivley.
Primary labels for each of these dimensions are stored as `obs_names` and `var_names`.
`layers` elements of the same shape as `X` to allow for multiple representations (e.g. different normalization strategies).
Simple annotations (e.g. 1d vectors of labels or statistics) for each dimension are stored in dataframes `obs` and `var`.
`obsm`, `varm` contain multidimensional arrays whose first dimension are aligned to their respective dimension.
Pairwise relationships within each dimension can be stored in `obsp` and `varp`.
Data which doesn’t fit this model, but should stay associated to the dataset can be stored in `uns`.
As examples, *(b)* the response variable ŷ learned from the data is stored as an annotation of it’s observations.
*(c)* Reduced dimensional representations PCA are stored with observation/ variable loadings aligned to the main dimensions.
*(d)* A K nearest neighbor representation of this PCA space is represented as an adjacency matrix, constituting a pairwise relationship of the observations, fitting in `obsp`.
*(e)* Subsetting the `AnnData` object by observations produces a view subsetting all elements aligned to this dimension.
\label{fig:overview}
](figures/overview.pdf)


## Data analysis workflow - iteratively learning representations and scalar annotations

Let us discuss a few canonical examples of the exploratory data analysis workflow. Fitting a classification, regression, or clustering to high dimensional data gives rise to a response variable ŷ learned from the data is stored as an annotation of its observations \autoref{fig:overview}[b].
Reduced dimensional representations PCA are stored with observation/ variable loadings aligned to the main dimensions \autoref{fig:overview}(c). A K nearest neighbor representation of this PCA space is represented as an adjacency matrix, constituting a pairwise relationship of the observations, fitting in `obsp` \autoref{fig:overview}(d). Subsetting the `AnnData` object by observations produces a view subsetting all elements aligned to this dimension. \autoref{fig:overview}(e). 


## Efficient data operations for efficient data analysis workflows

Due to the ever increasing scale of data AnnData is working with, emphasis has been placed on providing efficient data handling operations. There has been an overall emphasis on having low memory and runtime overhead. This is accomplished in a number of ways. To this end, AnnData offers sparse data support, out of core conversions between dense and sparse data, lazy subsetting, per element operations for low total memory usage, in place subsetting, combining AnnData objects with various merge strategies, and a backed out-of-memory mode.

Deep support for sparse data. `AnnData` takes great pains to support efficient operations with sparse data. While there currently is no equivalent API for working with sparse and dense data in the python ecosystem, `AnnData` abstracts over this making it much easier for novices to handle each. As mentioned above, on-disk formats for sparse data have also been defined, along with operations for out of core access to this data.

Subsetting anndata objects is lazy. This takes advantage of the fact that a great deal of the exploratory data analysis process is read-only, and that data is often sliced just for access to a subset of one element. For typical use cases of tidy-data (and for data frames), data storage is columnar (or "variable major"). Our access patterns to X are typically row based, so we use CSR and C order arrays (or "observation major"), which allows efficiently accessing batches of the dataset, to meet the needs of batched learning algorithms.

Datasets can be joined along variables or observations.
That is, from multiple individual dataset can be combined to have a superset of either the observations or variables, depending on the direction of concatenation.
The other dimensions are merged to contain either the union or intersection of labels.


## An on disk representation for sharing data analysis results

An AnnData object captures a useful unit (the dataset) in the data analysis workflow.
Providing a stable, and standard on disk format for this unit relieves the pain of working with many competing formats for each individual element.

Another big advantage is the on-disk representation, which even for pandas DataFrames is not yet resolved in a canonical way.
For instance, there is none of the binary persistent formats are able to represent all entry types of AnnData.
For instance, even such a key data type a categorical data types are not yet represented in the HDF5 format.
Pickled dataframes are stable, but they are non-persistent.

In the R ecosystem, in-memory objects are serialized and written to disk.
This is problematic since that data cannot be read by another tool, and may become inaccessible even after software updates.
If one chooses to use standard formats to represent all elements of a dataset, a set of standards has to be chosen.
AnnData has chosen self-describing hiearchical data formats such as HDF5 and `zarr` [https://doi.org/10.5281/zenodo.3773449] for this purpose.
AnnData objects can be efficiently saved to disk using standardized formats \autoref{fig:ecosystem}.
This means the data is accessible from other programming environments, as opposed to a serialized format like `pickle` or `Rdata`.

By choosing standardized formats, stored data can be accessed from a variety of ecosystems including `python`, `julia`, `R`, `java`, and `javascript`.
While the project has tried to stick to standardized formats, there are a few cases where no standards existed within our models.
An especially important example of this is sparse array formats, which are critical for efficient processing of scRNA-seq data.
To account for this, we define schemas for these types, which specify how these elements can be read from disk to memory.
These specifications are versioned and stored in an internal registry.
Versioning allows the specifications to evolve with the project while maintaining the ability to access older data.

Like the AnnData object itself, the on-disk representations of these objects closely mirrors their in-memory representation.
Compressed sparse matrices (CSR and CSC format) are stored as a collection of three arrays, `data`, `indices`, and `indptr`, while tabular data is stored in a columnar format.


# How the ecosystem uses the AnnData object

Transcriptional data is stored in a large variety of formats. The distributed nature of research can lead to fractured ecosystems without consortia for organization. 

## AnnData provides conventions for data handling that are used by many tools

AnnData provides a common format and set of conventions for handling numeric datasets (like those generated in scRNA-seq). This consists of an in memory model, which operates as the core data model for a number of tools \autoref{fig:ecosystem}. The data can be moved back and forth to disk. Since it is stored in standardized formats, the dataset is distributable a wider ecosystem of tools. These include data portals, viewers, and the ecosystem beyond python.

![**AnnData provides common conventions for data handling for a variety of tools.**
*(a)* Data flows using the `anndata` model. `AnnData` objects can be created from a number of formats, including common delimited text files, or domain/ tool specific formats like `loom` or `cellranger` outputs.
Once in memory, AnnData provides an api for handling annotated matrix objects, proving a common base object used by a range of analysis tools.
The in memory format has a one to one relationship with it's on disk format.
The on disk format for this model uses lanugage independent technologies, facilitating use by other tools and interchange with other ecosystems.
*(b)* The on disk schema for maps the schema to a hierarchical model (mapping of elements indicated by color).
Each element is annotated with a type and schema version to facilitate interchange.
*(c)* AnnData is widely used in the single cell RNA seq ecosystem.
\label{fig:ecosystem}
](figures/ecosystem.pdf)

## Examples of use for analysis of spatial transcriptomics, RNA velocity, and multiple modalities

AnnData is widely used in single cell analysis across spatial transcriptomics, RNA velocity, and multiple modalities.

![
**Examples across spatial transcriptomics, RNA velocity, and multiple modalities.**
\label{fig:examples}
](figures/examples.pdf)


# Future directions

The AnnData project is under active development and will have more features. These include, but are not limited to, more advanced out of core access, a split-apply-combine framework, integration with more of the python ecosystem, and interchange with more formats like apache Arrow.

Non-homogeneous data in `X`, to enable learning from Electronic Health Records.

# Author contributions

Isaac has led the anndata project since v0.7, and contributed as a developer before. He rewrote wide parts of the code introducing high efficiency, robustness, concatenation, ... Sergei made diverse contributions to the code base, in particular, ... Phil suggested to replace Scanpy's initial unstructured annotated data object to one mimicking R's ExpressionSet, and wrote AnnData's [first implementation](https://github.com/theislab/scanpy/commit/315859c5586116434ea3b7ce97512a5e2a1030e2) with indexing, slicing, ... and ascertained good software practices in the project. ... Alex led the project until 0.7, introduced the idea of centering data science workflows around an [initially unstructured annotated data object](https://github.com/theislab/scanpy/tree/c22e48abe45a6ccca5918bbf689637caa4b31250), guided the user design of the first implementation, implemented most of the subsequent early functionality and wrote the documentation until 0.7: reading & writing, paired in-memory manipulation and on-disk format, sparse data support on disk and fast loading, backed mode, views, AnnData's slots beyond `.obs` and `.var`, and numerous smaller design choices to integrate AnnData well into Scanpy's workflows.

# References
