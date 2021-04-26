---
title: 'anndata: annotated data'
authors:
  - name: Isaac Virshup
    orcid: 0000-0002-1710-8945
    affiliation: "1"    
  - name: Sergei Rybakov
    affiliation: "2"
  - name: Philipp Angerer
    affiliation: "2, 3"    
  - name: F. Alexander Wolf
    orcid: 0000-0002-8760-7838
    affiliation: "3, 2"
affiliations:
 - name: University of Melbourne.
   index: 1
 - name: Helmholtz Munich, Institute of Computational Biology, Munich, Germany.
   index: 2
 - name: Cellarity, Cambridge, MA.
   index: 3
date: 1 October 2020
bibliography: paper.bib
---

# Abstract

anndata is a python software package for handling annotated datasets.

# Statement of need

<!--
Broadly define problem
Specifics of how we handle it, what do we uniquely solve

* Structured representation of a machine learning dataset
* scRNA-seq, looks like these datasets
  * Specific needs for sparse data
  * Interchange between ecosystems, on disk format
* Previous/ parallel solutions (xarray, netCDF / CF conventions) while close do not solve all these problems

-->

<!--
Be more specific about iterative nature of exploratory analysis, especially when little is known about the data a-priori. Broader, less formal.

We want to capture most of this process.
-->
Generating insight from high-dimensional data often involves two steps: condensing them into useful representations and assigning semantic labels.
Both is achieved by learning patterns in high dimensions based on existing annotations.
Once insight is gained, knowledge is updated and assigned in form of learned annotations.
This defines core workflow within exploratory data analysis (/ data science?), which, which anndata, we strive represent in a data structure that integrates well with the pydata ecosystem.
<!-- This defines the general machine learning and data science paradigm, which, which anndata, we strive represent in a data structure that integrates well with the pydata ecosystem. -->

<!-- Put single cell motivation above? -->
The recent advances in single cell high throughput sequencing have brought new classes of analysis problems to the field.
While previous bulk studies had few samples with known labels, current datasets have many samples with little sample level information.
This has introduced a number of challenges, handling this large data and developing new methods to analyse it.
The size of the data has been a problem as it's very high dimensional (>20k genes in standard human genome annotation) and ever increasing sizes of datasets.
As a small fraction of genes are detected in a single cell, the data is highly sparse. <!-- This could go later, the main point is that tools like `xarray` and `pandas` don't handle sparse data, but it's needed -->

On the analytical tools side, the characteristics of these new large datasets are a good fit for modern machine learning methods, which are increasingly implemented in the python ecosystem.
Machine learning/ data analysis methods in python work well with data "shaped-like" data from scRNA-seq (e.g. methods in scikit-learn).
By making it easy to handle this data in python, we can more easily take advantage of these libraries.
However, siloing data in a single ecosystem benefits no one, and there is much to gain from making the data accesible from any ecosystem.

These problems are addressed by anndata, an object which provides a powerful model for representing numeric datasets and associated models, has efficient storage and operations on sparse data, and provides an accesable on disk format. <!-- There's probably a better term than accesbile -->

<!-- While these libraries have great computational tools, they frequently work with unlabelled numpy arrays/ scipy sparse arrays. -->

# Introduction/ background

Specifically, the central `AnnData` class stores observations (samples) of variables (features) in the rows of a matrix.
This is the convention of the modern classics of statistics [@Hastie2009] and machine learning [@Murphy2012], the convention of dataframes both in R and Python and the established statistics and machine learning packages in Python (statsmodels, scikit-learn).

<!-- Is this more "statement of need"? -->
NGS datasets end up being represented as a matrix of values.
These are scalar values of the probed variables for each observation in the dataset.
Recent advances in single cell methods mean that the number of observations in each study has exploded, along with a greater sparsity of values for each cell.

<!-- Move to body? -->
Within the pydata ecosystem the closest package that would be amenable to serve this paradigm is xarray [@Hoyer2017], which enables to deal with highly complex labelled data tensors of arbitrary dimensions.
On the other hand, there is the highly popular package pandas [@McKinney2010], which merely provides and operates on `DataFrames`, that is, single tables of data.
anndata is positioned in between by providing the minimal additional structure to enable storing compact annotations and representations of high-dimensional data, making the book keeping during learning from it much simpler.

With that, anndata perfectly integrates into scikit-learn [@Buitinck2013], statsmodels `[@statsmodels_author:YYYY]`, seaborn [@Waskom2021], and easily interfaces with pytorch and tensorflow.

AnnData is a structured representation of high dimensional datasets.
It was designed to efficiently represent large datasets in a user friendly way, and one that builds on existing standards rather than inventing it's own.

<!-- Why do this in python/ what's the difference from SingleCellExperiment. -->

Since the last advance in high throughput molecular profiling (micro-array -> rna-seq), Python has emerged as an extremely popular language for data analysis and machine learning.
While the analysis of bio-molecular data has previously been focussed on in the R language.

An AnnData object captures a useful unit (the dataset) in the data analysis workflow.
Providing a stable, and standard on disk format for this unit relieves the pain of working with many competing formats for each individual element.

### Figure:

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

> Figure one will be more of a "schematic". Basically the idea of `obs x vars`, how this is commonly used in the literature, and a layout of the object. This could also include what kinds of representations can be stored. E.g. data matrix, annotation, graph, and unstructured.


### Figure: Ecosystem/ usage

![**AnnData provides common conventions for data handling.**
*(a)* Data flows using the `anndata` model. `AnnData` objects can be created from a number of formats, including common delimited text files, or domain/ tool specific formats like `loom` or `cellranger` outputs.
Once in memory, AnnData provides an api for handling annotated matrix objects, proving a common base object used by a range of analysis tools.
The in memory format has a one to one relationship with it's on disk format.
The on disk format for this model uses lanugage independent technologies, facilitating use by other tools and interchange with other ecosystems.
*(b)* The on disk schema for maps the schema to a hierarchical model (mapping of elements indicated by color).
Each element is annotated with a type and schema version to facilitate interchange.
*(c)* AnnData is widely used in the single cell RNA seq ecosystem.
\label{fig:ecosystem}
](figures/ecosystem.pdf)

> A figure looking at usage of the object. Represents things like:
> * What projects use it (usage stats)
> * What problem does it solve in the ecosystem
>   * This would include at least idea of on-disk representation
> * Benchmarks?


<!-- 
I think a big part of the value proposition of AnnData is that the representation works well with the kinds of operations we want to do with single cell data. 
It fits the semantics of the problem well. How do I describe these semantics.
-->

# What is anndata?

## General features

* Representation of data
    * Labelled arrays
    * Multiple representations (layers/ obsm)
    * Associated metadata/ computed values
        * Kinds of representations
* Manipulation of datasets
    * Optimized for work with sparse, high dimensional data
    * Efficient subsetting
    * Combination
    * Out of core access
* Persistent representation
    * One-to-one on disk representation (for HDF5 and Zarr)
    * Uses standard on disk formats, does not invent it's own
    * Out of core access
    * Defined disk represenation, that will be consistent, and a way to make sure the data will be read in the future
        * This doesn't happen when you've got pickled or RDA formats
* Efficiency is key
    * Lengths have been taken to avoid copying uneccesary data
    * Lazyness in views
    * Sparse data

### Labels

Keeping labels on data is useful [@Hoyer2017].
Keeping those labels associated with the data as it moves through an analysis relieves a lot of cognitive burden on the scientist.
Basic numeric structures like arrays forgo tracking this information for efficiency.

Having a stuctured collection of objects which are aligned to the same set of labels allows for a number of higher order interactions.
This includes maintaining relationships between the objects through metadata (e.g. observation and variable loadings of a PCA, the distance matrix a weighted representation was derived from).
We can also keep further annotations on the dataset, e.g. colors associated with categorical labels, so these are preserved on subsetting \autoref{fig:overview}.

### Kinds of elements

* Pairwise elements/ graph operations
### On disk representation

<!-- figure out how to cite zarr and hdf5, zarr has zenodo entry here: https://doi.org/10.5281/zenodo.3773449 -->
<!--
Another big advantage is the on-disk represention, which even for pandas DataFrames is not yet resolved in a canonical way.
For instance, there is none of the binary persisent formats are able to represent all entry types of AnnData.
For instance, even such a key data type a categorical data types are not yet represented in the HDF5 format.
Pickled dataframes are able, but they are non-persistent. -->

<!-- Open with that we have a specified format, then discuss why this is important -->
A conflict in saving datasets and their annotations has been standards vs. ease of use. In the R ecosystem, ease of use has taken precidence. Objects are serialized and written to disk. This is problematic since that data cannot be read by another tool, and may become inaccessible even after software updates. If one chooses to use standard formats to represent all elements of a dataset, a set of standards has to be chosen. AnnData has chosen self-describing hierarchichal data formats such as HDF5 and `zarr` for this purpose. AnnData objects can be efficiently saved to disk using standardized formats \autoref{fig:ecosystem}. This means the data is accessible from other programming environments, as opposed to a serialized format like `pickle` or `Rdata`.

By choosing standardized formats, stored data can be accessed from a variety of ecosystems including `python`, `julia`, `R`, `java`, and `javascript`. While the project has tried to stick to standardized formats, there are a few cases where no standards existed within our models. An especially important example of this is sparse array formats, which are critical for efficient processing of scRNA-seq data. To account for this, we define schemas for these types, which specify how these elements can be read from disk to memory. These specifications are versioned and stored in an internal registry. Versioning allows the specifications to evolve with the project while maintaining the ability to access older data.

<!-- Sparse matrix citation? -->

Like the AnnData object itself, the on-disk representations of these objects closely mirrors their in-memory representation.
Compressed sparse matrices (CSR and CSC format) are stored as a collection of three arrays, `data`, `indices`, and `indptr`, while tabular data is stored in a columnar format.

* Reading and writing
    * Usability
        * Structured, common format which holds all information about the dataset
    * Efficiency
        * The field has typically used text based files for storing data. While there are many problems with this, that's a particularly bad way to hold large numeric data. Better solutions exist so we use them.
* Holding annotations
    * Associating semantic information with your dataset in an organized way (labelled arrays).
        * Also in a way which does not incur much overhead.
* Holding computed properties

### Operations

Mostly IO, subsetting, concatenation? Working on futher

## Examples of use

AnnData is widley used in single cell analysis.
### Figure: Examples

![
**Examples:**
\label{fig:examples}
](figures/examples.pdf)

> Illustrative examples of how anndata is used. Include case studies like:
> * scvelo and layers
> * different sets of variables (modality)
> * spatial information (use of unstructured)

# Future directions

The AnnData project is under active development and will have more features. These include, but are not limited to, more advanced out of core access, a split-apply-combine framework, integration with more of the python ecosystem, and interchange with more formats like apache Arrow.

# Author contributions

...

# References
