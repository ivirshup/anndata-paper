---
title: "anndata: Annotated data"
authors:
  - name: Isaac Virshup
    orcid: 0000-0002-1710-8945
    affiliation: "1,2,†"
  - name: Sergei Rybakov
    orcid: 0000-0002-4944-6586
    affiliation: "2"
  - name: Philipp Angerer
    orcid: 0000-0002-0369-2888
    affiliation: "2,†,‡"
  - name: F. Alexander Wolf
    orcid: 0000-0002-8760-7838
    affiliation: "2,†,‡"
affiliations:
 - name: University of Melbourne.
   index: 1
 - name: "Helmholtz Munich."
   index: 2
 - name: Corresponding authors.
   index: †
 - name: "Present address: Cellarity, Cambridge, MA."
   index: ‡   
date: November 10, 2021
bibliography: [./paper.bib]
---

# Summary

anndata is a Python package for handling annotated data matrices in memory and on disk.
It is positioned between pandas and xarray by providing structure that organizes data matrix annotations.
anndata offers a broad range of computationally efficient features including, among others, sparse data support, lazy operations, and a PyTorch interface.


# Statement of need

Generating insight from high-dimensional data typically works through training models that annotate observations and variables via low-dimensional representations.
In particular in exploratory workflows, this involves iteratively training models on pre- and post-learned annotations of a data matrix requiring to book-keep both its annotations and learned representations.
anndata offers a canonical data structure for this, which is neither addressed by pandas [@McKinney2010], nor xarray [@Hoyer2017], nor commonly-used modeling packages like sckit-learn [@Pedregosa2011].


# Introduction

Since its initial publication as part of Scanpy [@Wolf2018], anndata matured into an independent software project and became widely adopted (674k total PyPI downloads & 40k downloads/month, 220 GitHub stars & 551 dependent repositories).

So far, anndata has been particularly useful for data analysis in computational biology where advances in single-cell RNA sequencing (scRNA-seq) have given rise to new classes of analysis problems with a stronger adoption of Python over the traditional R ecosystem.
Previous bulk RNA datasets had few observations with dense measurements while more recent scRNA-seq datasets come with high numbers of observations and sparse measurements, both in dimensions of 20k and more.
These new data profit much from the application of scalable machine learning tools of the Python ecosystem.


# The AnnData object

`AnnData` was inspired by similar data structures within the R ecosystem, in particular, `ExpressionSet` [@Huber2015], and the more recent `SingleCellExperiment` [@amezquita2020].

Within the pydata ecosystem, the closest package amenable to store an annotated data matrix is xarray [@Hoyer2017], which enables to deal with labelled data tensors of arbitrary dimensions.
By contrast, the highly popular package pandas [@McKinney2010] operates on single data matrices represented as `DataFrame` objects.
anndata is positioned in between anndata and xarray by providing structure that organizes data matrix annotations.

![**Structure of the AnnData object.**
**a,** The AnnData object is a collection of arrays aligned to the common dimensions of observations (`obs`) and variables (`var`).
Here, color is used to denote elements of the object, with "warm" colors selected for elements aligned to the observations and "cool" colors for elements aligned to variables.
The object is centered around the main data matrix `X`, whose two dimensions correspond to observations and variables respectively.
Primary labels for each of these dimensions are stored as `obs_names` and `var_names`.
If needed, `layers` stores matrices of the exact same shape as `X`.
One-dimensional annotations for each dimension are stored in dataframes `obs` and `var`.
Multi-dimensional annotations are stored in `obsm` and `varm`.
Pairwise relationships are stored in `obsp` and `varp`.
Unstructured data which doesn’t fit this model, but should stay associated to the dataset are stored in `uns`.
**b,** Let us discuss a few examples. The response variable ŷ learned from X is stored as a one-dimensional annotation of observations.
**c,** Principal components and the transformed dimensionality-reduced data matrix obtained through PCA can be stored as multi-dimensional annotations of variables and observations, respectively.
**d,** A k-nearest neighbor graph of any desired representation is stored as a sparse adjacency matrix of pairwise relationships among observations in `obsp`. This is useful to have easy access to the connectivities of points on a low-dimensional manifold.
**e,** Subsetting the `AnnData` object by observations produces a view of data and annotations.
\label{fig:overview}
](figures/overview.pdf)


## The data structure

Standard data structures facilitate data science, with one of the most adopted standards being *tidy data* [@Wickham2014]. anndata defines a standard data structure that makes use of conserved dimensions between data matrix and annotations, similar to `ExpressionSet` [@Huber2015]. With that, `AnnData` makes a particular choice for data organization that has been left unaddressed by packages like scikit-learn or PyTorch [@Paszke2019], which model input and output for each computation as unstructured sets of tensors. Furthermore, `AnnData` offers an on-disk representation that allows sharing data and structured analysis results in form of learned annotations.

Being inspired by `ExpressionSet`, `SingleCellExperiment`, and machine-learning-ready provisioning of data, we note that `AnnData` also complies with *tidy data* standard [@Wickham2014]. However, ...

At the core of `AnnData` is the measured data matrix from which we wish to generate insight (`X`). Each data matrix element stores a value and belongs to an observation in a row (`obs_names`) and a variable in a column (`var_names`), following the *tidy data* standard.
One builds an understanding of the data matrix by annotating observations and variables using `AnnData`'s fields (\autoref{fig:overview}):

* Annotations that can be stored in a single vector get added to the main annotation `DataFrames` for each axis, `obs` and `var`.
* Multi-dimensional representations get added to `obsm` and `varm` and graph-like relations among observations are added to `obsp` and `varp`.

Prior annotations of observations will often denote experimental groups, while derived annotations of observations might be summary statistics, cluster assignments, low-dimensional representations or manifolds. Annotations of variables will often denote alternative names or feature importance measures.

In the context of *tidy data* [@Wickham2014; @Codd1990], one can think of `X` as grouping the data of a specific set of measured variables of interest in an analysis, typically high-dimensional *measured* data in an experiment. Other tables aligned to the observations axis in `AnnData` are then available to store *fixed* data of the experiment, often called metadata, or derived data.


## The data analysis workflow

Let us illustrate how `AnnData` supports workflows of iteratively learning representations and scalar annotations in exploratory data analysis.
For instance, training a clustering, classification or regression model on raw data in `X` produces an estimate of a response variable ŷ. This derived vector is conveniently kept track off by adding it as an annotation of observations (`obs`, \autoref{fig:overview}b).
A reduced dimensional representation obtained through, say Principal Component Analysis or any bottleneck layer of a machine learning model, would be stored as multi-dimensional annotation (`obsm`, \autoref{fig:overview}c).
Storing low-dimensional manifold structure within a desired reduced representation is achieved through a k-nearest neighbor graph in form of a sparse adjacency matrix: a matrix of pairwise relationships of observations (`obsp`, \autoref{fig:overview}d).
Subsetting the `AnnData` object by observations produces a memory-efficient view of data and annotations (\autoref{fig:overview}e).


## The efficiency of data operations

Due to the increasing scale of data, emphasis has been placed on providing efficient data handling operations with low memory and runtime overhead.
To this end, AnnData offers sparse data support, out of core conversions between dense and sparse data, lazy subsetting ("views"), per-element operations for low total memory usage, in-place subsetting, combining AnnData objects with various merge strategies, lazy concatentation, batching, and a backed out-of-memory mode.

In particular, `AnnData` takes great pains to support efficient operations with sparse data. While there currently is no equivalent API for working with sparse and dense data in the python ecosystem, `AnnData` abstracts over the differing existing APIs making it much easier for novices to handle each.
This concerns handling data both on-disk and in-memory with operations for out-of-core access.
A noteworthy design choice means is that we do not follow columnar (or "variable major") data storage as for tidy-data and `DataFrames`.
Our access patterns to X are typically row based, so we use CSR and C order arrays (or "observation major"), which allows efficiently accessing batches of the dataset, to meet the needs of batched learning algorithms.


## The on disk representation

An `AnnData` object captures a unit of the data analysis workflow that groups prior and derived data together.
Providing a persistent and standard on disk format for this unit relieves the pain of working with many competing formats for each individual element and aids reproducibility.
This is particularly needed as even pandas `DataFrames` have no canonical persistent data storage format. `AnnData` has chosen the self-describing hierarchical data formats HDF5 [@collette14] and zarr [@zarr] for this purpose (\autoref{fig:ecosystem}), which are compatible with non-Python programming environments. The broad compatibility and high stability of the format led to wide adoption, and initiatives like the Human Cell Atlas and HubMAP distribute their single-cell omics datasets through `.h5ad`.

![**AnnData provides broad interoperability with tools and platforms.**
`AnnData` objects can be created from a number of formats, including common delimited text files, or domain-specific formats like `loom` files or `CellRanger` outputs.
Once in memory, AnnData provides an API for handling annotated matrices, proving a common base object used by the Python APIs of a range of analytic computational biology tools and integrating well with the APIs of the established Python machine learning ecosystem.
The in memory format has a one-to-one relationship with its hierarchical on disk formats (mapping of elements indicated by color) and uses language-independent technologies, enabling the use by non-Python applications and interchange with other ecosystems.
\label{fig:ecosystem}
](figures/ecosystem.pdf)

Within HDF5 and zarr, we could not find a standard for sparse arrays and DataFrames.
To account for this, we define a schema for these types, which specify how these elements can be read from disk to memory.
These specifications are versioned and stored in an internal registry, which evolves with the project while maintaining the ability to access older data. Like the `AnnData` object itself, the on-disk representations of these types closely mirror their in-memory representation.
Compressed sparse matrices (CSR and CSC format) are stored as a collection of three arrays, `data`, `indices`, and `indptr`, while tabular data is stored in a columnar format.


# The ecosystem

Over the past 5 years, an ecosystem of packages that are built around anndata has grown. This ecosystem is highly focused on scRNA-seq (\autoref{fig:ecosystem}), and ranges from Python APIs [@Zappia2021] to user-interface-based applications [@Megill2021]. Also tools that are not designed around anndata, like scikit-learn and UMAP [@mcinnes2020], nonetheless integrate seamlessly with anndata-based workflows. Since releasing the PyTorch data loader interface `AnnLoader` and the lazy concatenation structure `AnnCollection`, `anndata` offers native ways of integrating into the Pytorch ecosystem in addition to the integration that scvi-tools [@Gayoso2021] offers.

Through the language-independent on-disk format `h5ad`, interchange of data with non-Python ecosytems is easily possible. For analysis of scRNA-seq data in R, this has been particularly simplified by anndata2ri, which allows conversion to `SingleCellExperiment` [@amezquita2020] and Seurat's format [@Hao2020].

![
**Examples of how AnnData is used by packages in the eco system.**
**a,** Squidpy uses AnnData objects for working with spatial data. The coordinates of each sample are stored as an array in `obsm`, an image to overlay the plot on is stored in `uns`, and spatial graph representation in `obsp`.
**b,** Multiple modalities can be represented in a single anndata objects. The variable axis now corresponds to the union of the features across the modalities, modality specific or joint embeddings are stored as separate elements in `obsm` or `obsm`, while inter-modality relations can be stored as graphs in `varp`.
**c,** The `AnnData` model allows for representing RNA velocity analyses by storing counts of different splicing states as separate layers, with velocity based directed graphs in `obsp`.
\label{fig:examples}
](figures/examples.pdf)

Let us discuss examples of anndata's ecosystem for spatial transcriptomics, multiple modalities, and RNA velocity (\autoref{fig:examples}).
In spatial transcriptomics, each high-dimensional observation is annotated with spatial coordinates.
Squidpy [@Palla2021] uses `AnnData` to model their data by storing spatial coordinates as an array (`obsm`) and a spatial neighborhood graph (`obsp`), which is used to find features which are spatially correlated (\autoref{fig:examples}a).
In addition, values from the high dimensional transcriptomic measurement can be overlayed on an image of the sampled tissue, where the image array is stored in `uns`, or handled externally.

`AnnData` can be used to model multimodal data beyond exploiting `AnnData`'s available fields.
For instance, analyses specific to each modality are carried out on separate `AnnData` objects (\autoref{fig:examples}b).
From this, analyses inferring or annotations interactions between modalities can be stored as graphs in `varp`.
Analyses using information from both modalities, like a joint manifold [@Hao2020], can be stored in `obsp`.
A related approach for modelling multimodal data has been utilized by the `muon` package [@Bredikhin2021].
Here, a new `MuData` object is defined, which is essentially a collection of `AnnData` objects, one for each modality measured.
Annotations shared across modalities are stored for the observations for the whole object.
This structure extends to the on-disk format where individual `AnnData` objects are stored as discrete elements inside the `MuData`'s `h5mu` files.
This approach significantly differs from the previous approach by allowing for disjoint sets of observations measured for each modality but is quite similar to `MultiAssayExperiment` within the Bioconductor ecosystem [@Ramos2017].

`AnnData` has been used to model data for fitting models of RNA velocity [@Bergen2020] exploiting the `layers` slot to store a set of matrices for different types of RNA counts (\autoref{fig:examples}c).


# Outlook

The anndata project is under active development towards a variety of features: more advanced out-of-core access, better cloud & relational database integration, a split-apply-combine framework, and interchange with more formats, like Apache Arrow. Furthermore, anndata engages with projects that aim at building out infrastructure for modeling multi-modal data and representing non-homogeneous data to enable learning from Electronic Health Records [@Heumos2021]. Finally, we aim at further extending anndata's data model by interfacing with scientific domain knowledge and data provenance tracking.


# Acknowledgements

I.V. is grateful to Christine Wells for consistent support and freedom to pursue work on anndata and Scanpy.
We are grateful to Ryan Williams and Tom White for contributing code related to zarr.
We thank Jonathan Bloom for contributing a comprehensive PR on group-by functionality.
We are grateful to Fabian Theis & lab for continuing dissemination along with Scanpy over the past years.
F.A.W. and P.A. thank Cellarity for supporting continued engagement with open source work.

# Author contributions

Isaac has led the anndata project since v0.7, and contributed as a developer before. He rewrote wide parts of the code introducing high efficiency, robustness, concatenation, ...
Sergei made diverse contributions to the code base, in particular, the first version of `layers`, benchmarking and improvement of the earlier versions of the IO code, the PyTorch dataloader `AnnLoader` and the lazy concatenation data structure `AnnCollection`.
Phil suggested to replace Scanpy's initial unstructured annotated data object to one mimicking R's ExpressionSet, and wrote AnnData's [first implementation](https://github.com/theislab/scanpy/commit/315859c5586116434ea3b7ce97512a5e2a1030e2) with indexing, slicing, ... and ascertained good software practices in the project. ...
Alex led the project until 0.7, introduced centering data science workflows around an [initially unstructured annotated data object](https://github.com/theislab/scanpy/tree/c22e48abe45a6ccca5918bbf689637caa4b31250), designed the API, wrote tutorials and documentation until 0.7, and implemented most of the early functionality, among others, reading & writing, the on-disk format `h5ad`, sparse data support, backed mode, views.

# References
