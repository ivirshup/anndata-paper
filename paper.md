---
title: "anndata: Annotated data matrices"
authors:
  - name: Isaac Virshup
    orcid: 0000-0002-1710-8945
    affiliation: "1,†"
  - name: Sergei Rybakov
    affiliation: "2"
  - name: Philipp Angerer
    affiliation: "2,†,‡"
  - name: F. Alexander Wolf
    orcid: 0000-0002-8760-7838
    affiliation: "2,†,‡"
affiliations:
 - name: University of Melbourne.
   index: 1
 - name: "Helmholtz Munich, Institute of Computational Biology."
   index: 2
 - name: Corresponding authors.
   index: †
 - name: "Present address: Cellarity, Cambridge, MA."
   index: ‡   
date: November 1st, 2021
bibliography: paper.bib
---

# Summary

anndata is a Python package for handling annotated data matrices in memory and on disk.


# Statement of need

In exploratory data analysis -- say, based on sckit-learn [@Pedregosa2011] -- generating insight from high-dimensional data is typically achieved through training models to learn patterns that allow to condense data into meaningful low-dimensional representations and assign meaning to observations and variables.
This involves workflows of iteratively training models on pre- and post-learned annotations of data, requiring to book-keep their representations and scalar annotations (labels and numericals).
anndata's purpose is to make such workflows as efficient as possible through a data structure that naturally integrates book-keeping with model training and analysis, well-integrated into the pydata ecosystem. Neither pandas nor xarray meet this need.

anndata turned out to be particularly useful for data analysis in computational biology, where advances in single-cell RNA sequencing (scRNA-seq) have given rise to new classes of analysis problems.
While previous bulk RNA datasets typically have few observations with dense measurements, scRNA-seq datasets typically come with high numbers of observations with very sparse measurements, all with dimensions of 20k and more.
Generating insight from these new data profits much from the application of scalable machine learning tools, which are more abundant in the Python ecosystem than in the R ecosystem.
To nontheless enable use of the wealth of computational biology tools in the R ecosystem, anndata offers a cross-ecosystem on-disk format: h5ad.

Since its initial publication as a part of Scanpy [@Wolf2018], anndata matured into an independent software project with significant progress having been made in the past years, and long-term expansion plans for upcoming years. The present paper introduces anndata in this entirety.


# The AnnData object

`AnnData` was inspired by similar data structures within the R ecosystem, in particular, `ExpressionSet`, and the more recent `SingleCellExperiment`.

Within the pydata ecosystem, the closest package that would be amenable to store an annotated data matrix is xarray [@Hoyer2017], which enables to deal with highly complex labelled data tensors of arbitrary dimensions.
By contrast, the highly popular package pandas [@McKinney2010] operates on single data matrices represented as `DataFrame` objects.
anndata is positioned in between xarray and pandas by providing the minimal additional structure that enables storing annotations of a data matrix, through sharing their conserved dimensions.


## The data structure

Consistent data structures facilitate both exploratory data analysis and sharing of data by saving the time to translate between differing data structures.
For instance, the tidyverse project [@Wickham2014] of the R ecosystem defined a successful consistent data standard for an entire field.

By making use of conserved dimensions between data matrix and annotations, `AnnData` makes a particular choice for data organization that has been left unaddressed by packages like scikit-learn or PyTorch, which model input and output for each computation as unstructured sets of tensors. Furthermore, `AnnData` offers an on-disk representation that allows sharing data and structured analysis results in form of learned annotations.

At the core of `AnnData` is the measured data matrix from which we wish to generate insight (`X`). Each data matrix element belongs to an observation (`obs_names`) and a variable (`var_names`) and contains a value (which can be "missing", like `nan`).
We build our understanding of the data matrix by adding annotated and derived values onto observations and variables \autoref{fig:overview} using `AnnData`'s canonical locations:
Simple annotations and derived values that can be stored in a single vector are added to the main annotation `DataFrames` for each axis, `obs` and `var`.
Multi-dimensional representations are added to `obsm` and graph-like relations among observations are added to `obsp`.
Annotations of variables include values like alternative names (e.g. different identifier mappings) or categories for each variable.
Annotations of observations at dataset creation may list experimental groups, whereas derived annotations could be descriptive statistics (e.g. mean and variance), cluster assignments, low-dimensional representations (`obsm`) or manifolds (`obsp`).

![**Structure of the AnnData object.**
*(a)* The AnnData object is a collection of arrays aligned to the common dimensions of observations (`obs`) and variables (`var`).
Here, color is used to denote elements of the object, with "warm" colors selected for elements aligned to the observations and "cool" colors for elements aligned to variables.
The object is centered around the main data matrix `X`, whose two dimensions correspond to observations and variables respectively.
Primary labels for each of these dimensions are stored as `obs_names` and `var_names`.
If needed, `layers` stores matrices of the exact same shape as `X`.
One-dimensional annotations for each dimension are stored in dataframes `obs` and `var`.
Multi-dimensional annotations are stored in `obsm` and `varm`.
Pairwise relationships are stored in `obsp` and `varp`.
Unstructured data which doesn’t fit this model, but should stay associated to the dataset are stored in `uns`.
*(b)* The response variable ŷ learned from X is stored as a one-dimensional annotation of observations.
*(c)* Principal components and the transformed dimensionality-reduced data matrix obtained through PCA can be stored as multi-dimensional annotations of variables and observations, respectively.
*(d)* A k-nearest neighbor graph of any desired representation is represented as a sparse adjacency matrix, constituting a pairwise relationship of observations in `obsp`.
*(e)* Subsetting the `AnnData` object by observations produces a view of data and annotations.
\label{fig:overview}
](figures/overview.pdf)


## The data analysis workflow

Let us illustrate how `AnnData` supports workflows of iteratively learning representations and scalar annotations in exploratory data analysis.
For instance, training a clustering, classification or regression model on raw data in `X` produces an estimate of a response variable ŷ. This derived vector is conveniently kept track off by adding it as an annotation of observations (`obs`, \autoref{fig:overview}[b]).
A reduced dimensional representation obtained through, say Principal Component Analysis or any bottleneck layer of a machine learning model, would be stored as multi-dimensional annotation (`obsm`, \autoref{fig:overview}(c)).
Storing low-dimensional manifold structure within a desired reduced representation is achieved through a k-nearest neighbor graph in form of a sparse adjacency matrix: a matrix of pairwise relationships of observations (`obsp`, \autoref{fig:overview}(d)).
Subsetting the `AnnData` object by observations produces a memory-efficient view of data and annotations (\autoref{fig:overview}(e)).


## The efficiency of data operations

Due to the increasing scale of data, emphasis has been placed on providing efficient data handling operations with low memory and runtime overhead.
To this end, AnnData offers sparse data support, out of core conversions between dense and sparse data, lazy subsetting ("views"), per-element operations for low total memory usage, in-place subsetting, combining AnnData objects with various merge strategies, lazy concatentation, batching, and a backed out-of-memory mode.

In particular, `AnnData` takes great pains to support efficient operations with sparse data. While there currently is no equivalent API for working with sparse and dense data in the python ecosystem, `AnnData` abstracts over the differing existing APIs making it much easier for novices to handle each.
This concerns handling data both on-disk and in-memory with operations for out-of-core access.
A noteworthy design choice means is that we do not follow columnar (or "variable major") data storage as for tidy-data and `DataFrames`.
Our access patterns to X are typically row based, so we use CSR and C order arrays (or "observation major"), which allows efficiently accessing batches of the dataset, to meet the needs of batched learning algorithms.


## The on disk representation

An `AnnData` object captures a unit of the data analysis workflow that groups original and derived data together.
Providing a persistent and standard on disk format for this unit relieves the pain of working with many competing formats for each individual element and aids reproducibility.
This is particularly needed as even pandas `DataFrames` have no canonical persistent data storage format, yet, which only starts to get addressed by an improved Parquet interface. Also the R ecosystem has not yet arrived at a fully satisfactory solution, and many tools still serialize in-memory objects to disk.
This is problematic since it prohibits reading data by another tool and is highly non-persisent, meaning, it may become inaccessible even after software updates.

If one chooses to use standard formats to represent all elements of a dataset, a set of standards has to be chosen.
`AnnData` has chosen the self-describing hierarchical data formats HDF5 and zarr [https://doi.org/10.5281/zenodo.3773449] for this purpose (\autoref{fig:ecosystem}), which are compatible with many programming environments.

![**AnnData provides common conventions for data handling for an ecosystem of tools.**
`AnnData` objects can be created from a number of formats, including common delimited text files, or domain-specific formats like `loom` files or `CellRanger` outputs.
Once in memory, AnnData provides an API for handling annotated matrices, proving a common base object used by the Python APIs of a range of analytic tools.
The in memory format has a one to one relationship with its hierarchical on disk formats (mapping of elements indicated by color) and uses language-independent technologies, facilitating use by non-Python applications and interchange with other ecosystems.
\label{fig:ecosystem}
](figures/ecosystem.pdf)

anndata has adopted standardized formats where possible, but could not find a standard for sparse arrays and DataFrames.
To account for this, we define a schema for these types, which specify how these elements can be read from disk to memory.
These specifications are versioned and stored in an internal registry, which allows the specifications to evolve with the project while maintaining the ability to access older data.

Like the `AnnData` object itself, the on-disk representations of these objects closely mirrors their in-memory representation.
Compressed sparse matrices (CSR and CSC format) are stored as a collection of three arrays, `data`, `indices`, and `indptr`, while tabular data is stored in a columnar format.


# The ecosystem

Over the past 5 years, an ecosystem of packages that are built around anndata has grown (as of today: 674k PyPI downloads, 40k downloads/month and 220 GitHub stars). This ecosystem is highly focused on scRNA-seq (\autoref{fig:ecosystem}), and ranges from Python APIs [@Gayoso2021; @Palla2021; Bergen2020; Bredikhin2021] to user-interface-based applications [@Megill2021].

![
**AnnData is used to model multiple data types.**
Examples of how AnnData is used by packages in the eco system.
*(a)* Squidpy uses AnnData objects for working with spatial data. The coordinates of each sample are stored as an array in `obsm`, an image to overlay the plot on is stored in `uns`, and spatial graph representation in `uns`.
*(b)* Multiple modalities can be represented in a single anndata objects. The variable axis now corresponds to the union of the features across the modalities, modality specific or joint embeddings are stored as seperate elements in `obsm` or `obsm`, while inter-modality relations can be stored as graphs in `varp`.
*(c)* The `AnnData` model allows for representing rna velocity analyses by storing counts of different splicing states as separate layers, with velocity based directed graphs in `obsp`.
\label{fig:examples}
](figures/examples.pdf)

Let us illustrate the compatibility for spatial transcriptomics, RNA velocity, and generally work with multiple modalities \autoref{fig:examples}.

In spatial transcriptomics, each observation has spatial coordinates associated with it.
Squidpy [@Palla2021] uses `AnnData` to model their data by storing spatial coordinates as an array in `obsm`.
These coordinates are used to create a spatial graph, which is used to find features which are spatially coorelated.
Values from the high dimensional experiment can be overlayed on an image of the sampled tissue (where the image array is stored in `uns`, or externally handled). AnnData has been used to model the data for RNA velocity [@Bergen2020].


`AnnData` can also be used to model multimodal data, though multiple approaches for this can be used.
In one approach, analyses specific to each modality are carried out separate `AnnData` objects.
These can then be combined along the variable axis (using the shared observations) for multimodal analyses.
Here an indicator array can be used to separate which variables belong to which modality.
From this, analyses inferring or annotations interactions between modalities can be stored as graphs in `varp`.
Analyses using information from both modalities, like a joint manifold [@Hao2020], can be stored in `obsp`.

Another approach for modelling multimodal data has been utilized by the `muon` package [@Bredikhin2021].
Here, a new `MuData` object is defined, which is essentially a collection of `AnnData` objects, one for each modality measured.
Annotations shared across modalities are stored for the observations for the whole object.
This structure extends to the on-disk format where individual `AnnData` objects are stored as discrete elements inside the `MuData`'s `h5mu` files.
This approach significantly differs from the previous approach by allowing for disjoint sets of observations measured for each modality.
This approach is quite similar to the MultiAssayExperiment from the bioconductor ecosystem [@Ramos2017].

# Conclusions

The AnnData project is under active development and will have more features. These include, but are not limited to, more advanced out of core access, a split-apply-combine framework, integration with more of the python ecosystem, and interchange with more formats like apache Arrow. Beyond further building out infrastructure for modeling multi-modal data, the ecosystem also works towards being able to represent non-homogeneous data to enable learning from Electronic Health Records.


# Acknowledgements

We are grateful to Ryan Williams for contributing code related to zarr.


# Author contributions

Isaac has led the anndata project since v0.7, and contributed as a developer before. He rewrote wide parts of the code introducing high efficiency, robustness, concatenation, ...
Sergei made diverse contributions to the code base, in particular, the first version of `layers`, benchmarking and improvement of the earlier versions of the IO code, the PyTorch dataloader `AnnLoader` and the lazy concatenation data structure `AnnCollection`.
Phil suggested to replace Scanpy's initial unstructured annotated data object to one mimicking R's ExpressionSet, and wrote AnnData's [first implementation](https://github.com/theislab/scanpy/commit/315859c5586116434ea3b7ce97512a5e2a1030e2) with indexing, slicing, ... and ascertained good software practices in the project. ...
Alex led the project until 0.7, introduced the idea of centering data science workflows around an [initially unstructured annotated data object](https://github.com/theislab/scanpy/tree/c22e48abe45a6ccca5918bbf689637caa4b31250), guided the user design of the first implementation, implemented most of the subsequent early functionality and wrote the documentation until 0.7: reading & writing, paired in-memory manipulation and on-disk format, sparse data support on disk and fast loading, backed mode, views, AnnData's slots beyond `.obs` and `.var`, and numerous smaller design choices to integrate AnnData well into Scanpy's workflows.

# References
