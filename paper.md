---
title: 'anndata: Annotated Data Matrices'
bibliography: paper.bib
---

# Summary

<!-- A high level overview of the paper for a non-technical audience. -->

AnnData is a structured representation of high dimensional datasets. It was designed to efficiently represent large datasets in a user friendly way, and one that builds on existing standards rather than inventing it's own.

NGS datasets end up being represented as a matrix of values. These are scalar values of the probed variables for each observation in the dataset.
Recent advances in single cell methods mean that the number of observations in each study has exploded, along with a greater sparsity of values for each cell.

<!-- Why do this in python/ what's the difference from SingleCellExperiment. -->

Since the last advance in high throughput molecular profiling (micro-array -> rna-seq), Python has emerged as an extremely popular language for data analysis and machine learning.
While the analysis of bio-molecular data has previously been focussed on in the R language,

Machine learning/ data analysis methods in python work well with data "shaped-like" data from scRNA-seq (e.g. methods in scikit-learn).
By making it easy to handle this data in python, we can more easily take advantage of these libraries. While these libraries have great computational tools, they frequently work with unlabelled numpy arrays/ scipy sparse arrays.

What does AnnData provide:

<!-- 
I think a big part of the value proposition of AnnData is that the representation works well with the kinds of operations we want to do with single cell data. 
It fits the semantics of the problem well. How do I describe these semantics.
-->

## Manipulating scRNA-seq datasets

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
* Efficiency is key
    * Lengths have been taken to avoid copying uneccesary data
    * Lazyness in views
    * Sparse data

# Statement of need

* Reading and writing
    * Usability
        * Structured, common format which holds all information about the dataset
    * Efficiency
        * The field has typically used text based files for storing data. While there are many problems with this, that's a particularly bad way to hold large numeric data. Better solutions exist so we use them.
* Holding annotations
    * Associating semantic information with your dataset in an organized way (labelled arrays).
        * Also in a way which does not incur much overhead.
* Holding computed properties
