---
title: 'anndata: annotated data'
authors:
  - name: Isaac Virshup
    orcid: 0000-0002-1710-8945
    affiliation: "1"    
  - name: Sergei Rybakov
    affiliation: "2"
  - name: Philipp Angerer
    affiliation: "2, 3, †"    
  - name: F. Alexander Wolf
    orcid: 0000-0002-8760-7838
    affiliation: "2, 3, †"
affiliations:
 - name: University of Melbourne
   index: 1
 - name: Helmholtz Munich, Institute of Computational Biology, Munich, Germany.
   index: 2
 - name: Cellarity, Cambridge, MA.
   index: 3
 - name: Equal contributions
   index: †
date: 1 October 2020
bibliography: paper.bib
---

# Summary

anndata is a python software package for handling annotated datasets.

# Introduction

Generating insight from high-dimensional data typically requires two steps. Condense them into efficient representations and assign semantic labels. Both is achieved by learning patterns in high dimensions based on existing annotations. Once insight is gained, knowledge is updated and assigned in form of learned annotations. This defines the general machine learning and data science paradigm, which, which anndata, we strive represent in a data structure that integrates well with the pytdata ecosystem.

Within the pydata ecosystem the closest package that would be amenable to serve this paradigm is xarray `[@xarray_author:YYYY]`, which enables to deal with highly complex labelled data tensors of arbitrary dimensions. On the other hand, there is the highly popular package pandas`[@pandas_author:YYYY]`, which merely provides and operates on `DataFrames`, that is, single tables of data. anndata is positioned in between by providing the minimal additional structure to enable storing compact annotations and representations of high-dimensional data, making the book keeping during learning from it much simpler.

With that, anndata perfectly integrates into scikit-learn `[@scikit_author:YYYY]`, statsmodels `[@statsmodels_author:YYYY]`, seaborn `[@waskom:2016]`, and easily interfaces with pytorch and tensorflow.

Specifically, the central `AnnData` class stores observations (samples) of variables (features) in the rows of a matrix. This is the convention of the modern classics of statistics [Hastie09] and machine learning [Murphy12], the convention of dataframes both in R and Python and the established statistics and machine learning packages in Python (statsmodels, scikit-learn).

Machine learning/ data analysis methods in python work well with data "shaped-like" data from scRNA-seq (e.g. methods in scikit-learn).
By making it easy to handle this data in python, we can more easily take advantage of these libraries. While these libraries have great computational tools, they frequently work with unlabelled numpy arrays/ scipy sparse arrays.

AnnData is a structured representation of high dimensional datasets. It was designed to efficiently represent large datasets in a user friendly way, and one that builds on existing standards rather than inventing it's own.

NGS datasets end up being represented as a matrix of values. These are scalar values of the probed variables for each observation in the dataset.
Recent advances in single cell methods mean that the number of observations in each study has exploded, along with a greater sparsity of values for each cell.

<!-- Why do this in python/ what's the difference from SingleCellExperiment. -->

Since the last advance in high throughput molecular profiling (micro-array -> rna-seq), Python has emerged as an extremely popular language for data analysis and machine learning.
While the analysis of bio-molecular data has previously been focussed on in the R language.

<!-- 
I think a big part of the value proposition of AnnData is that the representation works well with the kinds of operations we want to do with single cell data. 
It fits the semantics of the problem well. How do I describe these semantics.
-->

# What is anndata?

## General features

* Reading and writing
    * Usability
        * Structured, common format which holds all information about the dataset
    * Efficiency
        * The field has typically used text based files for storing data. While there are many problems with this, that's a particularly bad way to hold large numeric data. Better solutions exist so we use them.
* Holding annotations
    * Associating semantic information with your dataset in an organized way (labelled arrays).
        * Also in a way which does not incur much overhead.
* Holding computed properties

## scRNA-seq related features

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

# References