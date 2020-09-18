---
title: 'anndata: Annotated Data Matrices'
bibliography: paper.bib
---

# Summary

AnnData is a structured representation of high dimensional datasets. It was designed to efficiently represent large datasets in a user friendly way.

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
