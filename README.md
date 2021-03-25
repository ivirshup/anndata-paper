# anndata-paper

anndata joss paper

## Building PDF

A PDF is built on every commit as a github action. To get a pdf of the most recent commit click the ✅, and download the artifact.

PDF can be built locally using a docker image:

```sh
docker run --rm \
    --volume $PWD:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/paperdraft
```

## Formatting

[JOSS submission guide](https://joss.readthedocs.io/en/latest/index.html)

The paper should be in `paper.md`, written in markdown. Citations for in `paper.bib` as they would for `bibtex`. We can check if the paper compiles using [this online tool].

For submission, the paper must be in the same repository as the source code. I think using a seperate repo for the writing process will be easier.
