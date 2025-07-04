# Example use of Precedent finder code

Examples are taken from the Open Reaction Database: https://open-reaction-database.org/client/id/ord-006fb353766e418abd7ffca7f96dda71
and from RXN-Insight publication: https://doi.org/10.1186/s13321-024-00834-z

## Command-line tool

This is an example with the command-line tool that requires an atom-mapped reaction SMILES as input

and the reaction was mapped with the `rxnmapper` tool.

Run Precent finder using

    bash run.sh

this will produce as tab-separated CSV file with the precedents.

## Notebook

The Jupyter notebook is provided an example for viewing the results.

It will use `rxnmapper` to assigned atom-mapping within the notebook.

One of the examples are provided in the publication.