# Contributing

To help other PICRUSt users with the software, see the [picrust-users](https://groups.google.com/forum/?#!forum/picrust-users) mailing list. If you'd like to contribute code and/or documentation, this document will outline how to get started.

## Getting Started

First, [fork](https://help.github.com/articles/fork-a-repo/) the repository and clone it (replace `<your-username>` with your Github username:

```
$ git clone git@github.com:<your-username>/picrust.git
$ cd picrust/
```

Then, set up a development environment using [`conda`](https://conda.io/miniconda.html):

```
$ conda env create -f dev-environment.yml
$ source activate picrust-dev
$ pip install --no-deps --editable .
```

Run the tests to verify:

```
$ pytest
```
