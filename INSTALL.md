# Quick Install

All installation instructions are in more detail [here](https://github.com/picrust/picrust2/wiki/Installation).

First, clone the repository and enter the directory:

```
git clone https://github.com/picrust/picrust2.git
cd
```

Next, install the required dependencies listed [here](https://github.com/picrust/picrust2/wiki/Installation#pre-requisites).


Install PICRUSt2 into it's own environment using [`conda`](https://conda.io/miniconda.html):

```
conda env create -f picrust2-env.yaml
source activate picrust2
pip install --no-deps --editable .
```

Finally, run the tests to verify the install was successful:

```
pytest
```
