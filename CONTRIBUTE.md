# GlyLES Contribution Guide

## Code

In order to fix bugs or request features, you have two options. Creating an issue in the issue-tracker of the GitHub 
repository or submitting a pull request with a fix for the bug of an implementation of the feature you'd like to have. 
Both options have advantages and disadvantages

| Method       | Advantage     | Disadvantage                  |
|--------------|---------------|-------------------------------|
| Issue        | easy          | fix not immediately available |
| Pull request | fast solution | includes coding               |

### Open an Issue

When opening an issue in the [Issue-Tracker](https://github.com/kalininalab/GlyLES/issues), please indicate what is 
wrong and provide a small example to reproduce your problem. In non-obvious cases of wrong behaviour, please also 
provide/link an explanation of why tools behaviour is wrong.

You can also open an issue if you have a feature-request. In this case, please provide a short example of how the 
requested feature would look like for an exemplary glycan.

### Create a Pull Request

The pull request should contain a comment on what feature has been added or which bug has been removed. For this, 
please include a pytest-based test-case that highlights the bug. That means, you provide a test that fails in the 
current implementation but passes after applying your changes. In non-obvious cases of wrong behaviour, 
please also provide/link an explanation of why the new behaviour is correct and why the old behaviour was wrong.

You can also open a pull request for a new feature implementation. Feature requests should always be opened to the 
`dev` branch of the repository. In case your bug fix or feature implementation requires a new dependency, please update 
the `pyproject.toml` accordingly or state that package and its used version in the comment of the pull request.

If you need help with opening pull requests: These two articles might provide help for you:
  - [Creating a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request)
  - [Creating a pull request from a fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork)

#### Circumvent poetry

If you want to contribute to this project but don't want to install another package manager, you can use poetry within 
a conda environment by executing the following command in the root directory of this project:

```shell
conda create --name glyles python=3.7
conda activate glyles
conda install poetry
poetry install
```

Now, you have a conda environment with all the dependencies you need for development of GlyLES.

## Documentation

In order to update the documentation, please you need to install `sphinx` and `sphinx_rtd_theme` to build the 
documentation by executing

```shell
rm -rf build & sphinx-build ./ ./build -a
```

in the `docs` folder. This will first remove the old documentation build and create a completely new one. Removing is 
necessary as in some cases malfunctions accumulate and the number of warnings increase and the documentation even might 
not build successfully.