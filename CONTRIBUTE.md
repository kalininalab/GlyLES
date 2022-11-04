# GlyLES Contribution Guide

## Code

In order to fix bugs, you can have two options. Creating an issue in the issue-tracker of the github repository or 
submitting a pull request with a fix for the bug. Both options have an advantage and a disadvantage

| Method       | Advantage    | Disadvantage                  |
|--------------|--------------|-------------------------------|
| Issue        | easy         | fix not immediately available |
| Pull request | fast solving | includes coding               |

### Open an Issue

When opening an issue in the [Issue-Tracker](https://github.com/kalininalab/GlyLES/issues), please indicate what is 
wrong and provide a small example to reproduce your problem. In non-obvious cases of wrong behaviour, please also 
provide/link an explanation of why tools behaviour is wrong.

### Create a Pull Request

The pull request should contain a comment on what feature has been added or which bug has been removed. For this, 
please include a pytest-framework based testcase that highlights the bug. In non-obvious cases of wrong behaviour, 
please also provide/link an explanation of why the new behaviour is correct and why the old behaviour was wrong.

## Documentation

In order to update the documentation, please you need to install `sphinx` and `sphinx_rtd_theme` to build the 
documentation by executing

```shell
rm -rf build & sphinx-build ./ ./build -a
```

in the `docs` folder. This will first remove the old documentation build and create a completely new one. Removing is 
necessary as in some cases malfunctions accumulate and the number of warnings increase and the documentation even might 
not build successfully.