# LucasStokey_PBSet
Problem set solution for [PhD course in Fiscal Policy](https://sites.google.com/a/nyu.edu/axelleferriere/teaching) taught by [Axelle Ferriere](https://sites.google.com/a/nyu.edu/axelleferriere/home) [PSE, fall 18]

## Dependencies
```
julia v1.0.2
  Dierckx v0.4.1
  Distributions v0.16.4
  IJulia v1.14.1
  LaTeXStrings v1.0.3
  NLsolve v3.0.1
  PGFPlots v3.0.1
  Plots v0.21.0
  SpecialFunctions v0.7.2
python v3.7
  Jupyter v1.0.0
Pandoc 2.2.1
XeTeX v3-2.6
```

To run the minimal code under the `bin` folder, make sure you have dependencies by running the following line in your Julia REPL
```julia
using Pkg; Pkg.add(["Dierckx", "Distributions", "LaTeXStrings", "NLsolve", "PGFPlots", "Plots", "SpecialFunctions"])
```

You can use another plot backend by changing the first line ```julia using Plots; pgfplots()```. For the default backend, write ```julia using Plots; gr()```.

## Reproducing results
You can download and run the julia notebook in the `src` folder. To build web pages and LaTeX report, run the makefile in the root directory (only on MacOS and GNU/Linux)
