# nlmixr2 2.0.7

* nlmixr2 now re-exports `logit` so that certain models will work with
  a simple `library(nlmixr2)` instead of
  `library(rxode2);library(nlmixr2)`

* `vpcSim()` now exports the new `nretry` option for more robust
  control of `vpcSim()`

* Update documentation to mention the package names that work with
  nlmixr2, like `xpose.nlmixr2` instead of `xpose.nlmixr`

* Manual hard re-export of `nlmixr2est::nlmixr2` to allow `pkgdown` to
  document this function.

# nlmixr2 2.0.6

* nlmixr2 is an umbrella package to include the lower level packages
  `rxode2`, `nlmixr2est`, `nlmixr2extra`, and `nlmixr2plot`

* Added a `NEWS.md` file to track changes to the package.
