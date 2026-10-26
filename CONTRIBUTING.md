# Contributing to MRIReco.jl

Thanks for helping to improve MRIReco.jl! This guide explains how to ask for help, how to report problems and how to contribute changes. Everyone taking part in this project is expected to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Questions and Bug Reports

If you have a question about using the packages, please ask in [GitHub Discussions](https://github.com/MagneticResonanceImaging/MRIReco.jl/discussions).

If you find a bug or want to suggest a feature, please search the [existing issues](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues) first and open a new one if nothing matches. A good bug report contains

- a short example that we can run to reproduce the problem, ideally with simulated data or a small public dataset,
- the complete error message and stack trace,
- the versions you used, which you can print with

```julia
using Pkg
println("Julia ", VERSION)
for (_, pkg) in Pkg.dependencies()
    startswith(pkg.name, "MRI") && println(pkg.name, " ", pkg.version)
end
```

## Repository Layout

MRIReco.jl is a monorepo. The root directory contains the MRIReco package, and each subpackage lives in its own directory with its own `Project.toml`:

| Package | Purpose |
|:---|:---|
| MRIBase | Data types for raw data, acquisition data and trajectories |
| MRIOperators | Encoding operators |
| MRIFiles | Reading and writing raw data files |
| MRISampling | Undersampling patterns |
| MRICoilSensitivities | Coil sensitivity estimation and coil compression |
| MRISimulation | Simulation of MRI data |
| MRIReco | Image reconstruction |

The subpackages are registered as separate packages, but all of them except MRIBase depend on other packages in this repository. During development, use the local copies of these packages instead of the registered versions, as described below.

## Development Setup

1. Fork the repository on GitHub and clone your fork.
2. Start Julia 1.10 or newer in the repository root and activate the package you want to work on, for example `MRISampling`.
3. Develop the local copies of the packages it needs and install the remaining dependencies:

```julia
using Pkg
Pkg.activate("MRISampling")
Pkg.develop([PackageSpec(path = "MRIBase")])
Pkg.instantiate()
```

The packages to develop are the ones that the CI uses (`deps` in `.github/workflows/ci.yml`). They include packages that are only needed for the tests:

| Package | Path to activate | Paths to develop |
|:---|:---|:---|
| MRIReco | `.` | `MRIBase`, `MRIOperators`, `MRISimulation`, `MRISampling` |
| MRIBase | `MRIBase` | none |
| MRIFiles | `MRIFiles` | `.`, `MRIBase`, `MRICoilSensitivities` |
| MRIOperators | `MRIOperators` | `MRIBase`, `MRISimulation` |
| MRISampling | `MRISampling` | `MRIBase` |
| MRISimulation | `MRISimulation` | `MRIBase`, `MRIOperators` |
| MRICoilSensitivities | `MRICoilSensitivities` | `.`, `MRIBase`, `MRISimulation`, `MRISampling` |

For MRIOperators, developing MRISimulation creates a dependency cycle, because MRISimulation itself depends on MRIOperators. Pkg warns about the cycle and skips precompiling these two packages; the CI uses the same setup.

`Pkg.develop` changes the `Manifest.toml` of the activated package and may also change its `Project.toml`, for example by adding the developed packages. Manifest files are ignored by git, but please don't commit changes that `Pkg` makes to a `Project.toml` unless your pull request is about dependencies.

## Running the Tests

With the environment of a package activated as above, run its tests with

```julia
Pkg.test("MRISampling")
```

The CI tests a pull request in three ways:

- The "All Packages (latest)" workflow tests every package with the local copies from the table above, on the latest Julia release and on Julia 1.10. It shows whether the packages work together in their current state.
- Each package also has its own workflow, which runs when the package changes. It tests the package on Linux, macOS and Windows with the released versions of its dependencies, which shows whether the package can be released on its own.
- Changes to MRIBase and MRIOperators are also tested against the released versions of the packages that depend on them, which shows whether the changes would break them.

The steps of these workflows are in `.github/workflows/`, so you can also run them locally. GPU tests run on Buildkite with NVIDIA and AMD GPUs; the commands are in `.buildkite/pipeline.yml`, and the tests are in `test/gpu/` and `MRIOperators/test/gpu/`.

## Building the Documentation

The documentation is built with [Documenter.jl](https://documenter.juliadocs.org/stable/), and the tutorials and how-to guides in `docs/src/literate` are generated with [Literate.jl](https://fredrikekre.github.io/Literate.jl/v2/). To build it locally, develop all packages into the `docs` environment, as the CI does:

```julia
using Pkg
Pkg.activate("docs")
Pkg.develop([PackageSpec(path = "."), PackageSpec(path = "MRIBase"), PackageSpec(path = "MRIOperators"),
             PackageSpec(path = "MRIFiles"), PackageSpec(path = "MRISampling"),
             PackageSpec(path = "MRICoilSensitivities"), PackageSpec(path = "MRISimulation")])
Pkg.instantiate()
```

Then run `julia --project=docs docs/make.jl` from the repository root and open `docs/build/index.html`. The build runs all examples, so it takes several minutes. Don't commit the changes that `Pkg.develop` makes to `docs/Project.toml`.

## Changing Several Packages

Some changes affect more than one package, for example a new function in MRIBase that MRIFiles uses. Work in the environment of the package that depends on the other one, here MRIFiles. MRIBase is developed in this environment (see the table above), so MRIFiles uses your local copy of MRIBase, and the tests of MRIFiles run with your changes to both packages. Then run the tests of MRIBase in its own environment:

```julia
using Pkg
Pkg.activate("MRIFiles")
Pkg.develop([PackageSpec(path = "."), PackageSpec(path = "MRIBase"),
             PackageSpec(path = "MRICoilSensitivities")])
Pkg.instantiate()
Pkg.test("MRIFiles")

Pkg.activate("MRIBase")
Pkg.instantiate()
Pkg.test("MRIBase")
```

The documentation build also uses the local copies of all packages, so it includes both changes. In the CI, the workflow of MRIFiles tests it with the released version of MRIBase, so it fails until the changes in MRIBase are released. Please mention in the pull request that MRIFiles needs a new release of MRIBase, so the packages can be released in the right order.

## Pull Requests

If you are looking for something to work on, the issues labeled [good first issue](https://github.com/MagneticResonanceImaging/MRIReco.jl/labels/good%20first%20issue) are a good place to start.

- Keep a pull request focused on one fix or feature; smaller pull requests are easier to review.
- Add tests for new functionality and bug fixes.
- Add a docstring to every new exported function or type, and check that the documentation still builds.
- There is no enforced code formatter. Please follow the style of the surrounding code.
- Mention the issue that the pull request addresses, e.g. `Fixes #123`.

Please make sure that the CI passes, or explain in the pull request why a check fails.

## License

By contributing, you agree that your contributions are licensed under the [MIT License](LICENSE.md) of this repository.
