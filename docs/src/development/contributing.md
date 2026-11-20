# Development and documentation

Activate the repository root for package work. Keep the documentation environment
separate so tutorial-only dependencies do not become runtime dependencies.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

`Pkg.test()` activates test extras, including OptimizationOptimJL. Running
`test/runtests.jl` directly from the root environment does not activate those extras.

## Build the site

From the repository root, set up the docs environment once:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
```

Install Node.js 20 or newer with npm, then build and preview:

```sh
julia --project=docs docs/make.jl
npm --prefix docs run docs:preview
```

The build executes Markdown examples and doctests, checks exported docstrings,
generates Markdown, runs `npm ci`, and builds VitePress. It does not deploy by
default. `DOCS_DEPLOY=true` enables deployment and is set only in the documentation
workflow; Documenter also checks the GitHub event before publishing.
For a production dev artifact without publishing, run
`DOCS_DEPLOY_LINKS=true julia --startup-file=no --project=docs docs/make.jl`.
This includes version-navigation scripts while leaving deployment disabled.

For live editing, run these commands from the repository root:

```sh
npm --prefix docs ci
npm --prefix docs run docs:dev
```

The active single-source tutorial scripts live in `docs/examples/tutorials`.
`docs/make.jl` generates their pages under `docs/src/tutorials` with Literate.
Full workflows are rendered as non-executing Julia; small explicitly named
`@example` waveform plots run during normal builds. Running either script directly
executes its original-size model, scans and calibration, so do so deliberately.
To execute both complete tutorials and embed every calibration, population and
Floquet figure, stop `docs:dev` and run:

```sh
DOCS_RUN_TUTORIALS=true julia --startup-file=no --project=docs docs/make.jl
npm --prefix docs run docs:preview
```

This explicitly runs the original-size models, long trajectories, scans and
calibrations. `DOCS_RUN_TUTORIALS=true` cannot be combined with `DOCS_DRAFT=true`:
draft mode skips examples, so the build rejects that combination. The default
build remains inexpensive and includes only the small waveform tutorial plots.
The earlier shortened archives remain unchanged in `docs/tutorials`.

`docs:dev` watches `docs/src` and active tutorial Julia sources and regenerates Markdown in Documenter's draft mode
before refreshing the browser. Julia formats the documentation, but examples and
doctests do not execute; calculated outputs and plots are omitted. Save a source
file and allow generation to finish before expecting the update. Stop with Ctrl+C.
From inside `docs/`, use `npm run docs:dev` instead.

To serve the last generated docs without regeneration, use
`npm --prefix docs run docs:serve`. This watches `docs/build/.documenter` only.
Run the full Julia build above when you need updated example results or docstrings
from package source. `DOCS_MARKDOWN_ONLY=true` skips the website build but still
executes examples unless you also set `DOCS_DRAFT=true`.

## Source organization

Keep explicit overloads in `overloads.jl` beside the structures they extend.
`src/gates/evolution` owns solver-based evolution, and gate calibration/metrics
live in `src/gates/calibration`. Hamiltonians use frequency units; apply `2π` at
direct solver calls using `make_hamiltonian(...; scalar=2pi)`; built-in `evolve` supplies it internally.

## Adding documentation

- Add API docstrings next to their definitions, including units and return values.
- Render each supported API in one reference `@docs` block. Build coverage is strict.
- Put runnable tutorial code in named `@example` blocks. Use assertions for expected behavior and small deterministic models.
- Keep output plots and temporary bundles in the build or temporary directories.
- Add pages to the navigation in `docs/make.jl`.
- Run the production build, package tests for source changes, and `git diff --check`.

CI builds pull requests and deploys the `dev` branch and version tags to
`gh-pages`. Every valid version tag appears in the version picker, including
prereleases. The npm lockfile is committed so site dependencies are reproducible.

## Release workflow

Ongoing work happens on `dev`; `main` records released commits. To publish a
release, merge the validated `dev` work into `main`, set `Project.toml` to the
release version, and run the package tests and docs build. Tag that exact `main`
commit with `v` followed by the project version (for example,
`v0.1.0-alpha.1`), push the tag, and create a GitHub Release. Mark alpha and
beta releases as prereleases. The tag build publishes its own docs, while `dev`
continues to publish at `/dev/`. After release, advance `dev` to the next
development version.

Docs watcher verification marker 20261008.
