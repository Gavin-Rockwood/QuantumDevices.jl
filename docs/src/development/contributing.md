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

For editing, regenerate the Markdown then start the VitePress development server:

```sh
DOCS_MARKDOWN_ONLY=true julia --project=docs docs/make.jl
npm --prefix docs ci
npm --prefix docs run docs:dev
```

Re-run generation after changing source Markdown or docstrings. The scripts all
use `docs/build/.documenter`; do not run VitePress against the source directory.

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
