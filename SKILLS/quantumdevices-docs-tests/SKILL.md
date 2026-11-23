---
name: quantumdevices-docs-tests
description: Work on QuantumDevices documentation, tutorial generation, and feature-focused or full package test workflows.
---

Read [documentation and builds](../../docs/src/development/contributing.md) and inspect `test/runtests.jl` and the relevant feature test before selecting checks; new feature test files must be included in `test/runtests.jl`. `Pkg.test()` activates test extras; an ordinary root environment may lack them. For the full four-thread package suite, use:

```sh
julia --startup-file=no --threads=4 --project=. -e 'using Pkg; Pkg.test()'
```

Active tutorial sources are `docs/examples/tutorials/*.jl`; `docs/make.jl` generates their Markdown. A default build executes cheap examples, not full tutorial workflows. `DOCS_MARKDOWN_ONLY=true DOCS_DRAFT=true julia --project=docs docs/make.jl` skips site building and example execution; normal builds use `julia --project=docs docs/make.jl`. The full tutorial build, `DOCS_RUN_TUTORIALS=true`, reuses current whole-tutorial execution caches and runs missing/stale expensive models, scans, and calibrations; use only when requested and ensure draft and Markdown-only flags are unset. `DOCS_EXECUTION=refresh` forces both full tutorials to execute. `DOCS_EXECUTION=preview` uses available full cached results with current/stale/missing labels, executes no examples or doctests, and cannot deploy. The cache covers only the two full tutorials; guides still execute normally outside preview/draft. Caches live in ignored `docs/.execution-cache`, outside `docs/build`. Tiny checks are `julia --startup-file=no --project=docs docs/test/tutorial_cache.jl` and `julia --startup-file=no --project=docs docs/test/tutorial_render.jl`; neither runs real tutorials. Derive the expected figure inventory from the active tutorial and guide sources.

`npm --prefix docs run docs:dev` watches sources and regenerates draft output, so it omits executed plots. Avoid a concurrent watcher rewriting the same output: coordinate an owned watcher or use an isolated build. Use `npm --prefix docs run docs:preview` to inspect that result. Do not confuse previewing with publication.
