---
name: quantumdevices-fit-avoided-crossing
description: Use fit_avoided_crossing to fit already sampled parameter-gap data without a Floquet calculation or model build.
---

Read [spectral reference](../../docs/src/reference/spectral_tools.md) and [implementation](../../src/utils/spectral_tools/avoided_crossing.jl).

`fit_avoided_crossing(parameters, gaps; fit_kwargs=(;))` fits `sqrt(minimum_gap^2 + slope^2*(parameter-center)^2)`. Supply matching real vectors with at least five distinct finite parameter values and finite nonnegative gaps. Samples are sorted together. The sampled minimum must be interior and data nonflat; the fitted center must remain bracketed. Inputs may describe any suitable isolated avoided crossing, not only Floquet spectra.

```julia
# A small synthetic isolated crossing; replace x and gaps with existing samples.
x = collect(range(-1.0, 1.0; length=9))
gaps = sqrt.(0.1^2 .+ 0.4^2 .* x.^2)
fit = fit_avoided_crossing(x, gaps)
fit.center, fit.minimum_gap, fit.rms_residual
```

`AvoidedCrossingFit` contains `center`, `minimum_gap`, `slope`, sorted `parameters`/`gaps`, `fitted_gaps`, `residuals` (fit minus data), `rms_residual`, normalization data, and `diagnostics` (convergence, Jacobian rank, singular values). Outputs retain input units; slope has gap units per parameter unit. Inspect residuals and interval dependence for competing crossings. Nonconvergence, nonfinite results, or rank deficiency throw errors. `fit_kwargs` forwards LsqFit options such as `maxIter`; `lower`, `upper`, and `inplace` are reserved. Do not resample or run dynamics when supplied gaps already answer the task.
