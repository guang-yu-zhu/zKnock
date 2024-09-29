## ZKnock 1.0.0 (2024-05-01)
- `knockoff.filter` now accepts a `family` argument, which is passed to statistics functions like `stat.glmnet_coefdiff`.

## ZKnock 1.0.1 (2024-08-15)
- Added `create.knockoff`, a function to generate multiple knockoff samples.
- Added `create.pc.knockoff`, a function for generating knockoffs using a principal component construction algorithm.
- Added `stat.SHAP` function.

## ZKnock 1.0.2 (2024-09-17)
- Added `create.pls.knockoff` and `create.zpls.knockoff` functions.

## ZKnock 1.0.3 (2024-09-27)
- Added `agg_Avg`, `agg_BH`, and `agg_Freq` functions to perform multiple knockoff procedures.
- Revised `knockoff.filter` to accept aggregate functions as an argument.

## ZKnock 1.0.4 (2024-09-28)
- Added `generate_X`, `generate_y`, `generate_lp` and `generate_Weib` functions to generate simulation data.
- Added `create.seq` and `create.sparse_seq` functions to generate Sequential knockoffs
- Added `perf_eval` function to evaluate tdr and fdr.
