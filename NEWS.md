# ZKnock 1.0.4 (2024-09-28)

-   Added
    [generate_X](https://guang-yu-zhu.github.io/zKnock/reference/generate_X.html),
    [generate_y](https://guang-yu-zhu.github.io/zKnock/reference/generate_y.html),
    [generate_lp](https://guang-yu-zhu.github.io/zKnock/reference/generate_lp.html),
    and
    [generate_Weib](https://guang-yu-zhu.github.io/zKnock/reference/generate_Weib.html)
    functions to generate simulation data.
-   Added
    [create.seq](https://guang-yu-zhu.github.io/zKnock/reference/create.seq.html)
    and
    [create.sparse_seq](https://guang-yu-zhu.github.io/zKnock/reference/create.sparse_seq.html)
    functions to generate Sequential knockoffs.
-   Added
    [perf_eval](https://guang-yu-zhu.github.io/zKnock/reference/perf_eval.html)
    function to evaluate TDR and FDR.

# ZKnock 1.0.3 (2024-09-27)

-   Added
    [agg_Avg](https://guang-yu-zhu.github.io/zKnock/reference/agg_Avg.html),
    [agg_BH](https://guang-yu-zhu.github.io/zKnock/reference/agg_BH.html),
    and
    [agg_Freq](https://guang-yu-zhu.github.io/zKnock/reference/agg_Freq.html)
    functions to perform multiple knockoff procedures.
-   Revised
    [knockoff.filter](https://guang-yu-zhu.github.io/zKnock/reference/knockoff.filter.html)
    to accept aggregate functions as an argument.

# ZKnock 1.0.2 (2024-09-17)

-   Added
    [create.pls.knockoff](https://guang-yu-zhu.github.io/zKnock/reference/create.pls.knockoff.html)
    and
    [create.zpls.knockoff](https://guang-yu-zhu.github.io/zKnock/reference/create.zpls.knockoff.html)
    functions.

# ZKnock 1.0.1 (2024-08-15)

-   Added
    [create.knockoff](https://guang-yu-zhu.github.io/zKnock/reference/create.knockoff.html),
    a function to generate multiple knockoff samples.
-   Added
    [create.pc.knockoff](https://guang-yu-zhu.github.io/zKnock/reference/create.pc.knockoff.html),
    a function for generating knockoffs using a principal component
    construction algorithm.
-   Added
    [stat.SHAP](https://guang-yu-zhu.github.io/zKnock/reference/stat.SHAP.html)
    function.

# ZKnock 1.0.0 (2024-05-01)

-   [knockoff.filter](https://guang-yu-zhu.github.io/zKnock/reference/knockoff.filter.html)
    now accepts a family argument, which is passed to statistics
    functions like
    [stat.glmnet_coefdiff](https://guang-yu-zhu.github.io/zKnock/reference/stat.glmnet_coefdiff.html).
