# fdid 0.2.0

* Initial public release.
* Implements six estimation methods: `ols1`, `ols2`, `did`, `ebal`, `ipw`, `aipw`.
* Supports robust, bootstrap, and jackknife variance estimation.
* Parallel computation via `foreach`/`doParallel` for bootstrap and jackknife.
* `fdid_prepare()` for reshaping long-format panel data.
* `fdid_list()` for comparing multiple estimates.
* S3 methods: `print`, `summary`, `plot` for `fdid` and `fdid_list` objects.
* Built-in `mortality` dataset (Cao et al. 2022, China Great Famine).
* Reference: Xu, Zhao & Ding (2026) <doi:10.1080/01621459.2026.2628343>.
