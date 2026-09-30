# 14b_albuminuria_models
albuminuria_models.tsv: rank-OLS (HC3) of albuminuria/proteinuria band rank on each iPT score + diagnosis, then + eGFR band rank, + cortical IF% rank, + RAAS, + all; both orientations; standardized (z) coefficients, n per complete-case sample; estimate_raw = 12c raw-unit score ~ diagnosis + pct rank.
slope_difference.tsv: marker (DKDiPT_up) minus common injury (SharedInjury_up) standardized coefficient; separate models with 2000 diagnosis-stratified participant bootstrap resamples (seeded), and joint model with HC3 Wald test.
Scopes: all = DKD+HKD+HKD_withDM; DKD_only (no diagnosis term); no_HKD_withDM = DKD+HKD. sn primary, sc secondary.
Bands are ordinal percentile ranks (no midpoints); p values are exploratory (analysis_plan change 4/5 family).
