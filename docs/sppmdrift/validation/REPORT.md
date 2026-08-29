# M1 validation report

*Populated as cases are built. Each row states the metric, the tolerance, the
published value, the model value and the verdict — including honest fails.*

| Tier | Case | Metric | Tolerance | Published | Model | Verdict |
|------|------|--------|-----------|-----------|-------|---------|
| C | SG2000 Table 2 | max abs error, 15 entries | 0.15 | see PROVENANCE | 0.000 | **pass** |
| C | Settling vs Maggi (2013) | log10-RMSE, bb16 csf=0.7 | < 0.25 | — | 0.202 | **pass** |
| C | Settling vs Maggi (2013) | log10-RMSE, dietrich sphere | < 0.25 | — | 0.219 | **pass** |

Tier A, and Tier B (S1/S2/S3, W1/W2): not yet built.
