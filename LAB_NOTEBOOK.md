# California Water Data Challenge

## Purpose

This project explores how sampled greenness on developed land varies across public water system (PWS) service areas, alongside estimated outdoor residential water use, 2022 growing-season climate, developed-land fraction, and an ACS income proxy. It combines California water-system boundaries, supplier-reported monthly water use, block-group ACS estimates, Sentinel-2 imagery, and GRIDMET climate summaries processed in Google Earth Engine (GEE).

The figures and multivariable models are exploratory associations, not water-efficiency rankings or causal estimates. The adjusted models use 2022 growing-season climate and do not account for other years or all relevant confounders.

## Project layout

| Path | Purpose |
|---|---|
| `R/00_config.R` | Shared paths, years, source columns, and residential water categories. |
| `R/01`–`R/05` download scripts | PWS boundaries, monthly water use, ACS, hydrologic regions, and conservation reports. |
| `python/04_greenness_gee.py` | Sentinel-2 greenness and developed-fraction summaries. |
| `python/05_climate_gee.py` | 2022 GRIDMET precipitation and reference ET by regional block group. |
| `python/06_landscape_space_gee.py` | NLCD developed and impervious-area checkpoints. |
| `R/06_build_dataset.R` | Block-group/PWS intersections and supplier analysis tables. |
| `R/07_plot.R`, `R/08_plot_all_regions.R` | Regional figures and combined scatterplot. |
| `R/09_fit_models.R` | Pooled fractional-logit models and adjusted plots. |
| `R/10_fit_san_francisco_bay.R` | Single-region models; accepts any `hr_` slug. |
| `R/11`–`R/13` experiment scripts | Archived shortage-stage, housing-density, and landscape-space comparisons. |
| `R/14_make_infographic.R` | Builds the static statewide greenness and model-progression infographic. |
| `R/15_spatial_residual_analysis.R` | Maps final-model residuals by PWS-region and tests spatial autocorrelation. |
| `data/raw/`, `data/processed/` | Source/checkpoint data and analysis tables. |
| `outputs/models/`, `outputs/figures/` | Model tables, fitted objects, and figures; historical outputs are under `archived/`. |

Archived artifacts are kept under `outputs/figures/archived/` (conservation experiments, duplicate no-conservation figures, and the 06067 pilot), `data/processed/archived/06067/` (county pilot tables), and `data/raw/archived/` (the pilot and 06067 greenness checkpoints). Regional and current all-regions figures remain active. Historical pilot commands below expect the original active input paths; restore the archived inputs before replaying them. Running extraction without restoring its checkpoint can repeat already completed pilot work.

Run commands from the project root so `here::here()` resolves paths correctly.

## Setup and credentials

Install the R packages used by the scripts:

```r
install.packages(c(
  "here", "sf", "dplyr", "readr", "lubridate", "ggplot2",
  "httr2", "purrr", "tidycensus", "scales", "sandwich", "readxl", "ggrepel"
))
```

Install the Python packages:

```powershell
python -m pip install -r requirements.txt
```

The ACS download requires a Census API key. In R, configure it once with `tidycensus::census_api_key("YOUR_KEY", install = TRUE)`, then restart R. Do not commit the key.

The greenness extraction requires access to the GEE project configured as `ca-water-data-challenge` in `python/04_greenness_gee.py`, plus Earth Engine authentication. Authenticate in a terminal with:

```powershell
earthengine authenticate
```

The account must have access to that GEE project. If `Rscript` is not on `PATH` in Windows, run the installed `Rscript.exe` by its full path.

## Workflow, top to bottom

### 1. Configure shared settings

Review `R/00_config.R` before downloading or processing data. It sets the data directories, ACS year (2023), water-use year (2022), PWS and water-use column names, and residential categories. The 2022 water year is kept in sync with `NDVI_YEAR` in the Python script.

The scripts create the raw, processed, and figure directories automatically.

### 2. Download PWS boundaries

```powershell
Rscript R/01_download_pws_boundaries.R
```

The script pages through the ArcGIS feature service and writes `data/raw/pws_boundaries.gpkg`. The source has multiple boundary types; dataset construction keeps features whose type is `Water Service Area` and which have a PWS ID. A boundary-layer ID count is not necessarily a count of every regulated PWS in California.

### 3. Select and download monthly water use

```powershell
Rscript R/02_download_water_use.R
```

This first prints candidate data.ca.gov datasets and resource URLs. Inspect the candidates, set `water_use_url` in `R/02_download_water_use.R` to the chosen CSV resource URL, then run the script again. The download is intentionally skipped while `water_use_url` is `NULL`. The resulting file must contain the columns mapped in `R/00_config.R`; confirm those names and residential category values against the downloaded CSV before building the dataset.

### 4. Download ACS block-group data

```powershell
Rscript R/03_download_acs.R
```

This fetches California ACS 5-year block-group estimates for median household income (B19013), population (B01003), and housing units (B25001), and writes geometries and estimates to `data/raw/acs_bg.gpkg`. With ACS year 2023, the estimates cover the 2019–2023 period. This step requires the Census API key described above.

### 5. Download hydrologic-region boundaries

```powershell
Rscript R/04_download_hydrologic_regions.R
```

This retrieves the California Department of Water Resources' `i03 Hydrologic Regions` boundary layer and writes `data/raw/hydrologic_regions.gpkg`. The source has 12 polygon records but 10 distinct region names; disconnected pieces are dissolved by region name so the output has one multipolygon feature for each of the ten regions. Download these boundaries before running a region extraction.

### 6. Extract greenness in Earth Engine

The Python script sends block-group geometries to GEE; it does **not** download the satellite raster to the computer. GEE filters Sentinel-2 Surface Reflectance Harmonized imagery to June 1 through October 1, 2022, excludes scenes with cloudiness of 20% or more, masks pixels using Sentinel-2 scene classification, and forms a median composite. It computes NDVI from bands B8 and B4, then marks pixels with NDVI greater than 0.3 as green. NLCD 2021 classes 21–24 define developed land. **NLCD defines where to measure; Sentinel-2 NDVI defines what counts as green.**

`reduceRegions` applies a mean reducer at 10 m. The mean of the binary green band is the fraction of unmasked developed pixels classified as green; the mean of the developed band is the developed fraction of the block group. Only summary values are returned to the local CSV. The CSV is checkpointed as work completes so a later run can skip already-saved GEOIDs.

County mode supports small incremental runs:

```powershell
python python/04_greenness_gee.py --county 06067 --limit 10
```

This processes up to 10 pending, compact block groups in the county. Run it again to add another batch. For a complete hydrologic region, use its exact DWR name; geometries are clipped to the selected region and checkpointed to a region-specific CSV. For example:

```powershell
python python/04_greenness_gee.py --region "Sacramento River"
```

This region contains 2,317 intersecting ACS block groups. Region geometries are processed smallest-first, in batches of 25 by default, and checkpointed so reruns skip completed GEOIDs. Use `--chunk-size 5` for smaller requests if a large geometry times out, `--dry-run` to inspect the selection without contacting GEE, or `--limit 50` to process only the smallest 50 pending features. The default run is a one-block-group pilot; `--sample` selects one block group in each of six configured counties; `--full-state` explicitly requests statewide extraction and can consume substantially more GEE quota. The default pilot and `--sample` write separate CSV names and are not directly selected by the current R dataset builder.

### 7. Extract 2022 climate in Earth Engine

Run the climate extractor once for each hydrologic region, using its exact DWR name:

```powershell
python python/05_climate_gee.py --region "Sacramento River" --continue-on-error
```

The script uses the Earth Engine collection `IDAHO_EPSCOR/GRIDMET` and sums daily `pr` precipitation and `eto` grass-reference evapotranspiration from June 1 through October 1, 2022 (the end date is exclusive), producing seasonal totals in millimeters. GRIDMET is approximately 4 km resolution, so each region-clipped block group is assigned the values at its representative point. Results are checkpointed independently in `data/raw/climate_bg_hr_<slug>.csv`; failed requests are logged beside the CSV and can be retried. Use `--dry-run` to inspect pending block groups without contacting GEE.

NLCD 2021 developed fraction is already calculated by the greenness extraction: `developed` is the mean of a binary mask for classes 21–24. No separate NLCD download or extraction is needed for this covariate.

### 8. Build analysis datasets

For a county with greenness data, pass its five-digit California county FIPS:

```powershell
Rscript R/06_build_dataset.R 06067
```

For a hydrologic region, use the `hr_` prefix and lowercase underscore slug:

```powershell
Rscript R/06_build_dataset.R hr_sacramento_river
```

The script reads the matching greenness and climate CSVs, ACS block groups, PWS service areas, and monthly water-use records. Regional builds require the climate checkpoint from step 7. In region mode, block groups are clipped to the selected HR before intersecting with PWS service areas. Each overlap is weighted by its area fraction of the original block group times ACS population. This assumes population is distributed uniformly within each block group. Region outputs remain one row per PWS; PWSs crossing region boundaries contribute only their within-region spatial overlaps to that region's greenness, climate, developed-fraction, and income estimates. **Water-use values remain supplier-wide totals and are not apportioned to the HR.** A regional water-use aggregate would require an explicit allocation method for suppliers crossing region boundaries.

For water use, residential categories are summed by supplier and month. A monthly value is accepted only when both configured residential categories are present. Monthly gallons are converted to gallons per person per day. A supplier is marked complete only with all 12 valid months, all 3 winter baseline months (December–February), and positive finite annual use. Estimated outdoor use is annual daily use minus the winter daily-use baseline, floored at zero.

For each supplier, the builder calculates overlap-population-weighted means of sampled greenness, NLCD developed fraction, seasonal precipitation, reference ET, and block-group median household income. The latter is an **income proxy**, not a supplier-level median income. Regional supplier tables contain `developed_fraction_sampled`, `precip_gs_mm`, and `eto_gs_mm`; block-group/PWS tables retain the corresponding block-group values. Climate fields are populated for hydrologic-region builds. The builder also reports sampled population, sampled block groups, and the percentage of service-area area represented by the sampled overlaps.

The legacy county pilot tables are archived at `data/processed/archived/06067/pws_green_water_income.csv` and `data/processed/archived/06067/block_group_pws.csv`. Replaying the historical county build writes to the original active `data/processed/06067/` directory.
For Sacramento River HR, outputs are written to `data/processed/hr_sacramento_river/`.

For statewide processing, first create `data/raw/greenness_bg.csv` using the explicit `--full-state` extraction, then run:

```powershell
Rscript R/06_build_dataset.R statewide
```

The default scope is `statewide`, so omitting the argument also selects statewide inputs. Do not run that build until the statewide greenness CSV exists and has been checked.

### 9. Make figures

```powershell
Rscript R/07_plot.R 06067
```

For the regional figures:

```powershell
Rscript R/07_plot.R hr_sacramento_river
```

To exclude supplier records above a maximum outdoor GPCD from both scatterplots, pass a second argument. For the Sacramento River region, 132 GPCD is the rounded Tukey 1.5-IQR upper fence:

```powershell
Rscript R/07_plot.R hr_sacramento_river 132
```

For regional extraction, build, and screened figures, use the exact DWR region name and matching lowercase slug:

| Region | Greenness extraction | Notes |
|---|---|---|
| San Joaquin River | `python python/04_greenness_gee.py --region "San Joaquin River"` |  |
| Tulare Lake | `python python/04_greenness_gee.py --region "Tulare Lake"` |  |
| North Lahontan | `python python/04_greenness_gee.py --region "North Lahontan"` |  |
| South Lahontan | `python python/04_greenness_gee.py --region "South Lahontan" --chunk-size 5` | Smaller requests for large geometries. |
| Colorado River | `python python/04_greenness_gee.py --region "Colorado River" --chunk-size 5 --continue-on-error` | Retry failed requests as needed. |
| North Coast | Already extracted | Build and plot from the checkpoint. |
| Central Coast | `python python/04_greenness_gee.py --region "Central Coast" --continue-on-error` | Five large features timed out; retry with `--chunk-size 5`. |
| San Francisco Bay | `python python/04_greenness_gee.py --region "San Francisco Bay" --continue-on-error` |  |
| South Coast | `python python/04_greenness_gee.py --region "South Coast" --continue-on-error` |  |

For each row, follow with `Rscript R/06_build_dataset.R hr_<region_slug>` and `Rscript R/07_plot.R hr_<region_slug> 132` (replace the placeholder with the slug, e.g. `hr_san_joaquin_river`).

This is an outlier screen, not a physical definition of plausible use and not winsorization: values above the maximum are omitted from the scatterplots, not replaced with the maximum. The source and processed CSVs remain unchanged. To reduce label collisions, the scatterplots label only the ten highest retained GPCD suppliers; all supplier IDs remain in the CSV. Screened figures go to `outputs/figures/<scope>_gpcd_max_<maximum>/`. The sample-footprint map is regenerated without the GPCD filter because it shows where greenness was measured, not supplier water-use values.

The archived county pilot figures are in `outputs/figures/archived/06067/`; rerunning the pilot plot command recreates the active directory.

Use `Rscript R/07_plot.R statewide` only after statewide inputs and processed data exist. The plots filter to suppliers with complete water data and finite greenness and income values.

To combine all hydrologic-region supplier records into one scatterplot, run:

```powershell
Rscript R/08_plot_all_regions.R
```

The figure is saved to `outputs/figures/all_regions_gpcd_max_132/all_regions_green_vs_outdoor_use.png`. Points are colored by region; suppliers crossing regions may appear more than once. Water use remains supplier-wide, not apportioned, and 12 block-group/region records remain unprocessed.

### 10. Fit adjusted models

Fit the regional fractional-logit models and generate the result plots with:

```powershell
Rscript R/09_fit_models.R
```

The default run fits the retained no-conservation models on 1,282 supplier-region records at or below 132 GPCD, plus a no-GPCD-cutoff sensitivity fit on 1,365 proxy-covered records. Run the potential landscape-space extraction described below before fitting; `R/landscape_space_helpers.R` aggregates its checkpoints directly for the pooled and regional workflows, without depending on experimental model outputs. Main plots are in `outputs/figures/all_regions_model/` and model tables are in `outputs/models/`. `--baseline` is an alias for this retained analysis. Conservation trials are archived unsuccessful experiments and run only with explicit experiment flags; their historical M4 definitions are unchanged.

The revised "Who Gets the Green?" progression is:

1. M1, water use: does more residential outdoor use correspond to more green?
2. M2, + climate: how much do 2022 reference ET and precipitation explain?
3. M3, + development context: does developed-land fraction explain additional variation?
4. M4, + potential landscape space: does the non-impervious fraction of developed land explain additional variation?
5. M5, + income: does log income explain additional variation after landscape space? Conservation and housing density are not included.

The retained analysis fits M1-M5 without conservation on the same complete-case sample, including zero outdoor use. Outdoor use is log1p-transformed, income is log-transformed, and continuous predictors are standardized. Every pooled model includes region fixed effects and one common outdoor-use slope. Retained-model uncertainty and folds are grouped by PWS; original folds are preserved despite sensitivity-coverage exclusions. Archived conservation trials group systems sharing a reporting supplier and use their own common eligible samples. Developed fraction describes development context; potential landscape space is a separate imperviousness-based proxy, not actual irrigable or residential landscape area.

#### Failed Conservation Experiments

**Decision (2026-11-02): exclude conservation from the retained model.** Three exploratory encodings failed to improve predictive performance on their own matched samples and reporting-supplier folds:

- Reported-action count, 466 records: adding it to development-only M3 reduced CV R² from 0.509 to 0.500. Its M5 OR was 1.04 (95% CI 0.96–1.12).
- Assumed binary weekly restrictions, 430 records: adding the binary reduced CV R² from 0.505 to 0.498 before income and from 0.567 to 0.556 after income. Its M5 assumed Yes-versus-No OR was 0.95 (95% CI 0.76–1.18).
- Standardized shortage stage, 488 records: adding categorical peak stage reduced CV R² from 0.506 to 0.492 before income and from 0.578 to 0.564 after income. Sparse high stages limit precision.

These are unsuccessful predictive experiments with the available 2022 reporting data, not evidence that real conservation policies are ineffective. The count combines actions of unequal strength; the binary depends on an unverified unchecked-as-No assumption; standardized stage describes a shortage-plan response rather than actual watering restrictions or compliance. All restrict coverage. Data processing, assumptions, and outputs below are preserved for reproducibility, not used in the default fit or main plots.

Fetch and process the official [State Water Board conservation-actions resource](https://data.ca.gov/dataset/urws-conservation-supply-demand) with:

```powershell
Rscript R/05_download_conservation.R
```

The script downloads the per-column Excel workbook, documentation PDF, and companion per-supplier monthly supply/demand CSV to `data/raw/`. Downloads are cached; pass `--refresh` to retrieve them again. A source manifest records the URLs, file timestamps, and MD5 hashes in `data/processed/conservation/source_manifest.csv`. No API key is needed.

Processing preserves the original supplier reports, selects June-September 2022, expands the workbook's explicit comma-separated PWS lists, and retains the original reporting-supplier key. This is a reporting-unit membership link, not proof that each listed system had identical restrictions. One supplier-month report is required per PWS; conflicting duplicates stop processing. Reporting dates are monthly labels, not exact restriction effective dates.

Page 1 of the official documentation defines `Y` as a supplier reporting that it implemented an action in that month. Page 4 says the conservation-actions data dictionary is forthcoming. Consequently, `Weekly watering restrictions = Y` is coded as 1, and blanks remain unknown. A blank is not coded as 0, even when other conservation actions were reported. Local and workbook `DWR Stage` fields are preserved but not interpreted as interchangeable measures. The separate standardized-stage experiment below uses the PDF-defined `DWR_STANDARD_LEVEL` from the companion CSV instead. These are reports of actions, not measured compliance, watering-day limits, or restriction intensity.

Supplier summaries use 30, 31, 31, and 30 days for June-September (122 days). The seasonal `weekly_restriction_fraction` is calculated only when all four monthly values are known; otherwise it stays missing. The script also writes monthly coverage, positive-report days, unknown days, and regional match summaries. It never treats a missing report as no restrictions.

The verified snapshot contains 1,620 source reports from 409 reporting suppliers, listing 531 distinct PWS IDs after expansion. Of the 1,234 modeled suppliers, 480 match at least one report, 469 match all four reporting months, and 340 have four known weekly-restriction values. All four-known restriction histories equal 1, so that binary restriction variable still cannot be used to estimate a restriction association. Membership expansion improves coverage beyond the earlier direct-ID-only match of 339 suppliers.

For exploratory modeling, the processor instead counts checked `Y` responses across 47 explicitly selected action columns: demand reduction, use restrictions, communication/outreach, and enforcement actions. Supply augmentation, metadata, numerical incident counts, and reported types of water waste are excluded. The complete field list is saved in `counted_actions.csv`. Each selected checkbox contributes one, including communication channels and generic "Other" checkboxes; this is reporting breadth, not a standardized count of equally strong policies or residential-only measures.

Monthly counts are accepted when `Enacted any measures = Yes` and at least one selected action is checked, or when the response is explicitly `No` and no selected action is checked (count zero). Blank checkboxes do not become verified absences of individual measures. Empty Yes responses, missing reports, and No responses with positive checked counts remain unknown. The source-report audit identifies 1,464 accepted positive reports, 93 consistent explicit No reports, 57 conflicting No reports, and six empty/unknown Yes reports. These exclusions are recorded in `measure_response_audit.csv`.

The predictor `mean_reported_measures` is the day-weighted average monthly count over the 122-day window, calculated only when all four monthly counts are known. `model_readiness.csv` now evaluates this count and reports `model_ready = TRUE`; the binary weekly-restriction summary remains available for auditing but is not used in the model.

Reproduce the archived reported-action-count experiment explicitly:

```powershell
Rscript R/09_fit_models.R --conservation
```

This mode joins by PWS ID, requires four valid monthly counts and an unambiguous reporting-supplier key, and fits M1-M5 on the same eligible <=132 GPCD sample. M4 adds standardized `mean_reported_measures` to development-only M3; M5 then adds income. All five experiment steps use 466 supplier-region records from 446 PWSs and 350 reporting suppliers across ten regions. Systems sharing a report stay together for clustered uncertainty and CV folds. `sample_by_region.csv` records regional counts, including only seven North Lahontan records. Archived outputs go to `outputs/models/conservation_subset/` and `outputs/figures/archived/conservation_subset_model/`, never the main folders. Regional-only conservation fits are not enabled.

| Reported-action subset model | Records | CV R² | RMSE (percentage points) |
|---|---:|---:|---:|
| M1: water use + region | 466 | 0.364 | 16.72 |
| M2: add climate | 466 | 0.478 | 15.15 |
| M3: add development | 466 | 0.509 | 14.69 |
| M4: add reported actions | 466 | 0.500 | 14.83 |
| M5: add income | 466 | 0.561 | 13.89 |
| Sensitivity: M5, all GPCD | 478 | 0.539 | 14.18 |

Reported-action count does not improve CV performance when added to M3 in this fold split. In M5 its OR per SD is 1.04 (95% CI 0.96–1.12), an uncertain association. These scores cannot be directly compared with the larger baseline scores because sample membership and inference groups differ. The action count may reflect reporting practices, agency resources, drought severity, and mixtures of weak and strong measures; it is not evidence that restrictions did or did not work. The coefficient, prediction-comparison, calibration, and marginal-effect plots in the subset figure folder include this count.

Clarifying unchecked responses or obtaining independently verified restriction intensity is still required for a verified restriction-specific analysis. The separate sensitivity below explicitly tests an unverified unchecked-as-No assumption; it does not change the conservative source coding. Neither reported-action counts nor a later water-use-by-restriction interaction establish that drought policies caused a change.

The archived count experiment has six prediction panels. The main faceted predicted-versus-observed plot instead shows M1, M2, M3, M5 without conservation, and unrestricted sensitivity, with five panels. It is saved to `outputs/figures/all_regions_model/fractional_logit_models_predicted_vs_observed.png`; the main coefficient plot excludes conservation.

#### Failed Binary Restriction Experiment

```powershell
Rscript R/05_download_conservation.R
Rscript R/09_fit_models.R --binary-restrictions
```

This analysis assumes an unchecked `Weekly watering restrictions` field means No only in otherwise usable monthly reports (positive reported actions or a consistent explicit No). Checked Y means reported Yes. Missing/conflicting reports remain unknown. A seasonal binary is assigned only when all four months have the same assumed status; mixed and incomplete histories are excluded. This is not verified restriction-in-effect status, and reporting-month labels do not establish actual policy effective dates.

The <=132 GPCD sensitivity sample has 430 supplier-region records: 353 assumed Yes and 77 assumed No, representing 414 PWSs and 324 reporting suppliers across ten regions. All primary fits and the income-adjusted no-binary benchmark use these exact records and reporting-supplier-grouped folds. The binary predictor is not standardized: its OR compares assumed Yes with assumed No. Continuous predictors remain standardized.

| Binary sensitivity model | Records | CV R² | RMSE (percentage points) |
|---|---:|---:|---:|
| M1: water use + region | 430 | 0.419 | 15.92 |
| M2: add climate | 430 | 0.484 | 15.00 |
| M3: add development | 430 | 0.505 | 14.70 |
| M4: add assumed weekly restriction | 430 | 0.498 | 14.80 |
| M5: add income | 430 | 0.556 | 13.92 |
| Benchmark: income, no binary | 430 | 0.567 | 13.75 |
| Sensitivity: M5, all GPCD | 441 | 0.542 | 14.08 |

Adding the binary reduces CV R² by 0.007 before income and 0.011 after income in this fold split. The M5 assumed Yes-versus-No OR is 0.95 (95% CI 0.76–1.18), an uncertain association. Thus this coding does not improve measured predictive performance; it does not establish that actual restrictions are ineffective. The same-sample benchmark is the valid predictive comparator, not scores from the differently selected action-count or larger baseline analyses.

Archived assumption-labeled plots are saved under `outputs/figures/archived/binary_restriction_sensitivity_model/`, including model coefficients, predicted-versus-observed comparisons, and `binary_adjusted_predictions.png`. Model tables, group counts, assumptions, and incremental CV comparisons are under `outputs/models/binary_restriction_sensitivity/`. Coverage and excluded-status counts are in `data/processed/conservation/binary_restriction_coverage.csv`. This flag does not change the retained no-conservation results.

#### Standardized Shortage-Stage Experiment

```powershell
Rscript R/05_download_conservation.R
Rscript R/11_fit_shortage_stage.R
```

Page 4 of `data/raw/conservation_documentation.pdf` defines `DWR_STANDARD_LEVEL` as a standardized numeric level 0-6 corresponding to a supplier's local Water Shortage Contingency Plan stage; its accompanying flag identifies, among other cases, suppliers without a plan. The processor uses this documented field from the per-supplier monthly CSV, not a guessed mapping of free-text local stages or the workbook's `DWR Stage`. Blank, unrecognized, and nonblank-flagged stages remain unknown. Stage 0 is an explicit no-shortage-level-invoked category, not proof of no conservation activity.

The source has 1,592 June-September 2022 reports from 401 reporting suppliers. All four valid monthly stages and a single reporting organization are required per PWS. The predefined seasonal summary is the maximum reported stage during those four months, treated categorically with stage 0 as reference. No equal-spacing assumption or numeric severity slope is imposed. The eligible <=132 GPCD sample contains 488 supplier-region records from 465 PWSs and 365 reporting suppliers across ten regions. Peak-stage 0-5 counts are 34, 23, 381, 34, 12, and 4; stages 4 and 5 have only nine and four reporting suppliers. No eligible record reaches stage 6. This peak summary does not measure duration at a stage or exact policy implementation dates.

Five-fold CV is grouped by reporting organization (`ORG_ID`), stratified by peak stage, and checked to retain every observed stage and region in each training fold with a full-rank design. Fold assignments depend on predictor/category support only, not outcomes or predictive scores. M1-M5 and the no-stage income benchmark use identical observations and folds. Coefficient and marginal-prediction intervals are clustered by reporting organization.

| Standardized-stage experiment | Records | CV R² | RMSE (percentage points) |
|---|---:|---:|---:|
| M1: water use + region | 488 | 0.377 | 16.34 |
| M2: add climate | 488 | 0.480 | 14.93 |
| M3: add development | 488 | 0.506 | 14.54 |
| M4: add categorical peak stage | 488 | 0.492 | 14.75 |
| M5: add income | 488 | 0.564 | 13.68 |
| Benchmark: income, no stage | 488 | 0.578 | 13.45 |

Adding stage lowers CV R² by approximately 0.014 both before and after income in this fold split. A joint clustered test of the M5 stage coefficients has p = 0.042, but this exploratory association is not a predictive improvement or a causal policy finding. High categories are sparse, multiple exploratory encodings have been tested, and adjusted stage predictions may extrapolate across covariate combinations. The retained model therefore remains unchanged without conservation or stage.

Stage preprocessing and coverage audits are in `data/processed/conservation/stage_*.csv`. Model samples, folds, category support, coefficients, adjusted predictions, joint test, analysis definition, and paired CV comparisons are saved under `outputs/models/shortage_stage_experiment/`. The three experimental figures are under `outputs/figures/archived/shortage_stage_experiment/`: `predicted_vs_observed.png`, `stage_coefficients.png`, and `adjusted_greenness_by_stage.png`.

#### Retained Model Results

The default command (or its explicit `--baseline` alias) reproduces the retained larger-sample analysis without the reporting-coverage restriction:

```powershell
Rscript R/09_fit_models.R --baseline
```

Its coefficient, prediction, calibration, marginal-effect, and cumulative CV R² plots are saved under `outputs/figures/all_regions_model/`; model tables are under `outputs/models/`. `cumulative_cv_r_squared.png` shows R² in M1-to-M5 order, with each point labeled by its change from the previous step. Values are in `outputs/models/incremental_cv_r_squared.csv`. This comparison has five retained fits plus proxy-covered unrestricted sensitivity. These are the scores shown on the default main plots.

To test adding income immediately after the water-use-plus-region baseline, run `Rscript R/17_test_income_order.R`. This isolated experiment reuses the saved retained model inputs and exact PWS-grouped folds for all 1,282 records (1,234 unique suppliers), and verifies that baseline and final predictions reproduce the retained results. The alternative sequence is water use + region, income, climate, development, then potential landscape space. Income increases CV R² from 0.397 to 0.414 when added second (+0.0166), versus 0.720 to 0.724 when added last (+0.00387); RMSE reductions are 0.280 versus 0.0954 greenness percentage points. Income's earlier incremental gain is about 4.3 times larger, but the final model is unchanged (CV R² 0.724). Sequential gains depend on which predictors are already included and are not unique or causal measures of importance; this single-fold-assignment experiment does not establish statistical significance. Tables and predictions are under `outputs/models/income_order_experiment/`, and the comparison plot is under `outputs/figures/archived/income_order_experiment/`. Retained outputs and the infographic are not changed.

`final_model_observed_vs_predicted_anomalies.png` shows only final M5 out-of-fold predictions. Points with absolute residuals at least twice the pooled CV RMSE are highlighted by direction; the three largest positive and negative residuals are labeled by supplier name. This review threshold is not a formal outlier test. All record-level residuals and flags are in `outputs/models/all_regions_m5_cv_residuals.csv`.

Map residuals and test for geographic clustering with:

```powershell
Rscript R/15_spatial_residual_analysis.R
```

This clips PWS service areas to HRs and maps observed-minus-out-of-fold-predicted M5 greenness. Spatial neighbors are symmetrized 8-nearest representative points within each HR. The pooled test permutes residuals within HR; regional p-values are Benjamini–Hochberg adjusted. The residual map and Moran results are saved to `outputs/figures/all_regions_model/m5_residual_map.png` and `outputs/models/all_regions_m5_residual_morans_i.csv`. The infographic requires this script to run before `R/14_make_infographic.R`.

The stratified statewide Moran's I is 0.418 (9,999-permutation p = 0.0001). All ten regional I values are positive (0.084–0.677), and each is significant after BH adjustment (q ≤ 0.031). This suggests the final model leaves geographic structure unexplained; it is not evidence that location itself causes greenness. Neighbor choice is exploratory, and residuals describe sampled service-area overlaps rather than every parcel in the mapped boundary.

Screen the strongest positive or negative M5 residuals against CAL FIRE historical fire perimeters from 2015–2022 with:

```powershell
Rscript R/16_check_wildfire_exposure.R greener
Rscript R/16_check_wildfire_exposure.R less_green
```

The script searches within 10 km of each mapped PWS service area and reports direct perimeter intersections and unioned overlap area. The three strongest positive residual systems (Green Valley MWC, Valley View Park MWC, and Majestic Pines Community SD) had 4–15 perimeters within 10 km, but none directly intersected its mapped service area; their nearest perimeters were 2.92–5.22 km away. This does not support direct mapped-boundary wildfire overlap as a shared explanation for their positive residuals. The check uses full service-area boundaries, not only sampled block groups, and CAL FIRE's historical perimeter layer is incomplete and may be generalized; absence of a mapped overlap is not proof of no fire exposure. Direction-specific tables and perimeter GeoPackages are saved under `outputs/models/`.

| Model | Records | CV R² | RMSE (greenness percentage points) |
|---|---:|---:|---:|
| M1: outdoor use + region | 1,282 | 0.397 | 20.19 |
| M2: add climate | 1,282 | 0.509 | 18.23 |
| M3: add development | 1,282 | 0.524 | 17.94 |
| M4: add potential landscape space | 1,282 | 0.720 | 13.76 |
| M5: add income | 1,282 | 0.724 | 13.66 |
| Sensitivity: M5, proxy-covered, no GPCD cutoff | 1,365 | 0.721 | 13.65 |

R² is `1 - SSE/SST` on held-out predictions, not a likelihood-based pseudo-R². Retained out-of-fold predictions and metrics are saved under `outputs/models/`.

The default coefficient comparison and adjusted greenness-versus-outdoor-use plots remain in `outputs/figures/all_regions_model/` and exclude the conservation predictors. Their coefficients, sample counts, scaling parameters, adjusted predictions, and fitted models are saved under `outputs/models/`.

Build the infographic after extracting the ten regional greenness checkpoints and running the pooled model workflow:

```powershell
Rscript R/14_make_infographic.R
```

It writes `outputs/figures/california_greenness_infographic.png`, including the final-model out-of-fold residual screen and a separate key for the three largest residuals in each direction. These are exploratory flags, not supplier rankings or formal outlier tests.

For standalone PowerPoint figures, run `Rscript R/09_fit_models.R`, `Rscript R/15_spatial_residual_analysis.R`, and `Rscript R/14_make_infographic.R`. The scripts export four same-size (6.5 × 6.5 in), 300 dpi PNGs to `outputs/figures/powerpoint_panels/`:

| File | Source |
|---|---|
| `01_california_greenness_map.png` | `R/14_make_infographic.R` |
| `02_final_model_residual_map.png` | `R/15_spatial_residual_analysis.R` |
| `03_final_model_observed_vs_predicted.png` | `R/09_fit_models.R` |
| `04_cumulative_cv_r_squared.png` | `R/09_fit_models.R` |
| `05_standout_hydrologic_regions.png` | `R/15_spatial_residual_analysis.R` |

These are separate, title-free plot exports for arranging in slides. The first four are 6.5-inch square at 300 dpi and use 18 pt typography; the fifth is a portrait California map with 18 pt G1–G3/L1–L3 labels. The greenness map includes HR boundaries with its legend outside the map; the residual map uses gray only inside California. The scatter uses predicted greenness on x and observed greenness on y, with anomaly colors but no supplier labels. R² bars use green shades, M1–M5 predictor labels, and no “baseline” text. The residual choropleth remains unlabeled; use the fifth map for the standout IDs. The composite infographic retains its own panel titles and supplier standout key.

Two additional diagnostic figures are generated by the same script. In the main folder, `m5_adjusted_marginal_effects.png` shows retained M5 predictions across precipitation, reference ET, developed fraction, potential landscape space, and income, averaging over the other observed covariates and regions with PWS-clustered intervals. `fractional_logit_models_calibration.png` compares each prediction decile's mean prediction with mean observed greenness. These plots describe adjusted associations and calibration, not causal effects.

The retained M5 includes 1,282 supplier-region records from 1,234 PWSs. Potential landscape space has OR 2.65 per SD (95% CI 2.43-2.89), and log income has OR 1.13 per SD (95% CI 1.08-1.19). These are observational associations conditional on the other predictors, not causal effects or evidence of irrigation. The analysis also omits 12 region-feature greenness records listed in To Do.

The extraction targeted the original primary sample. Existing checkpoint coverage supports 1,365 of the former 1,404 unrestricted records; the remaining 39 are excluded, not imputed. Eligibility and population-weight coverage for all candidates are saved in `outputs/models/all_regions_landscape_space_coverage.csv`. Sensitivity scores therefore are not directly comparable to the old unrestricted score. The previous magnitude-comparison claim is not carried forward: `all_regions_m5_effect_contrast.csv` now explicitly describes a signed coefficient sum, not an absolute-magnitude comparison.

Before promotion, all 43 pooled and three-region model/figure files were copied and SHA-256 verified under `outputs/models/archived/pre_landscape_space_20261102/` and `outputs/figures/archived/pre_landscape_space_20261102/`. That snapshot preserves the previous income-without-space results and plots; current primary samples, folds, and M1-M3 predictions are unchanged.

#### Housing-Density Experiment

```powershell
Rscript R/12_fit_housing_density.R
```

This experiment uses the existing ACS B25001 housing-unit estimates and saved block-group/PWS overlaps; no new API download or Earth Engine extraction is needed. Housing units in each overlap are estimated as `housing_units * overlap_area_m2 / bg_area`, where `bg_area` is the original block-group area. Supplier-region density is the sum of those allocated counts divided by the summed overlap area in square kilometers, not an unweighted mean of block-group densities. Missing or invalid housing counts are not filled with zero; a supplier-region requires complete housing coverage.

The footprint matches the existing saved overlaps with positive population weights. It is gross polygon area, including any water or nonresidential land, not Census land-only area or developed/residential landscape area. Allocation assumes housing is uniformly distributed within each original block group. ACS estimates span 2019-2023 and are not a precise 2022 housing inventory; their estimation and allocation uncertainty is not propagated into the model intervals.

All 1,282 retained records (1,234 PWSs) have complete density values, so no records are lost. The four comparison models reuse the exact retained five-fold PWS assignments. These are historical pre-landscape-space formulas: the experiment's M5 means income without space. Density enters as `log1p(housing_units_per_km2)`; its coefficient is reported per SD of transformed density. Region, climate, developed fraction, and outdoor use remain in each fit, with paired comparisons before and after income.

| Housing-density comparison | Records | CV R² | RMSE (percentage points) |
|---|---:|---:|---:|
| M3: development, no income | 1,282 | 0.524 | 17.94 |
| M3 + housing density | 1,282 | 0.533 | 17.78 |
| M5: income, no density | 1,282 | 0.550 | 17.44 |
| M5 + housing density | 1,282 | 0.557 | 17.30 |

Density adds a modest CV R² gain of 0.008 before income and 0.007 after income in this fold split. Its adjusted M5 OR is 1.34 per SD (95% CI 1.16-1.55), a positive conditional association, not a causal effect. This is a more promising predictor than the unsuccessful conservation encodings, but the predictive gain is small and has not been validated on an independent dataset. Density remains excluded from the selected models; landscape space has been promoted separately.

Model data, overlap-allocation components, eligibility, coefficients, CV metrics, and adjusted predictions are saved under `outputs/models/housing_density_experiment/`. Comparison and adjusted-density plots are under `outputs/figures/archived/housing_density_experiment/`.

#### Potential Landscape-Space Experiment

Extract the proxy for each of the ten exact DWR region names, then fit the comparison models:

```powershell
python python/06_landscape_space_gee.py --region "Sacramento River" --continue-on-error
Rscript R/13_fit_landscape_space.R
```

Repeat the Python command for each region before the R comparison; reruns skip checkpointed GEOIDs. The extractor selects only block groups contributing to suppliers in the retained model, clips them to their hydrologic region, and writes separate `data/raw/landscape_space_bg_hr_<slug>.csv` checkpoints. It does not modify the existing greenness or climate checkpoints. `--dry-run` inspects pending selection, and `--limit` or `--chunk-size` controls request size.

The source is [NLCD 2021 in Earth Engine](https://developers.google.com/earth-engine/datasets/catalog/USGS_NLCD_RELEASES_2021_REL_NLCD), asset `USGS/NLCD_RELEASES/2021_REL/NLCD/2021`. `landcover` classes 21-24 define the same developed mask used for greenness, and `impervious` gives percent imperviousness at 30 m. Reductions use the native NLCD projection and 30 m resolution. The checkpoint stores developed area, developed area with valid impervious data, and fractional impervious area. The potential-space proxy is `1 - impervious_area_m2 / valid_developed_area_m2`; it is not the fraction of all regional land that is open.

All 24,451 selected region-block-group records are checkpointed. Four have no valid developed pixels at the reduction grid and remain missing, not zero imperviousness or 100% space. At aggregation, a block group requires at least 99% valid impervious-data coverage within its developed mask. Supplier-region values use the same overlap-population weights as greenness over greenness-contributing overlaps, and require proxy coverage of at least 99% of that population weight. Thus this is a population-weighted block-group proxy, not a direct area census of residential parcels. All 1,282 original model records remain eligible.

The six models use identical observations and the exact saved retained-model PWS folds. This historical experiment uses the pre-promotion numbering: its M5 is income without space, and its M5 + space is the current selected M5. Housing density is calculated as in the preceding experiment. Potential-space fraction enters linearly on the logit scale, with no equal-area or residential-only interpretation; all models retain the existing developed fraction, climate, outdoor use, and region controls, with income and density added in the specified comparisons.

| Landscape-space comparison | Records | CV R² | RMSE (percentage points) |
|---|---:|---:|---:|
| M3: development | 1,282 | 0.524 | 17.94 |
| M3 + landscape space | 1,282 | 0.720 | 13.76 |
| M5: income | 1,282 | 0.550 | 17.44 |
| M5 + landscape space | 1,282 | 0.724 | 13.66 |
| M5 + housing density | 1,282 | 0.557 | 17.30 |
| M5 + space + density | 1,282 | 0.738 | 13.30 |

This is a substantial predictive improvement in the matched-fold comparison. After income, space adds 0.174 to CV R²; after density and income, it adds 0.181. The space coefficient in the combined model has OR 2.72 per SD (95% CI 2.50-2.97). These are adjusted associations, not causal effects. Space accounts for a structural opportunity for vegetation within the same developed-land footprint; it is not proof of irrigation or a water-efficiency measure. Conditioning on space also changes the income/water-use question to comparisons at similar estimated landscape availability.

Non-impervious space may be bare, unirrigated, inaccessible, or nonresidential, and parks/golf courses remain included. NLCD is from 2021 while imagery is from 2022. The 30 m proxy and 10 m greenness measurement differ in resolution and valid-pixel support. Shared block-group measurements and spatial dependence across suppliers can make PWS-grouped validation optimistic for transfer to new places. Many predictors have been explored on these folds, so the result still needs independent validation rather than being treated as an unbiased final test score.

Raw components, population coverage, matched sample, coefficients, metrics, and prediction tables are under `outputs/models/landscape_space_experiment/`. The comparison and adjusted-space figures are under `outputs/figures/archived/landscape_space_experiment/`. Potential landscape space is now selected as M4 before income in the pooled and three standalone analyses. The historical experiment and pre-promotion snapshots remain available for comparison; density is not selected.

#### Standalone Regional Models

For a standalone San Francisco Bay analysis, run:

```powershell
Rscript R/10_fit_san_francisco_bay.R
```

This fits M1-M5 with five-fold PWS-grouped validation on the same 105 supplier records at or below 132 GPCD. CV R² increases from 0.314 with development-only M3 to 0.738 with M4 landscape space and 0.740 with M5 income (RMSE 11.93 greenness percentage points). Income adds little after space in this split. The no-cutoff sensitivity includes 108 proxy-covered records out of 109 otherwise complete records. Conservation and density are excluded. Treat this regional analysis cautiously given its smaller sample. Outputs are saved under `outputs/models/hr_san_francisco_bay/` and `outputs/figures/hr_san_francisco_bay_model/`.

The regional script accepts another `hr_` slug. For Sacramento River, run:

```powershell
Rscript R/10_fit_san_francisco_bay.R hr_sacramento_river
```

Sacramento River has 204 otherwise complete supplier records, with the same 185 primary records at or below 132 GPCD. CV R² rises from 0.160 at M3 to 0.509 with landscape space at M4, then falls slightly to 0.501 with income at M5 (RMSE 14.37 greenness percentage points). Income does not improve prediction in this fold split. The no-cutoff sensitivity includes 192 proxy-covered records. Conservation and density are excluded. Outputs are saved under `outputs/models/hr_sacramento_river/` and `outputs/figures/hr_sacramento_river_model/`.

For South Coast, which includes Los Angeles County, run:

```powershell
Rscript R/10_fit_san_francisco_bay.R hr_south_coast
```

South Coast has 318 otherwise complete supplier records, with the same 301 primary records at or below 132 GPCD. CV R² rises from 0.247 at M3 to 0.637 with landscape space at M4 and 0.670 with income at M5 (RMSE 9.81 greenness percentage points). The no-cutoff sensitivity includes 315 proxy-covered records. Conservation and density are excluded. Outputs are saved under `outputs/models/hr_south_coast/` and `outputs/figures/hr_south_coast_model/`.

## Regional and archived pilot status

The archived Sacramento County (`06067`) pilot has 30 sampled block groups; its inputs, tables, and figures are in the archive locations listed above. The Sacramento River HR run has 2,317 intersecting block groups and matched 759 PWSs; 210 have complete monthly water records. Ten block groups have no valid greenness observation. The 132 GPCD plot screen retained 185 of 204 suppliers otherwise eligible for the scatterplots.

The San Joaquin River HR run has 1,524 intersecting block groups and matched 925 PWSs; 161 have complete monthly water records. Ten block groups have no valid greenness observation. Using the same 132 GPCD screen for comparability, 129 of 152 suppliers otherwise eligible for the scatterplots were retained. This threshold was derived from the Sacramento River distribution, not independently calibrated for San Joaquin River. Review coverage and missingness before interpreting regional comparisons.

The Tulare Lake HR run has 1,574 intersecting block groups and matched 658 PWSs; 125 have complete monthly water records. Six block groups have no valid greenness observation. Using the shared 132 GPCD comparison screen, 110 of 123 suppliers otherwise eligible for the scatterplots were retained. The screen was derived from the Sacramento River distribution, not independently calibrated for Tulare Lake.

The North Lahontan HR run has 122 intersecting block groups and matched 124 PWSs; 22 have complete monthly water records. Four block groups have no valid greenness observation. Using the shared 132 GPCD screen, 18 of 20 suppliers otherwise eligible for the scatterplots were retained. This is a small supplier sample; interpret regional patterns cautiously.

The South Lahontan extraction has checkpointed 604 of 606 intersecting block groups. Two very large block groups (`060710103001` and `060270008001`) remain unprocessed after Earth Engine computation timeouts; six of the saved 604 have no valid greenness. The partial build matched 282 PWSs, with 98 complete monthly water records. The 132 GPCD screen retained 82 of 96 suppliers otherwise eligible for the plots. These tables and figures do not represent the two unprocessed block groups.

The Colorado River extraction checkpoint has 585 of 590 intersecting block groups. Five very large block groups remain unprocessed after request/checkpoint failures; two of the saved 585 have no valid greenness. The partial build matched 131 PWSs, with 60 complete monthly water records. The shared 132 GPCD screen retained 58 of 60 suppliers otherwise eligible for plots. These results do not include the five unprocessed block groups.

The North Coast greenness checkpoint has all 608 intersecting block groups, with one missing greenness result and no failed-feature log. The build matched 672 PWSs, 116 with complete monthly water records. The shared 132 GPCD screen retained 107 of 114 suppliers otherwise eligible for the plots. No further Earth Engine extraction is needed for this region.

The Central Coast checkpoint has 1,036 of 1,041 intersecting block groups; five large features timed out, and one saved block group has no valid greenness. The partial build matched 682 PWSs, with 217 complete monthly water records. Of 208 suppliers otherwise eligible for scatterplots, 21 exceeded the shared 132 GPCD screen and 187 were plotted. The three figures are in `outputs/figures/hr_central_coast_gpcd_max_132/`.

The San Francisco Bay checkpoint has all 4,640 intersecting block groups, with three missing greenness results and no failed-feature log. The build matched 361 PWSs, with 109 complete monthly water records. Of 109 suppliers otherwise eligible for scatterplots, four exceeded the shared 132 GPCD screen and 105 were plotted. The three figures are in `outputs/figures/hr_san_francisco_bay_gpcd_max_132/`.

The South Coast checkpoint has all 12,945 intersecting block groups, with one missing greenness result and no failed-feature log. The build matched 492 PWSs, with 322 complete monthly water records. Of 318 suppliers otherwise eligible for scatterplots, 17 exceeded the shared 132 GPCD screen and 301 were plotted. The three figures are in `outputs/figures/hr_south_coast_gpcd_max_132/`.

## Interpretation and limitations

- Greenness is a thresholded NDVI proxy on NLCD-developed pixels, not a direct measure of tree canopy, irrigation, lawn area, or landscape quality.
- The greenness outcome covers all NLCD developed classes (21–24), not residential land only; parks, golf courses, and other commercial/institutional landscapes can contribute. The outdoor-water-use predictor includes single- and multi-family residential use only, so the two measures have mismatched land-use coverage and greenness cannot be attributed solely to residential watering. Isolating residential greenness would require a residential parcel or land-use mask and recalculating the outcome.
- The satellite composite is from 2022 while the developed-land mask is NLCD 2021 and ACS estimates span 2019–2023.
- Supplier-wide water-use estimates are compared with greenness and income summarized only over sampled block-group/service-area overlaps.
- The winter-use baseline is a rough indoor-use proxy; the water data are supplier-reported, and supplier coverage/completeness varies.
- Climate adjustment uses only 2022 growing-season precipitation and reference ET; it does not represent long-term climate or other years.
- The reported-action count is exploratory: actions and communication channels receive equal weight, measures can target nonresidential users, and reporting coverage is incomplete. It measures reported activity, not restriction intensity, enforcement effectiveness, or compliance. Restricted-sample model performance is not directly comparable to the larger baseline sample.
- Binary restriction sensitivity assumes unchecked responses mean No in otherwise usable reports and excludes mixed histories. This assumption is unverified, can misclassify exposure, and must not be presented as confirmed policy-in-effect status or a causal policy evaluation.
- Twelve region-feature greenness records remain unprocessed, and the analysis should not be used to rank suppliers by efficiency.

## To Do

Conservation is excluded from the retained model; all three unsuccessful predictive trials are recorded above. The standardized-stage trial uses a documented field but still provides no held-out improvement. Reconsider these predictors with better temporal alignment, verified restriction exposure, or additional years rather than searching for a favorable R². Preserve PWS/reporting-unit identifiers and source citations; exact effective/end dates and residential watering-day limits would improve exposure measurement. Missing restrictions remain unknown outside the explicitly labeled binary sensitivity. One-year observational data cannot establish that drought restrictions caused a change.

Retry the remaining Earth Engine block groups. There are 12 region-feature records across three regions (11 unique GEOIDs; `060710103001` intersects both South Lahontan and Colorado River):

- Central Coast: `060790103012`, `060530114002`, `060690008021`, `060530114001`, `060790127081` (timed out; retry with `--chunk-size 5 --continue-on-error`).
- South Lahontan: `060710103001`, `060270008001` (timed out; retry with `--chunk-size 5`).
- Colorado River: `060710104252`, `060710103001`, `060710251001`, `060250124001`, `060650469003` (request/checkpoint failures; retry with `--chunk-size 5 --continue-on-error`).