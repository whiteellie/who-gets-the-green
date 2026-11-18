"""Checkpoint NLCD impervious area within developed pixels of modeled regional block groups."""
import argparse
import re
from pathlib import Path

import ee
import geopandas as gpd
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
ASSET = "USGS/NLCD_RELEASES/2021_REL/NLCD/2021"
AREA_COLUMNS = ["developed_area_m2", "valid_developed_area_m2", "impervious_area_m2"]
COLUMNS = ["GEOID", *AREA_COLUMNS, "impervious_fraction_developed", "potential_space_fraction", "impervious_coverage", "status"]


def build_image():
    nlcd = ee.Image(ASSET)
    landcover = nlcd.select("landcover")
    developed = landcover.gte(21).And(landcover.lte(24))
    impervious = nlcd.select("impervious")
    valid = developed.And(impervious.gte(0)).And(impervious.lte(100))
    area = ee.Image.pixelArea()
    return (
        area.updateMask(developed).rename("developed_area_m2")
        .addBands(area.updateMask(valid).rename("valid_developed_area_m2"))
        .addBands(area.multiply(impervious.divide(100)).updateMask(valid).rename("impervious_area_m2"))
    )


def load_checkpoint(path):
    if not path.exists():
        return pd.DataFrame(columns=COLUMNS)
    data = pd.read_csv(path, dtype={"GEOID": str})
    if not set(COLUMNS).issubset(data.columns) or data["GEOID"].duplicated().any():
        raise ValueError(f"Invalid checkpoint: {path}")
    return data.reindex(columns=COLUMNS)


def checkpoint(data, destination):
    temporary = destination.with_suffix(".tmp")
    data.to_csv(temporary, index=False)
    try:
        temporary.replace(destination)
    except PermissionError:
        print(f"Target locked; latest checkpoint remains in {temporary}", flush=True)


def summarize(features, requested_ids):
    properties = {str(feature["properties"]["GEOID"]): feature["properties"] for feature in features}
    rows = []
    for geoid in requested_ids:
        if geoid not in properties:
            raise ValueError(f"Earth Engine omitted requested GEOID {geoid}")
        row = {"GEOID": geoid, **{name: properties[geoid].get(name) for name in AREA_COLUMNS}}
        developed = row["developed_area_m2"] or 0
        valid = row["valid_developed_area_m2"] or 0
        impervious = row["impervious_area_m2"]
        row.update(impervious_fraction_developed=None, potential_space_fraction=None, impervious_coverage=None)
        if developed <= 0 or valid <= 0 or impervious is None:
            row["status"] = "no_valid_developed_pixels"
        else:
            fraction = impervious / valid
            coverage = valid / developed
            if not (0 <= fraction <= 1.000001 and 0 < coverage <= 1.000001):
                raise ValueError(f"Invalid NLCD area ratios for {geoid}: {fraction}, {coverage}")
            row.update(
                impervious_fraction_developed=min(1.0, fraction),
                potential_space_fraction=1 - min(1.0, fraction),
                impervious_coverage=min(1.0, coverage),
                status="ok",
            )
        rows.append(row)
    return pd.DataFrame(rows).reindex(columns=COLUMNS)


def main():
    parser = argparse.ArgumentParser(description="NLCD 2021 potential landscape-space proxy for retained-model block groups.")
    parser.add_argument("--region", required=True, help="Exact DWR hydrologic-region name.")
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    args = parser.parse_args()
    if args.chunk_size < 1 or (args.limit is not None and args.limit < 1):
        parser.error("--chunk-size and --limit must be positive")
    slug = "hr_" + re.sub(r"[^a-z0-9]+", "_", args.region.lower()).strip("_")
    regions = gpd.read_file(RAW / "hydrologic_regions.gpkg")
    region = regions.loc[regions["HR_NAME"] == args.region]
    if len(region) != 1:
        parser.error(f"Unknown region: {args.region}")
    folds = pd.read_csv(ROOT / "outputs" / "models" / "all_regions_fractional_logit_cv_predictions.csv")
    ids = set(folds.loc[folds["model"].str.startswith("M1:") & (folds["region"] == args.region), "pwsid"])
    if not ids:
        parser.error(f"No retained-model suppliers found for {args.region}")
    overlaps = pd.read_csv(PROCESSED / slug / "block_group_pws.csv", dtype={"GEOID": str}, usecols=["GEOID", "pwsid"])
    selected = set(overlaps.loc[overlaps["pwsid"].isin(ids), "GEOID"])
    bg = gpd.read_file(RAW / "acs_bg.gpkg")[["GEOID", "geometry"]]
    bg["GEOID"] = bg["GEOID"].astype(str)
    bg = gpd.clip(bg.loc[bg["GEOID"].isin(selected)].to_crs(4326), region.to_crs(4326))
    bg = bg.loc[~bg.geometry.is_empty].copy()
    if set(bg["GEOID"]) != selected:
        raise ValueError("Clipped selection does not match modeled block-group IDs.")
    bg["area_m2"] = bg.to_crs(3310).area
    bg = bg.sort_values(["area_m2", "GEOID"]).drop(columns="area_m2")
    destination = RAW / f"landscape_space_bg_{slug}.csv"
    temporary = destination.with_suffix(".tmp")
    data = load_checkpoint(destination)
    if temporary.exists():
        recovered = load_checkpoint(temporary)
        if not set(data["GEOID"]).issubset(set(recovered["GEOID"])):
            raise ValueError("Temporary checkpoint does not contain all committed rows.")
        if len(recovered) >= len(data):
            data = recovered
            if not args.dry_run:
                checkpoint(data, destination)
    pending = bg.loc[~bg["GEOID"].isin(data["GEOID"])]
    if args.limit is not None:
        pending = pending.head(args.limit)
    print(f"{args.region}: {len(selected)} modeled block groups; {len(pending)} pending", flush=True)
    if args.dry_run or pending.empty:
        return
    ee.Initialize(project="ca-water-data-challenge")
    ee.data.setDeadline(120_000)
    image = build_image()
    projection = ee.Image(ASSET).select("impervious").projection()
    for offset in range(0, len(pending), args.chunk_size):
        batch = pending.iloc[offset : offset + args.chunk_size]
        try:
            result = image.reduceRegions(
                collection=ee.FeatureCollection(batch.__geo_interface__),
                reducer=ee.Reducer.sum(),
                crs=projection,
                scale=30,
                tileScale=4,
            ).getInfo()
            rows = summarize(result["features"], batch["GEOID"].tolist())
        except ee.ee_exception.EEException as exc:
            if not args.continue_on_error:
                raise
            error_path = destination.with_name(destination.stem + "_errors.csv")
            errors = pd.read_csv(error_path, dtype={"GEOID": str}) if error_path.exists() else pd.DataFrame(columns=["GEOID", "error"])
            errors = pd.concat([errors, pd.DataFrame({"GEOID": batch["GEOID"].tolist(), "error": str(exc)})], ignore_index=True)
            errors.drop_duplicates("GEOID", keep="last").to_csv(error_path, index=False)
            print(f"Failed {len(batch)} features; left pending and logged to {error_path.name}", flush=True)
            continue
        data = rows.copy() if data.empty else pd.concat([data, rows], ignore_index=True)
        checkpoint(data, destination)
        print(f"Processed {min(offset + args.chunk_size, len(pending))}/{len(pending)}; checkpoint has {len(data)} rows", flush=True)
    print(f"Finished {args.region}; saved {len(data)} checkpoint rows", flush=True)


if __name__ == "__main__":
    main()