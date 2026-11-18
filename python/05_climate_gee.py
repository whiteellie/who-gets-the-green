"""Extract 2022 growing-season GRIDMET precipitation and reference ET by block group."""
import argparse
import re
from pathlib import Path

import ee
import geopandas as gpd
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
GEE_PROJECT = "ca-water-data-challenge"
START_DATE = "2022-06-01"
END_DATE = "2022-10-01"
GRIDMET_SCALE_M = 4638.3
CLIMATE_COLUMNS = ["precip_gs_mm", "eto_gs_mm"]
COLUMNS = ["GEOID", *CLIMATE_COLUMNS]


def build_seasonal_image():
    daily = (
        ee.ImageCollection("IDAHO_EPSCOR/GRIDMET")
        .filterDate(START_DATE, END_DATE)
        .select(["pr", "eto"])
    )
    return daily.sum().rename(COLUMNS[1:])


def append_errors(path, failed_rows):
    if not failed_rows:
        return
    error_path = path.with_name(f"{path.stem}_errors.csv")
    errors = pd.read_csv(error_path, dtype={"GEOID": str}) if error_path.exists() else pd.DataFrame(columns=["GEOID", "error"])
    errors = pd.concat([errors, pd.DataFrame(failed_rows)], ignore_index=True).drop_duplicates("GEOID", keep="last")
    errors.to_csv(error_path, index=False)


def main():
    parser = argparse.ArgumentParser(description="Extract 2022 June-September GRIDMET climate for one hydrologic region.")
    parser.add_argument("--region", required=True, help="Exact DWR hydrologic-region name, e.g. 'Sacramento River'.")
    parser.add_argument("--limit", type=int, help="Optional maximum pending block groups for this run.")
    parser.add_argument("--chunk-size", type=int, default=100, help="Block groups per Earth Engine request (default 100).")
    parser.add_argument("--continue-on-error", action="store_true", help="Log failed requests and continue with later batches.")
    parser.add_argument("--dry-run", action="store_true", help="Print pending selection without contacting Earth Engine.")
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.chunk_size < 1:
        parser.error("--chunk-size must be positive")

    regions = gpd.read_file(RAW / "hydrologic_regions.gpkg")
    region = regions.loc[regions["HR_NAME"] == args.region]
    if len(region) != 1:
        parser.error(f"Unknown hydrologic region: {args.region}")

    bg = gpd.read_file(RAW / "acs_bg.gpkg")[["GEOID", "geometry"]].to_crs(4326)
    bg["GEOID"] = bg["GEOID"].astype(str)
    bg = gpd.clip(bg, region.to_crs(bg.crs))
    bg = bg.loc[~bg.geometry.is_empty].copy()
    bg["_area_m2"] = bg.to_crs(3310).area
    bg = bg.sort_values(["_area_m2", "GEOID"]).reset_index(drop=True)

    points = bg[["GEOID", "geometry"]].to_crs(3310)
    points.geometry = points.representative_point()
    points = points.to_crs(4326)

    slug = re.sub(r"[^a-z0-9]+", "_", args.region.lower()).strip("_")
    destination = RAW / f"climate_bg_hr_{slug}.csv"
    temporary = destination.with_suffix(".tmp")
    df = pd.read_csv(destination, dtype={"GEOID": str}) if destination.exists() else pd.DataFrame(columns=COLUMNS)
    if not set(COLUMNS).issubset(df.columns) or df["GEOID"].duplicated().any():
        raise ValueError(f"Invalid climate checkpoint: {destination}")
    if temporary.exists():
        recovered = pd.read_csv(temporary, dtype={"GEOID": str})
        if not set(COLUMNS).issubset(recovered.columns) or recovered["GEOID"].duplicated().any():
            raise ValueError(f"Invalid temporary climate checkpoint: {temporary}")
        if not set(df["GEOID"]).issubset(set(recovered["GEOID"])):
            raise ValueError(f"Temporary checkpoint does not contain all committed rows: {temporary}")
        if len(recovered) > len(df):
            df = recovered
            print(f"Recovered {len(df)} rows from {temporary}", flush=True)
            try:
                temporary.replace(destination)
            except PermissionError:
                print(f"Checkpoint target is locked; continuing from {temporary}", flush=True)

    pending = points.loc[~points["GEOID"].isin(df["GEOID"])]
    if args.limit is not None:
        pending = pending.head(args.limit)
    plan = bg.loc[pending.index, ["GEOID", "_area_m2"]].rename(columns={"_area_m2": "area_m2"})
    print(f"REGION {args.region}: {len(pending)} pending block group(s)", flush=True)
    print(plan.head(20).to_string(index=False), flush=True)
    if len(plan) > 20:
        print(f"... and {len(plan) - 20} more block group(s)", flush=True)
    if args.dry_run:
        print("Dry run: no Earth Engine requests made.")
        return
    if pending.empty:
        checkpoint = temporary if temporary.exists() else destination
        print(f"All selected block groups already saved in {checkpoint}")
        return

    ee.Initialize(project=GEE_PROJECT)
    ee.data.setDeadline(120_000)
    image = build_seasonal_image()
    for offset in range(0, len(pending), args.chunk_size):
        features = pending.iloc[offset : offset + args.chunk_size]
        collection = ee.FeatureCollection(features.__geo_interface__)
        try:
            result = image.reduceRegions(
                collection,
                ee.Reducer.first(),
                scale=GRIDMET_SCALE_M,
                tileScale=4,
            ).getInfo()
        except ee.ee_exception.EEException as exc:
            if not args.continue_on_error:
                raise
            append_errors(
                destination,
                [{"GEOID": geoid, "error": str(exc)} for geoid in features["GEOID"].astype(str)],
            )
            print(f"Failed {len(features)} block group(s); logged to {destination.stem}_errors.csv", flush=True)
            continue

        batch = pd.DataFrame([feature["properties"] for feature in result["features"]]).reindex(columns=COLUMNS)
        valid = batch[CLIMATE_COLUMNS].notna().all(axis=1)
        failed_rows = [
            {"GEOID": geoid, "error": "No GRIDMET values returned for representative point."}
            for geoid in batch.loc[~valid, "GEOID"].astype(str)
        ]
        returned = set(batch["GEOID"].astype(str))
        failed_rows.extend(
            {"GEOID": geoid, "error": "No feature returned by Earth Engine."}
            for geoid in features["GEOID"].astype(str)
            if geoid not in returned
        )
        append_errors(destination, failed_rows)
        batch = batch.loc[valid]
        if not batch.empty:
            df = pd.concat([df, batch], ignore_index=True)
            df.to_csv(temporary, index=False)
            try:
                temporary.replace(destination)
            except PermissionError:
                print(f"Checkpoint target is locked; latest rows remain in {temporary}", flush=True)
        print(f"Processed {min(offset + args.chunk_size, len(pending))} of {len(pending)} pending block groups", flush=True)

    checkpoint = temporary if temporary.exists() else destination
    print(f"Saved {len(df)} climate records to {checkpoint}")


if __name__ == "__main__":
    main()