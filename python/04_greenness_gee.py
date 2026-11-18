"""Percent green (summer NDVI > threshold) per block group via Earth Engine, restricted to developed land (NLCD).

Python only because Earth Engine's Python API is the supported route. Needs data/raw/acs_bg.gpkg from R/03_download_acs.R.
"""
import argparse
import re
from pathlib import Path

import ee
import geopandas as gpd
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
GEE_PROJECT = "ca-water-data-challenge"
NDVI_YEAR = 2022
NDVI_GREEN_THRESHOLD = 0.3
SAMPLE_COUNTIES = ("06037", "06073", "06019", "06067", "06075", "06089")


def mask_clouds(image):
    classification = image.select("SCL")
    clear = classification.remap([4, 5, 6], [1, 1, 1], 0)
    return image.updateMask(clear)


def build_image(region):
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(region)
        .filterDate(f"{NDVI_YEAR}-06-01", f"{NDVI_YEAR}-10-01")
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
        .map(mask_clouds)
        .median()
    )
    green = s2.normalizedDifference(["B8", "B4"]).gt(NDVI_GREEN_THRESHOLD)
    nlcd = ee.Image("USGS/NLCD_RELEASES/2021_REL/NLCD/2021").select("landcover")
    developed = nlcd.gte(21).And(nlcd.lte(24))  # developed open -> high intensity
    return green.updateMask(developed).rename("green").addBands(developed.rename("developed"))


def main():
    parser = argparse.ArgumentParser(description="Default: one small block-group pilot, not a statewide run.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--full-state", action="store_true", help="Explicitly request statewide extraction; consumes substantially more quota.")
    mode.add_argument("--sample", action="store_true", help="Select one compact block group in each of six counties.")
    mode.add_argument("--county", help="Five-digit California county FIPS, e.g. 06067.")
    mode.add_argument("--region", help="DWR hydrologic region name, e.g. 'Sacramento River'.")
    parser.add_argument("--limit", type=int, help="Optional maximum pending block groups per run; county mode defaults to 10.")
    parser.add_argument("--chunk-size", type=int, default=25, help="Block groups per Earth Engine reduction request (default 25).")
    parser.add_argument("--continue-on-error", action="store_true", help="For regional runs, log failed requests and continue with later features.")
    parser.add_argument("--dry-run", action="store_true", help="Print selected areas without contacting Earth Engine.")
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.chunk_size < 1:
        parser.error("--chunk-size must be positive")
    if args.county and (len(args.county) != 5 or not args.county.isdigit() or not args.county.startswith("06")):
        parser.error("--county must be a five-digit California county FIPS")
    bg = gpd.read_file(RAW / "acs_bg.gpkg")[["GEOID", "geometry"]].to_crs(4326)
    if args.region:
        regions = gpd.read_file(RAW / "hydrologic_regions.gpkg")
        region = regions.loc[regions["HR_NAME"] == args.region]
        if len(region) != 1:
            parser.error(f"Unknown hydrologic region: {args.region}")
        bg = gpd.clip(bg, region.to_crs(bg.crs))
        bg = bg.loc[~bg.geometry.is_empty].sort_values("GEOID").reset_index(drop=True)
        if bg.empty:
            parser.error(f"No block groups intersect {args.region}")
    areas = bg.to_crs(3310).area
    if args.region:
        bg = bg.assign(_area_m2=areas).sort_values(["_area_m2", "GEOID"]).drop(columns="_area_m2").reset_index(drop=True)
        areas = bg.to_crs(3310).area
    if args.county:
        bg = bg.loc[bg["GEOID"].str.startswith(args.county) & areas.between(100_000, 2_000_000)].sort_values("GEOID")
        if bg.empty:
            parser.error("No compact block groups found in this county")
    elif args.sample:
        eligible = bg.loc[areas.between(100_000, 2_000_000)].copy()
        eligible["area_m2"] = areas.loc[eligible.index]
        selected = []
        for county in SAMPLE_COUNTIES:
            candidates = eligible.loc[eligible["GEOID"].str.startswith(county)].sort_values(["area_m2", "GEOID"])
            if candidates.empty:
                raise ValueError(f"No eligible block group in county {county}")
            selected.append(candidates.index[len(candidates) // 2])
        bg = bg.loc[selected]
    elif not args.full_state and not args.region:
        candidates = areas[areas >= 100_000]
        bg = bg.loc[[candidates.idxmin()]]
    label = f"COUNTY {args.county}" if args.county else f"REGION {args.region}" if args.region else "STATEWIDE" if args.full_state else "SAMPLE" if args.sample else "PILOT"
    filename = "greenness_bg.csv" if args.full_state else "greenness_bg_sample.csv" if args.sample else "greenness_bg_pilot.csv"
    if args.county:
        filename = f"greenness_bg_{args.county}.csv"
    if args.region:
        region_slug = re.sub(r"[^a-z0-9]+", "_", args.region.lower()).strip("_")
        filename = f"greenness_bg_hr_{region_slug}.csv"
    destination = RAW / filename
    columns = ["GEOID", "pct_green_of_developed", "developed"]
    df = pd.read_csv(destination, dtype={"GEOID": str}) if destination.exists() else pd.DataFrame(columns=columns)
    if not set(columns).issubset(df.columns) or df["GEOID"].duplicated().any():
        raise ValueError(f"Invalid checkpoint: {destination}")
    pending = bg.loc[~bg["GEOID"].isin(df["GEOID"])]
    if args.county:
        county_limit = args.limit if args.limit is not None else 10
        pending = pending.sample(n=min(county_limit, len(pending)), random_state=NDVI_YEAR).sort_values("GEOID")
    elif args.limit is not None:
        pending = pending.assign(_area_m2=areas.loc[pending.index]).sort_values(["_area_m2", "GEOID"]).head(args.limit).drop(columns="_area_m2")
    bg = pending
    plan = bg[["GEOID"]].copy()
    plan["area_km2"] = areas.loc[bg.index] / 1_000_000
    print(f"{label}: {len(bg)} block group(s)", flush=True)
    print(plan.head(20).to_string(index=False), flush=True)
    if len(plan) > 20:
        print(f"... and {len(plan) - 20} more block group(s)", flush=True)
    if args.dry_run:
        print("Dry run: no Earth Engine requests made.")
        return
    if pending.empty:
        print(f"All selected block groups already saved in {destination}")
        print(df.to_string(index=False))
        return
    ee.Initialize(project=GEE_PROJECT)
    ee.data.setDeadline(120_000)
    for offset in range(0, len(pending), args.chunk_size):
        features = pending.iloc[offset : offset + args.chunk_size]
        fc = ee.FeatureCollection(features.__geo_interface__)
        try:
            img = build_image(fc.geometry())
            res = img.reduceRegions(fc, ee.Reducer.mean(), scale=10, tileScale=4).getInfo()
        except ee.ee_exception.EEException as exc:
            if not args.continue_on_error:
                raise
            error_path = destination.with_name(f"{destination.stem}_errors.csv")
            errors = pd.read_csv(error_path, dtype={"GEOID": str}) if error_path.exists() else pd.DataFrame(columns=["GEOID", "error"])
            errors = pd.concat(
                [errors, pd.DataFrame({"GEOID": features["GEOID"].astype(str), "error": str(exc)})],
                ignore_index=True,
            ).drop_duplicates("GEOID", keep="last")
            errors.to_csv(error_path, index=False)
            print(f"Failed {', '.join(features['GEOID'].astype(str))}; logged to {error_path}", flush=True)
            continue
        batch = pd.DataFrame([feature["properties"] for feature in res["features"]]).rename(columns={"green": "pct_green_of_developed"}).reindex(columns=columns)
        df = pd.concat([df, batch], ignore_index=True)
        temporary = destination.with_suffix(".tmp")
        df.to_csv(temporary, index=False)
        temporary.replace(destination)
        print(f"Saved {min(offset + args.chunk_size, len(pending))} of {len(pending)} pending block groups", flush=True)
    print(df.to_string(index=False))
    print(f"Saved {destination}")


if __name__ == "__main__":
    main()
