#!/usr/bin/env python3
"""Create a four-topounit Marcell Experimental Forest surface dataset.

The input template is expected to already be a single-gridcell ELM surface
dataset with either three or four topounits. When a three-topounit template is
used, the bog slot is duplicated to create a bog hollow / bog hummock pair.
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np


DEFAULT_TEMPLATE = Path("surfdata_marcell_3topounit.nc")
DEFAULT_OUTPUT = Path("surfdata_marcell_4topounit.nc")
DEFAULT_S2_OUTPUT = Path("surfdata_marcell_s2_4topounit.nc")
NTOP = 4
NAT_PFT_COUNT = 20
UPLAND_BOREAL_ENF_PFT = 2
PEATLAND_BOREAL_ENF_PFT = 3
UPLAND_DECID_NEEDLE_BOREAL_PFT = 4
PEATLAND_DECID_NEEDLE_BOREAL_PFT = 5
BROADLEAF_DECID_BOREAL_TREE_PFT = 10
UPLAND_BOREAL_SHRUB_PFT = 13
PEATLAND_BOREAL_SHRUB_PFT = 14
MOSS_PFT = 18
SEDGE_PFT = 19
PFT_AXIS_SOURCES_17_TO_20 = np.array(
    [0, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15, 16],
    dtype="i8",
)
NLEVSOI = 10
PEAT_ORGANIC_DENSITY = 130.0
MINERAL_ORGANIC_PROFILE = np.array([35, 25, 15, 8, 3, 1, 0, 0, 0, 0], dtype="f8")


def elm_soil_interfaces(nlevsoi: int = NLEVSOI, scalez: float = 0.025) -> np.ndarray:
    """Return ELM standard soil-layer interface depths for the surface file."""
    layer_midpoints = np.array(
        [scalez * (np.exp(0.5 * (j - 0.5)) - 1.0) for j in range(1, nlevsoi + 1)],
        dtype="f8",
    )
    interfaces = np.zeros(nlevsoi + 1, dtype="f8")
    for j in range(1, nlevsoi):
        interfaces[j] = 0.5 * (layer_midpoints[j - 1] + layer_midpoints[j])
    interfaces[nlevsoi] = layer_midpoints[-1] + 0.5 * (
        layer_midpoints[-1] - layer_midpoints[-2]
    )
    return interfaces


def organic_profile_from_peat_depth(peat_depth_m: float) -> np.ndarray:
    """Mix peat and mineral organic density by the fraction of each layer in peat."""
    if peat_depth_m <= 0.0:
        return MINERAL_ORGANIC_PROFILE.copy()

    interfaces = elm_soil_interfaces()
    organic = np.empty(NLEVSOI, dtype="f8")
    for lev in range(NLEVSOI):
        top = interfaces[lev]
        bottom = interfaces[lev + 1]
        peat_thickness = max(0.0, min(peat_depth_m, bottom) - top)
        peat_fraction = min(1.0, peat_thickness / (bottom - top))
        organic[lev] = (
            peat_fraction * PEAT_ORGANIC_DENSITY
            + (1.0 - peat_fraction) * MINERAL_ORGANIC_PROFILE[lev]
        )
    return organic


def set_if_present(ds: nc.Dataset, name: str, values) -> None:
    """Assign values to a NetCDF variable if the template contains it."""
    if name in ds.variables:
        ds.variables[name][:] = values


def set_topounit_flag(ds: nc.Dataset, name: str, values, long_name: str) -> None:
    """Assign a topounit/gridcell integer flag, creating it when needed."""
    if name not in ds.variables:
        var = ds.createVariable(name, "i4", ("topounit", "gridcell"))
        var.long_name = long_name
        var.units = "1"
        var.flag_values = np.array([0, 1], dtype="i4")
        var.flag_meanings = "false true"
    ds.variables[name][:] = values


def set_topounit_real(ds: nc.Dataset, name: str, values, long_name: str, units: str) -> None:
    """Assign a topounit/gridcell real variable, creating it when needed."""
    if name not in ds.variables:
        var = ds.createVariable(name, "f8", ("topounit", "gridcell"))
        var.long_name = long_name
        var.units = units
    ds.variables[name][:] = values


def set_topounit_int(ds: nc.Dataset, name: str, values, long_name: str, units: str) -> None:
    """Assign a topounit/gridcell integer variable, creating it when needed."""
    if name not in ds.variables:
        var = ds.createVariable(name, "i4", ("topounit", "gridcell"))
        var.long_name = long_name
        var.units = units
    ds.variables[name][:] = values


def remap_pft_axis(data: np.ndarray, axis: int, new_count: int) -> np.ndarray:
    """Remap a template PFT axis into the local upland/peatland split layout."""
    old_count = data.shape[axis]
    if old_count == new_count:
        return data
    if old_count > new_count:
        raise ValueError(f"Cannot shrink PFT axis from {old_count} to {new_count}")
    if old_count != 17 or new_count != NAT_PFT_COUNT:
        raise ValueError(f"Cannot remap PFT axis from {old_count} to {new_count}")

    new_shape = list(data.shape)
    new_shape[axis] = new_count
    if np.ma.isMaskedArray(data):
        fill_value = data.fill_value
        expanded = np.ma.masked_all(new_shape, dtype=data.dtype)
        expanded.set_fill_value(fill_value)
    else:
        expanded = np.zeros(new_shape, dtype=data.dtype)

    src_index = [slice(None)] * data.ndim
    dst_index = [slice(None)] * data.ndim
    for dst_pft, src_pft in enumerate(PFT_AXIS_SOURCES_17_TO_20):
        if src_pft >= old_count:
            raise ValueError(f"Cannot seed PFT {dst_pft} from missing source PFT {src_pft}")
        src_index[axis] = src_pft
        dst_index[axis] = dst_pft
        expanded[tuple(dst_index)] = data[tuple(src_index)]

    return expanded


def copy_template_with_topounits(template: Path, output: Path, ntopounits: int) -> None:
    """Copy a surfdata file, expanding topounit dimension when needed."""
    with nc.Dataset(template) as src, nc.Dataset(output, "w", format=src.data_model) as dst:
        if "topounit" not in src.dimensions:
            raise ValueError("Template must have a topounit dimension")

        old_ntop = len(src.dimensions["topounit"])
        if old_ntop == ntopounits:
            topounit_index = np.arange(ntopounits)
        elif old_ntop == 3 and ntopounits == 4:
            topounit_index = np.array([0, 1, 1, 2])
        else:
            raise ValueError(f"Cannot map {old_ntop} template topounits to {ntopounits}")

        for name, dim in src.dimensions.items():
            if name == "topounit":
                dst.createDimension(name, ntopounits)
            elif name in ("natpft", "lsmpft"):
                dst.createDimension(name, NAT_PFT_COUNT)
            else:
                dst.createDimension(name, None if dim.isunlimited() else len(dim))

        for attr in src.ncattrs():
            dst.setncattr(attr, src.getncattr(attr))

        for name, src_var in src.variables.items():
            kwargs = {}
            if "_FillValue" in src_var.ncattrs():
                kwargs["fill_value"] = src_var.getncattr("_FillValue")
            dst_var = dst.createVariable(name, src_var.datatype, src_var.dimensions, **kwargs)
            for attr in src_var.ncattrs():
                if attr != "_FillValue":
                    dst_var.setncattr(attr, src_var.getncattr(attr))

            data = src_var[...]
            if "topounit" in src_var.dimensions:
                axis = src_var.dimensions.index("topounit")
                data = np.take(data, topounit_index, axis=axis)
            for pft_dim in ("natpft", "lsmpft"):
                if pft_dim in src_var.dimensions:
                    axis = src_var.dimensions.index(pft_dim)
                    data = remap_pft_axis(data, axis, NAT_PFT_COUNT)
            dst_var[...] = data


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a Marcell forest upland/bog/fen topounit surface dataset."
    )
    parser.add_argument(
        "--template",
        type=Path,
        default=DEFAULT_TEMPLATE,
        help=f"Template surfdata file. Default: {DEFAULT_TEMPLATE}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output surfdata file. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--preset",
        choices=("generic", "s2"),
        default="generic",
        help=(
            "Surface-data preset. 'generic' preserves the original Marcell test "
            "fractions; 's2' uses S2 watershed fractions and bog seepage targets."
        ),
    )
    parser.add_argument(
        "--update-coordinates",
        action="store_true",
        help=(
            "Set LONGXY/LATIXY to approximate Marcell Experimental Forest center. "
            "By default, preserve the template coordinates for domain compatibility."
        ),
    )
    args = parser.parse_args()

    if args.preset == "s2" and args.output == DEFAULT_OUTPUT:
        args.output = DEFAULT_S2_OUTPUT

    if not args.template.exists():
        raise FileNotFoundError(args.template)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    copy_template_with_topounits(args.template, args.output, NTOP)

    # Topounit order follows the branch hydrology assumptions: index 1 is the
    # low/outlet unit, with adjacent units increasing in elevation.
    if args.preset == "s2":
        class_names = ("lagg_low_outlet", "bog_hollow", "bog_hummock", "upland_high")
        # Verry et al. S2 watershed partition: lagg 4%, bog 29%, upland 67%.
        # Split the bog evenly into hollow and hummock microtopounits.
        topounit_frac = np.array([[0.04], [0.145], [0.145], [0.67]], dtype="f8")
    else:
        class_names = ("fen_low_outlet", "bog_hollow", "bog_hummock", "upland_high")
        topounit_frac = np.array([[0.10], [0.10], [0.10], [0.70]], dtype="f8")
    if args.preset == "s2":
        # Approximate S2 elevations from Figure 7.5: low poor-fen lagg near
        # 420 m, bog dome/water-table contours near 422 m, and upland contours
        # mostly between 424 and 428 m. Use an area-representative upland
        # elevation rather than the highest ridge contour.
        topounit_elev = np.array([[420.2], [421.95], [422.10], [426.0]], dtype="f8")
        # Refined from Figure 7.5's 100 m map scale and S2 component areas:
        # bog is a compact dome draining radially to a narrow lagg, and the
        # upland centroid is much closer to the bog edge than the old 300 m
        # placeholder implied. The hummock-hollow distance remains a
        # microtopography length scale.
        lateral_dist = np.array([[0], [75], [1], [50]], dtype="i8")
        regional_target = np.array([[0], [1], [2], [1]], dtype="i4")
    else:
        topounit_elev = np.array([[420.0], [422.925], [423.075], [433.0]], dtype="f8")
        lateral_dist = np.array([[0], [150], [1], [300]], dtype="i8")
        regional_target = np.array([[0], [1], [2], [3]], dtype="i4")
    topounit_is_bog = np.array([[0], [1], [1], [0]], dtype="i4")
    peat_depth = np.array([[2.0], [3.0], [3.0], [0.0]], dtype="f8")
    if args.preset == "s2":
        # S2 bog deep seepage is about 1 cm/yr. Store ksat in mm/s.
        bog_till_ksat_mm_day = 10.0 / 365.0
    else:
        bog_till_ksat_mm_day = 0.1
    till_ksat = np.array(
        [[0.0], [bog_till_ksat_mm_day / 86400.0], [bog_till_ksat_mm_day / 86400.0], [0.0]],
        dtype="f8",
    )
    weighted_elev = np.sum(topounit_frac[:, 0] * topounit_elev[:, 0])

    # Natural vegetation PFT mixes. This assumes the local peatlands parameter
    # file where PFTs 2/4/13 are upland tree/shrub classes, PFTs 3/5/14 are
    # peatland-tuned versions, and PFTs 18/19 are moss/sedge.
    pft_mix = np.zeros((NAT_PFT_COUNT, NTOP, 1), dtype="f8")
    fen = {
        PEATLAND_BOREAL_ENF_PFT: 10.0,
        PEATLAND_DECID_NEEDLE_BOREAL_PFT: 10.0,
        BROADLEAF_DECID_BOREAL_TREE_PFT: 20.0,
        PEATLAND_BOREAL_SHRUB_PFT: 25.0,
        MOSS_PFT: 35.0,
    }
    bog_hollow = {
        PEATLAND_BOREAL_ENF_PFT: 25.0,
        PEATLAND_DECID_NEEDLE_BOREAL_PFT: 10.0,
        PEATLAND_BOREAL_SHRUB_PFT: 30.0,
        MOSS_PFT: 35.0,
    }
    bog_hummock = {
        PEATLAND_BOREAL_ENF_PFT: 45.0,
        PEATLAND_DECID_NEEDLE_BOREAL_PFT: 15.0,
        PEATLAND_BOREAL_SHRUB_PFT: 25.0,
        MOSS_PFT: 15.0,
    }
    if args.preset == "s2":
        # The S2 hydrology chapter names MEF upland aspen, red pine, and
        # black spruce forests, but does not give an S2 upland PFT percentage.
        # Use an even deciduous/evergreen boreal forest mix.
        upland = {UPLAND_BOREAL_ENF_PFT: 50.0, BROADLEAF_DECID_BOREAL_TREE_PFT: 50.0}
    else:
        upland = {
            UPLAND_BOREAL_ENF_PFT: 35.0,
            UPLAND_DECID_NEEDLE_BOREAL_PFT: 5.0,
            BROADLEAF_DECID_BOREAL_TREE_PFT: 45.0,
            UPLAND_BOREAL_SHRUB_PFT: 15.0,
        }
    for topounit_index, mix in enumerate((fen, bog_hollow, bog_hummock, upland)):
        for pft_index, pct in mix.items():
            pft_mix[pft_index, topounit_index, 0] = pct

    # Fen/lagg and bogs are Histosols down to TopounitPeatDepth; upland is a
    # sandy/loamy mineral soil. ORGANIC is kg/m3 organic matter, so boundary
    # layers are mixed by the fraction of each ELM soil layer above peat depth.
    # SOIL_ORDER uses ELM's soilorder_varcon indices: 3=Histosols, 10=Spodosols.
    organic = np.zeros((NLEVSOI, NTOP, 1), dtype="f8")
    for topounit_index in range(NTOP):
        organic[:, topounit_index, 0] = organic_profile_from_peat_depth(
            float(peat_depth[topounit_index, 0])
        )

    pct_sand = np.zeros((10, NTOP, 1), dtype="f8")
    pct_clay = np.zeros((10, NTOP, 1), dtype="f8")
    pct_sand[:, 0, 0] = 2.0
    pct_sand[:, 1, 0] = 2.0
    pct_sand[:, 2, 0] = 2.0
    pct_sand[:, 3, 0] = [65, 65, 60, 55, 50, 50, 45, 45, 40, 40]
    pct_clay[:, 0, 0] = 2.0
    pct_clay[:, 1, 0] = 2.0
    pct_clay[:, 2, 0] = 2.0
    pct_clay[:, 3, 0] = [8, 8, 10, 12, 12, 14, 14, 16, 16, 18]

    with nc.Dataset(args.output, "r+") as ds:
        if len(ds.dimensions["topounit"]) != NTOP:
            raise ValueError(f"Output must have exactly {NTOP} topounits")
        if len(ds.dimensions["gridcell"]) != 1:
            raise ValueError("Template must have exactly one gridcell")

        set_if_present(ds, "topoPerGrid", np.array([NTOP], dtype="i8"))
        set_if_present(ds, "TopounitFracArea", topounit_frac)
        set_if_present(ds, "TopounitAveElv", topounit_elev)
        set_if_present(ds, "TopounitLateralDist", lateral_dist)
        set_topounit_int(
            ds,
            "TopounitRegionalTarget",
            regional_target,
            "local target topounit for regional lateral groundwater exchange",
            "local topounit index; 0 means no target",
        )
        set_topounit_flag(ds, "TopounitIsBog", topounit_is_bog, "topounit bog flag")
        set_topounit_real(ds, "TopounitPeatDepth", peat_depth, "topounit peat depth", "m")
        set_topounit_real(
            ds,
            "TopounitTillKsat",
            till_ksat,
            "restrictive till saturated hydraulic conductivity below topounit peat",
            "mm s-1",
        )
        set_if_present(ds, "TOPO2", np.array([weighted_elev], dtype="f8"))
        set_if_present(ds, "MaxTopounitElv", np.array([np.max(topounit_elev)], dtype="f8"))
        set_if_present(ds, "TOPO", topounit_elev)

        set_if_present(ds, "PCT_NATVEG", np.full((NTOP, 1), 100.0))
        set_if_present(ds, "PCT_CROP", np.zeros((NTOP, 1)))
        set_if_present(ds, "PCT_WETLAND", np.zeros((NTOP, 1)))
        set_if_present(ds, "PCT_LAKE", np.zeros((NTOP, 1)))
        set_if_present(ds, "PCT_GLACIER", np.zeros((NTOP, 1)))
        set_if_present(ds, "PCT_NAT_PFT", pft_mix)

        if "PCT_URBAN" in ds.variables:
            ds.variables["PCT_URBAN"][:] = 0.0

        set_if_present(ds, "peatf", np.array([[1.0], [1.0], [1.0], [0.0]], dtype="f8"))
        set_if_present(ds, "PCT_SAND", pct_sand)
        set_if_present(ds, "PCT_CLAY", pct_clay)
        set_if_present(ds, "ORGANIC", organic)
        set_if_present(ds, "SOIL_COLOR", np.array([[16], [16], [16], [10]], dtype="i4"))
        set_if_present(ds, "SOIL_ORDER", np.array([[3], [3], [3], [10]], dtype="i4"))

        if "FMAX" in ds.variables:
            fmax = np.array(ds.variables["FMAX"][:], dtype="f8")
            fmax[0, 0] = 0.95
            fmax[1, 0] = 0.85
            fmax[2, 0] = 0.65
            ds.variables["FMAX"][:] = fmax
        set_if_present(ds, "SLOPE", np.array([[0.05], [0.10], [8.53], [2.00]], dtype="f8"))
        set_if_present(ds, "STD_ELEV", np.array([[0.25], [0.05], [0.05], [3.00]], dtype="f8"))
        set_if_present(ds, "F0", np.array([[0.20], [0.12], [0.08], [0.01]], dtype="f8"))

        if args.update_coordinates:
            set_if_present(ds, "LATIXY", np.array([47.50639], dtype="f8"))
            set_if_present(ds, "LONGXY", np.array([266.54444], dtype="f8"))

        if args.preset == "s2":
            site_description = "Marcell Experimental Forest S2 upland/bog/lagg surface dataset"
            fraction_attr = "lagg=0.04, bog_hollow=0.145, bog_hummock=0.145, upland=0.67"
            elevation_attr = "lagg=420.2, bog_hollow=421.95, bog_hummock=422.10, upland=426.0"
            peat_depth_attr = "lagg=2.0, bog_hollow=3.0, bog_hummock=3.0, upland=0.0"
            till_ksat_attr = (
                f"lagg=0.0, bog_hollow={bog_till_ksat_mm_day:.6g}, "
                f"bog_hummock={bog_till_ksat_mm_day:.6g}, upland=0.0"
            )
            pft_note = (
                "S2 chapter names MEF upland aspen, red pine, and black spruce forest types, "
                "but does not provide S2 percentages; upland uses 50% upland needleleaf "
                "evergreen boreal tree (PFT 2) and 50% broadleaf deciduous boreal tree. "
                "Wetland topounits keep the peatland-tuned boreal evergreen tree, "
                "deciduous needleleaf tree, and shrub PFTs 3, 5, and 14."
            )
            microtopography_note = (
                "bog_hummock is 0.15 m above bog_hollow with 1 m lateral separation; "
                "S2 elevations are estimated from Figure 7.5 contours with bog mean "
                "about 1.83 m above the lagg and upland about 3.98 m above the bog mean; "
                "S2 lateral distances are estimated from Figure 7.5's 100 m map scale "
                "and S2 component areas"
            )
            regional_target_attr = "lagg=0, bog_hollow=lagg, bog_hummock=bog_hollow, upland=lagg"
            lateral_dist_attr = "lagg=0, bog_hollow_to_lagg=75, bog_hummock_to_hollow=1, upland_to_lagg=50"
        else:
            site_description = "Representative Marcell Experimental Forest upland/bog/fen surface dataset"
            fraction_attr = "fen=0.10, bog_hollow=0.10, bog_hummock=0.10, upland=0.70"
            elevation_attr = "fen=420.0, bog_hollow=422.925, bog_hummock=423.075, upland=433.0"
            peat_depth_attr = "fen=2.0, bog_hollow=3.0, bog_hummock=3.0, upland=0.0"
            till_ksat_attr = "fen=0.0, bog_hollow=0.1, bog_hummock=0.1, upland=0.0"
            pft_note = (
                "Vegetation is limited to tree, shrub, and moss PFTs. Upland tree/shrub "
                "fractions use upland PFTs 2/4/13; wetland topounits use peatland PFTs 3/5/14."
            )
            microtopography_note = (
                "bog_hummock is 0.15 m above bog_hollow with 1 m lateral separation; "
                "mean bog elevation is 3 m above fen and 10 m below upland"
            )
            regional_target_attr = "fen=0, bog_hollow=fen, bog_hummock=bog_hollow, upland=bog_hummock"
            lateral_dist_attr = "fen=0, bog_hollow_to_fen=150, bog_hummock_to_hollow=1, upland_to_hummock=300"

        ds.setncattr("site_description", site_description)
        ds.setncattr("topounit_order", ", ".join(f"{i+1}={name}" for i, name in enumerate(class_names)))
        ds.setncattr("topounit_fractions", fraction_attr)
        ds.setncattr("topounit_elevations_m", elevation_attr)
        ds.setncattr("topounit_is_bog", "lagg/fen=0, bog_hollow=1, bog_hummock=1, upland=0")
        ds.setncattr("topounit_peat_depth_m", peat_depth_attr)
        ds.setncattr("topounit_till_ksat_mm_per_day", till_ksat_attr)
        ds.setncattr("topounit_regional_targets", regional_target_attr)
        ds.setncattr("topounit_lateral_dist_m", lateral_dist_attr)
        ds.setncattr("topounit_soils", "fen=organic Histosol, bogs=organic Histosol, upland=mineral Spodosol")
        ds.setncattr(
            "topounit_organic_note",
            (
                "ORGANIC is kg/m3 organic matter density. Layers above TopounitPeatDepth "
                f"use {PEAT_ORGANIC_DENSITY:g}; non-peat portions use the upland mineral "
                "organic profile, mixed by layer thickness at the peat boundary."
            ),
        )
        ds.setncattr("topounit_pft_note", pft_note)
        ds.setncattr(
            "topounit_drainage",
            (
                "regular water table generally uses default ELM drainage; lagg/fen and upland are "
                "non-bog topounits; bog perched water over till uses Darcy leakage to the regular "
                "aquifer and lateral exchange with adjacent active perched bog topounits"
            ),
        )
        ds.setncattr(
            "topounit_microtopography",
            microtopography_note,
        )
        ds.setncattr(
            "construction_note",
            (
                "Built from Marcell three-topounit surface data by splitting the bog topounit; PCT_WETLAND remains zero "
                "so bog/fen are represented as natural vegetation topounits with peatf=1."
            ),
        )
        old_history = getattr(ds, "history", "")
        stamp = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
        entry = f"{stamp}: created Marcell four-topounit {args.preset} variant from {args.template}"
        ds.setncattr("history", f"{old_history}\n{entry}" if old_history else entry)

    print(args.output)


if __name__ == "__main__":
    main()
