"""Extract LIDC-IDRI nodule annotations with pylidc.

Writes data/raw/annotations.csv with one row per radiologist annotation.
Annotations are grouped into physical nodules with Scan.cluster_annotations().
Everything is computed from the XML contours bundled with pylidc, so no DICOM
images are needed.

Run from the repository root:
    uv run python python/extract.py
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pylidc as pl
from scipy.spatial import QhullError
from skimage.morphology import convex_hull_image

OUT = Path(__file__).resolve().parent.parent / "data" / "raw" / "annotations.csv"

RATINGS = ["subtlety", "internalStructure", "calcification", "sphericity",
           "margin", "lobulation", "spiculation", "texture", "malignancy"]


def shape_features(ann, scan):
    """Extra geometric features from the annotation's binary mask."""
    mask = ann.boolean_mask()
    dz = scan.slice_spacing if scan.slice_spacing > 0 else scan.slice_thickness
    sp = np.array([scan.pixel_spacing, scan.pixel_spacing, dz])
    pts = np.argwhere(mask) * sp                     # voxel coordinates in mm
    out = {"n_voxels": int(mask.sum()),
           "n_slices": int(mask.shape[2]),
           "bbox_x": mask.shape[0] * sp[0],
           "bbox_y": mask.shape[1] * sp[1],
           "bbox_z": mask.shape[2] * sp[2],
           "elongation": np.nan,
           "flatness": np.nan,
           "convex_volume": np.nan,
           "solidity": np.nan}
    if len(pts) >= 4:
        ev = np.sort(np.linalg.eigvalsh(np.cov(pts.T)))[::-1]
        if ev[0] > 0:
            out["elongation"] = np.sqrt(max(ev[1], 0) / ev[0])
            out["flatness"] = np.sqrt(max(ev[2], 0) / ev[0])
    # Convex hull measured in voxels, so solidity is in (0, 1]. A mask one voxel
    # thick along some axis (e.g. a nodule drawn on a single slice) is flat, so
    # take the hull in the remaining axes. Other flat masks give an empty hull.
    flat = np.squeeze(mask)
    if flat.ndim >= 2:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")      # Qhull's degenerate-input noise
                hull = convex_hull_image(flat)
        except (QhullError, ValueError):
            hull = None
        if hull is not None and hull.any():
            hull |= flat                             # guarantee solidity <= 1
            out["convex_volume"] = hull.sum() * np.prod(sp)
            out["solidity"] = flat.sum() / hull.sum()
    return out


def main():
    scans = pl.query(pl.Scan).all()
    rows, skipped = [], []
    for i, scan in enumerate(scans, 1):
        try:
            clusters = scan.cluster_annotations(verbose=False)
        except Exception as e:                       # noqa: BLE001
            skipped.append((scan.patient_id, scan.id, repr(e)))
            continue
        for nod_idx, anns in enumerate(clusters):
            for reader_idx, ann in enumerate(anns):
                row = {"patient_id": scan.patient_id,
                       "scan_id": scan.id,
                       "nodule_id": f"{scan.id}_{nod_idx}",
                       "ann_id": ann.id,
                       "reader_slot": reader_idx,    # NOT a stable reader identity
                       "n_readers": len(anns),
                       "pixel_spacing": scan.pixel_spacing,
                       "slice_thickness": scan.slice_thickness,
                       "slice_spacing": scan.slice_spacing,
                       "diameter": ann.diameter,
                       "surface_area": ann.surface_area,
                       "volume": ann.volume}
                row.update({r: getattr(ann, r) for r in RATINGS})
                row.update(shape_features(ann, scan))
                rows.append(row)
        if i % 50 == 0:
            print(f"{i}/{len(scans)} scans, {len(rows)} annotations", flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)

    print(f"\nwrote {OUT} ({len(df)} rows, {df.shape[1]} columns)")
    print(f"patients {df.patient_id.nunique()}, scans {df.scan_id.nunique()}, "
          f"nodules {df.nodule_id.nunique()}")
    print("annotations per nodule:")
    print(df.groupby("nodule_id").size().value_counts().sort_index().to_string())
    for pid, sid, err in skipped:
        print(f"skipped scan {sid} ({pid}): {err}")


if __name__ == "__main__":
    main()
