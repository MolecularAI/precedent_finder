"""Module containing a CLI tool for the Precent Finder"""

import argparse
from pathlib import Path
from typing import Tuple

import pandas as pd
import numpy as np
from paretoset import paretoset
from scipy import sparse
from scipy.spatial.distance import cdist
from rxnutils.chem.reaction import ChemicalReaction
from rxnutils.chem.features.reaction_centre_fp import (
    ReactionCentreFingerprint,
    reaction_center_similarity,
)
from rxnutils.chem.features.simple_rdkit import SimpleRdkitFingerprint
from rxnutils.chem.smartslib import SmartsLibrary


def _calc_center_similarity(
    reaction: ChemicalReaction,
    precomp_rcfp: np.ndarray,
    numbits: int = 1024,
    max_centers: int = 6,
):
    center_fps = ReactionCentreFingerprint(numbits=numbits)(reaction)[:max_centers]
    similarities = np.zeros(len(precomp_rcfp))
    for idx, fp_spec in enumerate(precomp_rcfp):
        # Un-pack the concatenated fingerprint
        rcfp2 = [
            fp_spec[1 + (numbits * idx) : 1 + (numbits * (idx + 1))]
            for idx in range(fp_spec[0])
        ]
        similarities[idx] = reaction_center_similarity(center_fps, rcfp2)
    return similarities


def _calc_global_similarity(
    reaction: ChemicalReaction,
    precomp_rfp: np.ndarray,
    numbits: int = 512,
):
    global_fp = SimpleRdkitFingerprint("fingerprint_mixed", numbits=numbits)(reaction)
    global_fp = np.asarray(global_fp)
    global_similarity = cdist(precomp_rfp, global_fp.reshape((1, -1)), "jaccard")
    return 1 - global_similarity.flatten()


def _calc_pareto(
    rxn_similarity: np.ndarray,
    center_similarity: np.ndarray,
    yields: np.ndarray,
    dates: np.ndarray,
    relative_cutoff: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray]:
    data = np.column_stack([rxn_similarity, center_similarity, yields, dates])

    sim_max = data[:, :2].max(axis=0)
    distances = cdist(data[:, :2], [sim_max], "sqeuclidean")
    sel = distances.flatten() < relative_cutoff
    if sel.sum() == 0:
        sel = distances.flatten() < min(1.5 * relative_cutoff, 1.0)
    if sel.sum() == 0:
        return [], []
    true_indices = np.where(sel)[0]
    data = data[sel, :]

    mask = paretoset(data, sense=["max"] * 4, distinct=False, use_numba=True)
    indices = np.where(mask)[0]
    subdata = data[indices]
    refpoint = [1.0, 1.0, 1.0, 1.0]
    weights = [1.0, 1.0, 0.8, 0.2]
    factor = np.sqrt(np.sum(np.power(weights, 2)))
    distances = cdist(subdata, [refpoint], "euclidean", w=weights).flatten()
    goodness = (factor - distances) / factor

    sort_indices = np.argsort(goodness)[::-1]
    return true_indices[indices[sort_indices]], goodness[sort_indices]


def find_precedents(
    reaction_smiles: str,
    metadata: pd.DataFrame,
    rfp_precomp: sparse.csr,
    rcfp_precomp: sparse.csr,
    smartslib_path: str,
    numbits: int = 1024,
    max_centers: int = 6,
) -> pd.DataFrame:
    """
    Find precedents of a query reaction using
    1. Sub-selection of data based on reactive function tag.
    2. Pareto-optimization in four dimensions (global reaction similarity, reaction center similarity, yield and date)

    The reaction SMILES need to contain atom-mapping

    :param reaction_smiles: the SMILES of the query reaction
    :param metadata: the processed metadata that includes columns ReactiveFunctions, CuratedYield and Year
    :param rfp_precomp: the pre-computed reaction fingerprints
    :param rcfp_precomp: the pre-computed reaction center fingerprints
    :param smartslib_path: the path to the SMARTS library
    :param numbits: number of bits for the fingerprint
    :param max_centers: the maximum number of reaction centers to consider
    :returns: a dataframe of the precedents, ranked by goodness
    """
    rxn = ChemicalReaction(reaction_smiles, clean_smiles=False)

    lib = SmartsLibrary.load(smartslib_path)
    function_keys, _ = lib.detect_reactive_functions(
        rxn,
        sort=True,
        add_none=True,
        target_size=None,
        max_reactants=5,
    )
    function_key = "|".join(function_keys) 

    sel = metadata["ReactiveFunction"] == function_key
    idx = metadata[sel].index
    sel_metadata = metadata[sel]
    sel_rfp = rfp_precomp[idx, :].toarray()
    sel_rcfp = rcfp_precomp[idx, :].toarray()

    global_similarities = _calc_global_similarity(rxn, sel_rfp, numbits)
    center_similarities = _calc_center_similarity(rxn, sel_rcfp, numbits, max_centers)
    yields = metadata.loc[idx, "CuratedYield"].values / 100.0
    dates = (
        pd.to_datetime(metadata.loc[idx, "Year"])
        .astype(int)
        .values
    )
    dates -= dates.min()
    dates = dates / dates.max()
    top_indices, goodness = _calc_pareto(
        global_similarities, center_similarities, yields, dates
    )

    results = sel_metadata.iloc[top_indices].to_dict(orient="list")
    results["global_similarity"] = global_similarities[top_indices].tolist()
    results["center_similarity"] = center_similarities[top_indices].tolist()
    results["goodness"] = goodness.tolist()
    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--smiles",
        help="the atom-mapped reaction SMILES query. Leaving it blank will start and interactive session",
    )
    parser.add_argument(
        "--output", default="precedents.csv", help="the output of the analysis"
    )
    parser.add_argument(
        "--data_folder",
        default=str(Path(__file__).parent / "data"),
        help="the folder with pre-processed USPTO data",
    )
    parser.add_argument("--fp_len", type=int, default=1024)
    parser.add_argument("--max_centers", type=int, default=6)
    args = parser.parse_args()

    rfp = sparse.load_npz(Path(args.data_folder) / "rfp_precomp.npz")
    rcfp = sparse.load_npz(Path(args.data_folder) / "rcfp_precomp.npz")
    data = pd.read_csv(Path(args.data_folder) / "metadata_cleaned.csv.gz", sep="\t")
    smartslib_path = str(Path(args.data_folder) / "group_smarts.json")

    if args.smiles:
        results = find_precedents(
            args.smiles, data, rfp, rcfp, smartslib_path, args.fp_len, args.max_centers
        )
        results.to_csv(args.output, sep="\t")
        print(f"Stored {len(results)} precedents to {args.output}")
        return

    n_pred = 1
    while True:
        print("Starting interactive session. Leave SMILES blank to stop.")
        try:
            rsmi = input("Enter a reaction SMILES: ")
        except (EOFError, KeyboardInterrupt):
            break
        if not rsmi:
            break
        results = find_precedents(
            rsmi, data, rfp, rcfp, smartslib_path, args.fp_len, args.max_centers
        )
        filename = args.output.replace(".csv", f"{n_pred}.csv")
        results.to_csv(filename, sep="\t")
        print(f"Stored {len(results)} precedents to {filename}\n")
        n_pred += 1


if __name__ == "__main__":
    main()
