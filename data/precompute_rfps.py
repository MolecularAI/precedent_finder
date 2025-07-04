"""Module containing routines to precompute fingerprints"""

import argparse
from typing import Optional, Sequence

from scipy import sparse
import pandas as pd
import numpy as np
from rdkit import RDLogger

from rxnutils.chem.reaction import ChemicalReaction
from rxnutils.chem.features.reaction_centre_fp import ReactionCentreFingerprint
from rxnutils.chem.features.simple_rdkit import SimpleRdkitFingerprint

rd_logger = RDLogger.logger()
rd_logger.setLevel(RDLogger.CRITICAL)


def _calc_global_fps(data: pd.Series, numbits: int) -> sparse.csr_array:
    fp_factory = SimpleRdkitFingerprint("fingerprint_mixed", numbits=numbits)
    inputs = [fp_factory(ChemicalReaction(rsmi, clean_smiles=False)) for rsmi in data]
    return sparse.lil_matrix(inputs, dtype=np.int8).tocsr()


def _calc_center_fps(data: pd.Series, numbits: int, max_centers: int) -> sparse.csr:
    fp_factory = ReactionCentreFingerprint(numbits=numbits, max_centers=max_centers)
    inputs = [fp_factory(ChemicalReaction(rsmi, clean_smiles=False)) for rsmi in data]
    return sparse.lil_matrix(inputs, dtype=np.int8).tocsr()


def main(args: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--fp_path", required=True)
    parser.add_argument("--rcfp_path", required=True)
    parser.add_argument("--rsmi_column", default="ReactionSmilesClean")
    parser.add_argument("--fp_len", type=int, default=1024)
    parser.add_argument("--max_centers", type=int, default=6)
    args = parser.parse_args(args)

    metadata = pd.read_csv(args.data, sep="\t", usecols=[args.rsmi_column])
    metadata = metadata.squeeze("columns")

    fp_global = _calc_global_fps(metadata, numbits=args.fp_len)
    sparse.save_npz(args.fp_path, fp_global, compressed=True)

    fp_center = _calc_center_fps(
        metadata, numbits=args.fp_len, max_centers=args.max_centers
    )
    sparse.save_npz(args.rcfp_path, fp_center, compressed=True)


if __name__ == "__main__":
    main()
