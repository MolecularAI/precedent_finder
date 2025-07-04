"""
Module containing pipeline for setting up the data necessary for the Precedent Finder
"""

from pathlib import Path

from rxnutils.data.batch_utils import (
    create_csv_batches,
    combine_csv_batches,
    combine_sparse_matrix_batches,
)
from rxnutils.pipeline.runner import main as validation_runner
from metaflow import step, FlowSpec, Parameter

from precompute_rfps import main as precompute_fps


class PrecedentFinderDataSetupFlow(FlowSpec):
    """Pipeline for setting up data for the Precedent Finder"""

    nbatches = Parameter("nbatches", type=int, required=True)
    folder = Parameter("folder", default=".")

    @step
    def start(self):
        """Setup data processing"""
        # pylint: disable=attribute-defined-outside-init
        self.partitions = create_csv_batches(
            str(Path(self.folder) / "uspto_data_mapped.csv"),
            self.nbatches,
            str(Path(self.folder) / "metadata_cleaned.csv"),
        )
        self.next(self.do_data_processing, foreach="partitions")

    @step
    def do_data_processing(self):
        """Perform data processing of USPTO data"""
        idx, start, end = self.input
        pipeline_path = str(Path(__file__).parent / "clean-up.yml")
        if idx > -1:
            validation_runner(
                [
                    "--pipeline",
                    pipeline_path,
                    "--data",
                    str(Path(self.folder) / "uspto_data_mapped.csv"),
                    "--output",
                    str(Path(self.folder) / f"metadata_cleaned.csv.{idx}"),
                    "--batch",
                    str(start),
                    str(end),
                    "--max-workers",
                    "1",
                    "--no-intermediates",
                ]
            )
            precompute_fps(
                [
                    "--data",
                    str(Path(self.folder) / f"metadata_cleaned.csv.{idx}"),
                    "--fp_path",
                    str(Path(self.folder) / f"rfp_precomp.{idx}.npz"),
                    "--rcfp_path",
                    str(Path(self.folder) / f"rcfp_precomp.{idx}.npz"),
                ]
            )
        self.next(self.join_batches)

    @step
    def join_batches(self, _):
        """Combine batches of data"""
        combine_csv_batches(
            str(Path(self.folder) / "metadata_cleaned.csv"), self.nbatches
        )
        combine_sparse_matrix_batches(
            str(Path(self.folder) / "rfp_precomp.npz"), self.nbatches
        )
        combine_sparse_matrix_batches(
            str(Path(self.folder) / "rcfp_precomp.npz"), self.nbatches
        )
        self.next(self.end)

    @step
    def end(self):
        """Final step, just print information"""
        print(
            f"Processed metadata file is locate here: {Path(self.folder) / 'metadata_cleaned.csv'}\n"
            f"Pre-computed global reaction FPs is locate here: {Path(self.folder) / 'rfp_precomp.npz'}\n"
            f"Pre-computed reaction center FPs is locate here: {Path(self.folder) / 'rcfp_precomp.npz'}"
        )


if __name__ == "__main__":
    PrecedentFinderDataSetupFlow()
