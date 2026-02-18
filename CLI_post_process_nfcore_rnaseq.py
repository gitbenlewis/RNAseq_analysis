#!/usr/bin/env python3
"""Post-process nf-core/rnaseq outputs with rnaseq_utils helpers
example usage ran from repo_root/scripts/ directory:
python post_process_nfcore_rnaseq.py \
    --nfcore-rnaseq-dir ../nfcore_outputs/rnaseq/results/ \
    --metadata-file ../data/sample_metadata.csv \
    --min-row-average 10 \
    --min-non-zeros-per-row 1 \
    --batch-key batch \
    --batch-label bulk_RNAseq_batch \
    --organism human \
    --h5ad-path ../nfcore_outputs/rnaseq/results/nfcore_bulk_rnaseq_adata.h5ad \
    --save-h5ad
."""

import argparse
from pathlib import Path
from typing import Optional, Sequence

from ._rnaseq_utils import (
    filter_zeros_n_rowaverage_of_nfcore_output_dir,
    make_adata_nfcore_rnaseq,
    make_gene_id2gene_name_file,
    extract_read_quantification_metrics_nfcore_rnaseq,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_RNASEQ_DIR = REPO_ROOT / "nfcore_outputs" / "rnaseq" / "results"
DEFAULT_H2M_FILE = REPO_ROOT / "config" / "nfcore_inputs" / "h2m_agg.csv"
DEFAULT_M2H_FILE = REPO_ROOT /"config"/ "nfcore_inputs" / "m2h_agg.csv"
DEFAULT_H5AD_PATH = REPO_ROOT / "nfcore_outputs" / "rnaseq" / "results" / "nfcore_bulk_rnaseq_adata.h5ad"


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate, filter, and convert nf-core/rnaseq outputs for downstream analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--nfcore-rnaseq-dir",
        type=Path,
        default=DEFAULT_RNASEQ_DIR,
        help="Directory containing nf-core/rnaseq results (expects star_salmon subdir).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Optional : Output directory for processed differential abundance results.If None, saves to nf-core RNAseq output dir.",
    )
    parser.add_argument(
        "--h2m-map-file",
        type=Path,
        default=DEFAULT_H2M_FILE,
        help="CSV mapping human gene IDs to mouse homologs (used when homologs enabled).",
    )
    parser.add_argument(
        "--m2h-map-file",
        type=Path,
        default=DEFAULT_M2H_FILE,
        help="CSV mapping mouse gene IDs to human homologs (used when homologs enabled).",
    )
    parser.add_argument(
        "--metadata-file",
        type=Path,
        default=None,
        help="Sample metadata CSV/TSV (required unless --skip-make-adata).",
    )
    parser.add_argument(
        "--h5ad-path",
        type=Path,
        default=DEFAULT_H5AD_PATH,
        help="Output path for the optional AnnData (.h5ad) file.",
    )
    parser.add_argument(
        "--min-row-average",
        type=float,
        default=10.0,
        help="Minimum row average threshold for filtering counts.",
    )
    parser.add_argument(
        "--min-non-zeros-per-row",
        type=int,
        default=1,
        help="Minimum number of non-zero samples per gene for filtering counts.",
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Column in metadata to use as the AnnData batch key.",
    )
    parser.add_argument(
        "--batch-label",
        default="bulk_RNAseq_batch",
        help="Value written to the AnnData batch column for this dataset.",
    )
    parser.add_argument(
        "--organism",
        choices=("human", "mouse"),
        default="human",
        help="Reference organism for homolog annotations.",
    )
    parser.add_argument(
        "--skip-also-save-post-process-results-to-src-dir",
        action="store_false",
        default=True,
        dest="also_save_post_process_results_to_src_dir",
        help="Do not also save gene_id2gene_name.csv and processed DESeq2 tables back to the nf-core source directories.",
    )
    parser.add_argument(
        "--skip-save-h5ad",
        action="store_false",
        default=True,
        dest="save_h5ad",
        help="Do not write the generated AnnData object to --h5ad-path.",
    )
    parser.add_argument(
        "--skip-make-gene-map",
        action="store_true",
        help="Skip rebuilding star_salmon/gene_id2gene_name.csv.",
    )
    parser.add_argument(
        "--skip-filter-counts",
        action="store_true",
        help="Skip filtering the nf-core count tables.",
    )
    parser.add_argument(
        "--skip-make-adata",
        action="store_true",
        help="Skip building the AnnData object from nf-core outputs.",
    )
    parser.add_argument(
        "--skip-extract-read-quantification-metrics-nfcore-rnaseq",
        action="store_true",
        help="Skip extracting read quantification metrics from nf-core RNAseq outputs.",
    )
    parser.add_argument(
        "--save-filtered-counts",
        action="store_true",
        help="Persist the filtered count tables (mirrors nf-core filenames with suffixes).",
    )

    parser.add_argument(
        "--no-add-homologs",
        action="store_false",
        dest="add_homologs",
        help="Disable merging homolog annotations when building the gene map.",
    )

    parser.set_defaults(add_homologs=True, save_h5ad=True, also_save_post_process_results_to_src_dir=True)
    return parser.parse_args(argv)


def run_pipeline(args: argparse.Namespace) -> None:
    nfcore_dir = str(args.nfcore_rnaseq_dir)
    output_dir = str(args.output_dir) if args.output_dir is not None else None
    h2m_file = str(args.h2m_map_file)
    m2h_file = str(args.m2h_map_file)

    if args.skip_make_gene_map:
        print("Skipping gene_id2gene_name generation (per flag).")
    else:
        print("Creating/refreshing gene_id2gene_name.csv file.")
        make_gene_id2gene_name_file(
            nfcore_output_rnaseq_dir=nfcore_dir,
            add_homologs=args.add_homologs,
            organism=args.organism,
            h2m_agg_file=h2m_file,
            m2h_agg_file=m2h_file,
            output_dir=output_dir,
            also_save_to_src_dir=args.also_save_post_process_results_to_src_dir,
        )

    if args.skip_filter_counts:
        print("Skipping count-table filtering (per flag).")
    else:
        print(
            "Filtering nf-core count tables (min_rowavg=%s, min_nonzero=%s)."
            % (args.min_row_average, args.min_non_zeros_per_row)
        )
        filter_zeros_n_rowaverage_of_nfcore_output_dir(
            nfcore_output_rnaseq_dir=nfcore_dir,
            min_rowaverage=args.min_row_average,
            min_non_zeros_perrow=args.min_non_zeros_per_row,
            save_output=args.save_filtered_counts,
            return_dfs=False,
            output_dir=output_dir,
            also_save_to_src_dir=args.also_save_post_process_results_to_src_dir,
        )

    if args.skip_make_adata:
        print("Skipping AnnData construction (per flag).")
    else:
        print("Building AnnData object from nf-core outputs.")
        make_adata_nfcore_rnaseq(
            nfcore_output_rnaseq_dir=nfcore_dir,
            metadata_file_path=str(args.metadata_file),
            batch_key=args.batch_key,
            batch_label=args.batch_label,
            save_h5ad=args.save_h5ad,
            output_h5ad_file_name=str(args.h5ad_path),
            output_dir=output_dir,
            also_save_to_src_dir=args.also_save_post_process_results_to_src_dir,
        )
        
    if args.skip_extract_read_quantification_metrics_nfcore_rnaseq:
        print("Skipping extraction of read quantification metrics (per flag).")
    else:
        extract_read_quantification_metrics_nfcore_rnaseq(
            nfcore_output_rnaseq_dir=nfcore_dir,
            save_output=True,
            output_dir=output_dir,
        )

def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    run_pipeline(args)


if __name__ == "__main__":
    main()
