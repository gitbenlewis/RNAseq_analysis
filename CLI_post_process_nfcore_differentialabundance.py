#!/usr/bin/env python3
"""Post-process nf-core/differentialabundance outputs via rnaseq_utils helpers.
example usage ran from repo_root/scripts/ directory:
python post_process_nfcore_differentialabundance.py \
    --nfcore-rnaseq-dir ~/projects/gitbenlewis/nfcore_outputs/rnaseq/results \
    --nfcore-da-dir ~/projects/gitbenlewis/nfcore_outputs/differentialabundance/results \
    --h2m-map-file ~/projects/gitbenlewis/data/h2m_agg.csv \
    --m2h-map-file ~/projects/gitbenlewis/data/m2h_agg.csv \
    --organism human    
verbose example usage with all options:
python post_process_nfcore_differentialabundance.py \
    --nfcore-rnaseq-dir ~/projects/gitbenlewis/nfcore_outputs/rnaseq/results \
    --nfcore-da-dir ~/projects/gitbenlewis/nfcore_outputs/differentialabundance/results \
    --h2m-map-file ~/projects/gitbenlewis/data/h2m_agg.csv \
    --m2h-map-file ~/projects/gitbenlewis/data/m2h_agg.csv \
    --raw-table-suffix .deseq2.results.tsv \
    --raw-table-separator '\t' \
    --processed-table-suffix .deseq2.results.csv \
    --concat-prefix differentialabundance \
    --organism human \
    --skip-make-gene-map \
    --skip-add-gene-names \
    --skip-rank-metric \
    --skip-concat \
    --no-add-homologs \
    --no-sort-rank \
    --no-save-concat       
"""


import argparse
from pathlib import Path
from typing import Optional, Sequence

from ._rnaseq_utils import (
    add_Rank_Metric_S_to_all_deseq2_tables_nfcore_differentialabundance,
    add_gene_names_all_to_deseq2_tables_nfcore_differentialabundance,
    concat_deseq2_files_nfcore_differentialabundance,
    make_gene_id2gene_name_file,
)




def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate and merge nf-core/differentialabundance DESeq2 tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--nfcore-rnaseq-dir",
        type=Path,
        #default=DEFAULT_RNASEQ_DIR,
        help="Directory containing nf-core/rnaseq results (expects star_salmon subdir).",
    )
    parser.add_argument(
        "--nfcore-da-dir",
        type=Path,
        #default=DEFAULT_DA_DIR,
        help="nf-core/differentialabundance results directory (expects tables/differential).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Optional : Output directory for processed differential abundance results.If None, saves to nf-core DA output dir.",
    )
    parser.add_argument(
        "--h2m-map-file",
        type=Path,
        #default=DEFAULT_H2M_FILE,
        help="CSV mapping human gene IDs to mouse homologs (used if homologs enabled).",
    )
    parser.add_argument(
        "--m2h-map-file",
        type=Path,
        #default=DEFAULT_M2H_FILE,
        help="CSV mapping mouse gene IDs to human homologs (used if homologs enabled).",
    )
    parser.add_argument(
        "--raw-table-suffix",
        default=".deseq2.results.tsv",
        help="Suffix for the raw nf-core DESeq2 result tables.",
    )
    parser.add_argument(
        "--raw-table-separator",
        default="\t",
        help="Delimiter used inside the raw nf-core DESeq2 tables.",
    )
    parser.add_argument(
        "--processed-table-suffix",
        default=".deseq2.results.csv",
        help="Suffix assigned to DESeq2 tables after adding gene names.",
    )
    parser.add_argument(
        "--concat-prefix",
        default="differentialabundance",
        help="Prefix for the concatenated DESeq2 CSV file.",
    )
    parser.add_argument(
        "--organism",
        choices=("human", "mouse"),
        default="human",
        help="Reference organism for homolog annotations.",
    )
    parser.add_argument(
        "--dataset-name",
        type=str,
        default="RNAseq Dataset",
        help="Name of the dataset for plot titles and annotations.",
    )
    parser.add_argument(
        "--skip-also-save-post-process-results-to-src-dir",
        action="store_false",
        default=True,
        dest="also_save_post_process_results_to_src_dir",
        help="Do not also save gene_id2gene_name.csv and processed DESeq2 tables back to the nf-core source directories.",
    )
    parser.add_argument(
        "--skip-make-gene-map",
        action="store_true",
        help="Assume gene_id2gene_name.csv already exists and skip rebuilding it.",
    )
    parser.add_argument(
        "--skip-add-gene-names",
        action="store_true",
        help="Skip adding gene name annotations to DESeq2 tables.",
    )
    parser.add_argument(
        "--skip-rank-metric",
        action="store_true",
        help="Skip computing the Rank_Metric_S column across DESeq2 tables.",
    )
    parser.add_argument(
        "--skip-concat",
        action="store_true",
        help="Skip concatenating the processed DESeq2 tables.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip generating volcano and MA plots from processed DESeq2 tables.",
    )
    parser.add_argument(
        "--no-add-homologs",
        action="store_false",
        dest="add_homologs",
        help="Disable merging homolog annotations when building the gene map.",
    )
    parser.add_argument(
        "--no-sort-rank",
        action="store_false",
        dest="sort_rank",
        help="Do not sort individual DESeq2 tables after computing Rank_Metric_S.",
    )
    parser.add_argument(
        "--no-save-concat",
        action="store_false",
        dest="save_concat",
        help="Do not write the concatenated CSV to disk (still returns dataframe).",
    )

    parser.set_defaults(add_homologs=True, sort_rank=True, save_concat=True, also_save_post_process_results_to_src_dir=True)
    return parser.parse_args(argv)


def run_pipeline(args: argparse.Namespace) -> None:
    nfcore_rnaseq_dir = str(args.nfcore_rnaseq_dir)
    nfcore_da_dir = str(args.nfcore_da_dir)
    output_dir = str(args.output_dir) if args.output_dir is not None else None
    h2m_file = str(args.h2m_map_file)
    m2h_file = str(args.m2h_map_file)

    if args.skip_make_gene_map:
        print("Skipping gene_id2gene_name generation (per flag).")
    else:
        print(f"Creating gene_id2gene_name.csv in {args.output_dir}")
        if args.also_save_post_process_results_to_src_dir:
            print(f"also Creating gene_id2gene_name.csv in {args.nfcore_rnaseq_dir / 'star_salmon'}")
        make_gene_id2gene_name_file(
            nfcore_output_rnaseq_dir=nfcore_rnaseq_dir,
            add_homologs=args.add_homologs,
            organism=args.organism,
            h2m_agg_file=h2m_file,
            m2h_agg_file=m2h_file,
            output_dir=output_dir,
            also_save_to_src_dir=args.also_save_post_process_results_to_src_dir
        )

    if args.skip_add_gene_names:
        print("Skipping gene name annotation of DESeq2 tables (per flag).")
    else:
        print("Adding gene names to nf-core/differentialabundance DESeq2 tables.")
        add_gene_names_all_to_deseq2_tables_nfcore_differentialabundance(
            nfcore_DA_output_dir=nfcore_da_dir,
            nfcore_output_rnaseq_dir=nfcore_rnaseq_dir,
            input_file_suffix=args.raw_table_suffix,
            input_file_separator=args.raw_table_separator,
            output_file_suffix=args.processed_table_suffix,
            run_make_gene_id2gene_name_file=True,
            add_homologs=args.add_homologs,
            organism=args.organism,
            h2m_map_file=h2m_file,
            m2h_map_file=m2h_file,
            output_dir=output_dir,
            also_save_to_src_dir=args.also_save_post_process_results_to_src_dir
        )

    if args.skip_rank_metric:
        print("Skipping Rank_Metric_S calculation (per flag).")
    else:
        print("Adding Rank_Metric_S to processed DESeq2 tables.")
        add_Rank_Metric_S_to_all_deseq2_tables_nfcore_differentialabundance(
            nfcore_DA_output_dir=nfcore_da_dir,
            input_file_suffix=args.processed_table_suffix,
            sort_local=args.sort_rank,
            output_dir=output_dir,
        )

    if args.skip_concat:
        print("Skipping concatenation of DESeq2 tables (per flag).")
    else:
        print("Concatenating annotated DESeq2 tables.")
        if output_dir is not None:
            print(f"Saving concatenated file to output dir: {output_dir}")

        concat_deseq2_files_nfcore_differentialabundance(
            nfcore_DA_output_dir=nfcore_da_dir,
            non_nfcore_DA_output_dir=output_dir,
            input_file_suffix=args.processed_table_suffix,
            output_file_prefix=args.concat_prefix,
            save_output=args.save_concat,
            output_dir=output_dir,
            also_save_to_src_dir=args.also_save_post_process_results_to_src_dir
        )
    if args.skip_plots:
        print("Skipping generation of volcano and MA plots (per flag).")
    else:
        from . import _rnaseq_utils as rnaseq_utils
        rnaseq_utils.make_volcano_MA_plots_from_nf_deseq_output_dir(
        output_dir,
        file_sufix='.deseq2.results.csv',
        plot_dataset_title=args.dataset_name,
        log2FoldChange_threshold=0.5, ylimit_volcano=None,xlimit_volcano_l2fc=None, xlimit_MAonly=None,
        label_gene_name=False,figsize=(12,12),
        save_plots=True,
        save_dir='std_plots'
        )


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    run_pipeline(args)

if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------
# Example usage:
# python scripts/post_process_nfcore_differentialabundance.py \
#   --nfcore-rnaseq-dir ~/projects/gitbenlewis/nfcore_outputs/rnaseq/results \
#   --nfcore-da-dir ~/projects/gitbenlewis/nfcore_outputs/differentialabundance/results \
#   --organism human
# ---------------------------------------------------------------------
# Verbose usage:
# python scripts/post_process_nfcore_differentialabundance.py \
#   --nfcore-rnaseq-dir ~/projects/gitbenlewis/nfcore_outputs/rnaseq/results \
#   --nfcore-da-dir ~/projects/gitbenlewis/nfcore_outputs/differentialabundance/results \
#   --h2m-map-file ~/projects/gitbenlewis/data/h2m_agg.csv \
#   --m2h-map-file ~/projects/gitbenlewis/data/m2h_agg.csv \
#   --raw-table-suffix .deseq2.results.tsv \
#   --raw-table-separator $'\\t' \
#   --processed-table-suffix .deseq2.results.csv \
#   --concat-prefix differentialabundance \
#   --organism human
# ---------------------------------------------------------------------
