import pathlib

from argparse import ArgumentParser

import pandas as pd


def create_comparison_table(base_table: pathlib.Path) -> pd.DataFrame:
    """Generate an overview table comparing bins to their closest reference genomes.
    Arguments:
        base_table: Path to the Excel file summarizing various QC parameters.
    Returns:
        A table directly comparing bins to closest reference genomes.
    """
    # Read in data
    bb_stats = pd.read_excel(base_table, sheet_name="BB_stats", index_col=0)
    checkm_stats = pd.read_excel(base_table, sheet_name="CheckM", index_col=0)
    drep_stats = pd.read_excel(base_table, sheet_name="dRep", index_col=0)

    summary_table = bb_stats[bb_stats["genome_type"] == "synthetic_MAG"][['bin_name']]

    # closest organism: SELECT drep_stats.query AS closest_genome, summary_table.bin_name, MAX(drep_stats.ani)
    # FROM summary table
    # LEFT JOIN drep_stats ON summary_table.bin_name = drep_stats.reference
    # GROUP BY summary_table.bin_name
    # closest organism: find highest ANI for each bin in drep_stats, filter down to only these
    # then left join (leaving all entries from the summary table intact) summary table on filtered drep_stats
    summary_table = pd.merge(summary_table, drep_stats[drep_stats.groupby("reference")['ani'].transform(max)
                                                       == drep_stats['ani']],
                            left_on='bin_name', right_on='reference', how='left')[["bin_name", "query"]]
    summary_table = summary_table.rename(columns={"query": "closest_genome"})
    # scaffold etc. difference
    # create temp columns - one with scaffold stats, one with source genome stats
    # TODO: this is ugly as sin, FIX IT
    summary_bb = pd.merge(summary_table, bb_stats, how="left", on="bin_name")
    summary_bb = pd.merge(summary_bb, bb_stats, left_on="closest_genome", right_on="bin_name", how="left",
                          suffixes=('_bin', '_ref'))
    summary_bb["scaffold_difference"] = summary_bb['n_scaffolds_bin'] - summary_bb["n_scaffolds_ref"]
    summary_bb["contig_difference"] = summary_bb['n_contigs_bin'] - summary_bb['n_contigs_ref']
    summary_bb["length_difference"] = summary_bb["scaf_bp_bin"] - summary_bb["scaf_bp_ref"]
    summary_bb["gc_difference"] = summary_bb["gc_avg_bin"] - summary_bb["gc_avg_ref"]

    summary_checkm = pd.merge(pd.merge(summary_table, checkm_stats, left_on="bin_name", right_on="Bin Id", how="left"),
                              checkm_stats, how="left", left_on="closest_genome", right_on="Bin Id",
                              suffixes=("_bin", "_ref"))
    summary_checkm["unique_markers_difference"] = summary_checkm["1_marker_bin"] - summary_checkm["1_marker_ref"]
    # to see how many marker genes were found at all: everything that wasn't *not* found (simplify?)
    summary_checkm["all_markers_difference"] = (summary_checkm["# markers_bin"] - summary_checkm["0_markers_bin"]) \
                                               - (summary_checkm["# markers_ref"] - summary_checkm["0_markers_ref"])

    summary_table = summary_bb[["bin_name_bin", "closest_genome", "scaffold_difference", "contig_difference",
                                "length_difference", "gc_difference"]]
    summary_table = summary_table.rename(columns={"bin_name_bin": "bin_name"})
    summary_table = pd.merge(summary_table, summary_checkm[["Bin Id_bin",
                                                            "unique_markers_difference",
                                                            "all_markers_difference"]],
                             how="left", left_on="bin_name", right_on="Bin Id_bin").drop("Bin Id_bin", axis=1)
    return summary_table

# TODO: this also needs an explanation sheet at some point


if __name__ == "__main__":
    parser = ArgumentParser(description="Generate a bin-focused overview from the general MAGICIAN summary file")
    parser.add_argument("infile", action="store", help="Path to summary .xlsx file")
    parser.add_argument("-o", "--outfile", action="store", help="Name of output file (default: bin_summary.xlsx)",
                        default="bin_summary.xlsx")
    args = parser.parse_args()

    infile = pathlib.Path(args.infile).resolve()
    outfile = pathlib.Path(args.outfile).resolve()

    summary = create_comparison_table(infile)

    with pd.ExcelWriter(outfile) as outfile_writer:
        summary.to_excel(outfile_writer, sheet_name="summary")
