import pathlib

from argparse import ArgumentParser
from typing import Tuple

import pandas as pd

# TODO: style pass!
def get_bb_stats(stats_file: pathlib.Path) -> pd.DataFrame:
    """Extract bin name, contig and scaffold counts, size, GC content, N50 and L50
    from a tab-separated file produced by bbtools's statswrapper.sh
    Arguments:
        stats_file: Path to the file containing assembly statistics
    Returns:
        A selection of the genome statistics.
    """
    stats = pd.read_csv(stats_file, sep="\t", usecols=['n_scaffolds', 'n_contigs', 'scaf_bp', 'scaf_N50', 'scaf_L50',
                                                       'gc_avg', 'filename'])[['filename', 'scaf_bp', 'gc_avg',
                                                                               'n_scaffolds', 'n_contigs', 'scaf_N50',
                                                                               'scaf_L50']]
    # N50 and L50 are switched around in BBTools - switch them back
    stats = stats.rename(columns={"filename": "bin_name", "scaf_N50": "scaffold_L50", "scaf_L50": "scaffold_N50"})
    # extract bin name from file name - TODO: make this more elegant than a double-replace?
    stats['bin_name'] = stats['bin_name'].apply(lambda x: (x.split('/')[-1].replace('.fa', '').replace('.', '_')))
    return stats


def get_checkm_stats(checkm_file: pathlib.Path) -> pd.DataFrame:
    """Extract taxonomy and marker genes from a tab-separated file produced by CheckM.
    Arguments:
        checkm_file: Path to the file containing CheckM results
    Returns:
        Taxonomy and absence/presence data for marker genes for all bins.
    """
    checkm_stats = pd.read_csv(checkm_file, sep="\t")
    # TODO: do we need to unify naming?
    # rename amount of marker genes for clarity
    checkm_stats = checkm_stats.rename(columns={'0': '0_markers', '1': '1_marker', '2': '2_markers', '3': '3_markers',
                                        '4': '4_markers', '5': '5_markers'})
    # drop all columns where every value is 0
    #checkm_stats = checkm_stats.loc[:, (checkm_stats != 0).any(axis=0)]
    return checkm_stats


def get_bins_and_genomes(mash_file: pathlib.Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Get list of genomes and bins from Mdb.csv
    TODO: improve documentation
    """
    # obtain a list of all genomes and bins by selecting the unique entries for genome_2 from Mdb.csv
    bins_and_genomes = pd.read_csv(mash_file, sep=",")
    # genomes are all that don't contain '.bin.'
    genomes = bins_and_genomes[~bins_and_genomes['genome2'].str.contains('.bin.',
                                                                         regex=False)].drop_duplicates(
        subset=["genome2"])[['genome2']].reset_index(drop=True)
    # bins are the opposite
    bins = bins_and_genomes[bins_and_genomes['genome2'].str.contains('.bin.',
                                                                     regex=False)].drop_duplicates(
        subset=["genome2"])[['genome2']].reset_index(drop=True)
    return genomes, bins


def get_drep_stats(mummer_file: pathlib.Path, all_bins: pathlib.Path) -> pd.DataFrame:
    """Extract ANI from comma-separated file produced by dRep's Mummer-based ANI step.
    Arguments:
        mummer_file:    Path to Ndb.csv
        all_bins:       Path to Mdb.csv
    Returns:
        Closest bin for each reference with ANI and coverage.
    """
    mummer_anis = pd.read_csv(mummer_file, sep=",")
    # fix typos?
    mummer_anis = mummer_anis.rename(lambda x: x.replace('querry', 'query'), axis="columns")
    # drop all lines that just compare a file with itself
    mummer_anis = mummer_anis.drop(mummer_anis[mummer_anis['query'] == mummer_anis['reference']].index)
    # for each cluster, take everything where the reference contains '.bin.'
    mummer_anis = mummer_anis[mummer_anis['reference'].str.contains('.bin.', regex=False)]
    # cut df down to query, reference, ANI, query and reference coverage
    mummer_anis = mummer_anis[['query', 'reference', 'ref_coverage', 'query_coverage', 'ani']].reset_index(drop=True)
    genomes, bins = get_bins_and_genomes(all_bins)
    # check which genomes are not in 'query' and append these with no closest bin, ....
    genomes_without_bin = genomes[~genomes['genome2'].isin(mummer_anis['query'])].reset_index(drop=True)
    genomes_without_bin = genomes_without_bin.rename(columns={'genome2': 'query'})
    mummer_anis = mummer_anis.append(genomes_without_bin, ignore_index=True)
    # check which bins are not in 'reference' and append these ditto
    bins_without_ref = bins[~bins['genome2'].isin(mummer_anis['reference'])].reset_index(drop=True)
    bins_without_ref = bins_without_ref.rename(columns={"genome2": "reference"})
    mummer_anis= mummer_anis.append(bins_without_ref, ignore_index=True)
    return mummer_anis


# final step: write it all to Excel
def write_summaries(bb_stats: pd.DataFrame, checkm_stats: pd.DataFrame, drep_stats: pd.DataFrame,
                    summary_file: str):
    """Collect all summary files in an Excel file.
    Arguments:
        bb_stats: Summary of bb_stats output
        checkm_stats: Summary of checkm output
        drep_stats: Summary of drep output
        summary_file: location to write summary file to
    """
    with pd.ExcelWriter(summary_file) as writer:
        bb_stats.to_excel(writer, sheet_name="BB_stats")
        checkm_stats.to_excel(writer, sheet_name="CheckM")
        drep_stats.to_excel(writer, sheet_name="dRep")


if __name__ == "__main__":
    # TODO: rework help!!!!
    parser = ArgumentParser(description="Extract genome statistics from CheckM, dRep and bbstats output.")
    parser.add_argument("stats", action="store", help="BBstats file for bins")
    parser.add_argument("genome_stats", action="store", help="BBstats file for input genomes")
    parser.add_argument("checkm", action="store", help="CheckM file")
    parser.add_argument("drep_mash", action="store", help="dRep Mash file")
    parser.add_argument("drep_mummer", action="store", help="dRep Mummer file")
    parser.add_argument("-o", "--outfile", action="store", help="name of excel file to write to",
                        default="samplestats.xlsx")
    args = parser.parse_args()

    stats = pathlib.Path(args.stats).resolve()
    original_stats = pathlib.Path(args.genome_stats).resolve()
    checkm = pathlib.Path(args.checkm).resolve()
    mash = pathlib.Path(args.drep_mash).resolve()
    mummer = pathlib.Path(args.drep_mummer).resolve()
    outfile = pathlib.Path(args.outfile).resolve()

    # Extract stats for both MAGs and original genomes
    bb_stats = get_bb_stats(stats)
    # mark the synthetic MAGs as such by adding a column
    bb_stats['genome_type'] = 'synthetic_MAG'
    # original genome stats for comparison
    reference_stats = get_bb_stats(original_stats)
    reference_stats['genome_type'] = 'reference'
    complete_stats = bb_stats.append(reference_stats, ignore_index=True)

    checkm_stats = get_checkm_stats(checkm)
    drep_stats = get_drep_stats(mummer, mash)

    write_summaries(complete_stats, checkm_stats, drep_stats, outfile)