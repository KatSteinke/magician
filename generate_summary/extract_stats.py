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
    # rename amount of marker genes for clarity
    checkm_stats = checkm_stats.rename(columns={'0': '0_markers', '1': '1_marker', '2': '2_markers', '3': '3_markers',
                                        '4': '4_markers', '5': '5_markers'})
    # unify bin naming for easier comparison
    checkm_stats["Bin Id"] = checkm_stats["Bin Id"].apply(lambda x: x.replace('.', '_'))
    return checkm_stats


def get_bins_and_genomes(mash_file: pathlib.Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract list of all reference genomes and bins from dRep's Mash clustering output (Mdb.csv).
    Arguments:
        mash_file: Path to Mash clustering file (Mdb.csv)
    Returns:
        All unique genomes, and all unique bins, used in the dRep clustering.
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

    # unify naming
    genomes['genome2'] = genomes['genome2'].apply(lambda x: x.replace('.fa', '').replace('.', '_'))
    bins['genome2'] = bins['genome2'].apply(lambda x: x.replace('.fa', '').replace('.', '_'))

    return genomes, bins


def get_drep_stats(mummer_file: pathlib.Path, all_bins: pathlib.Path) -> pd.DataFrame:
    """Extract ANI from comma-separated file produced by dRep's Mummer-based ANI step.
    Arguments:
        mummer_file:    Path to Ndb.csv
        all_bins:       Path to Mdb.csv
    Returns:
        Closest bin for each reference with ANI and coverage.
    """
    mummer_anis_original = pd.read_csv(mummer_file, sep=",")
    # fix typos
    mummer_anis_original = mummer_anis_original.rename(lambda x: x.replace('querry', 'query'), axis="columns")
    # unify genome/bin names for easier comparison with other stats
    mummer_anis_original['query'] = mummer_anis_original['query'].apply(lambda x:
                                                                        x.replace('.fa', '').replace('.', '_'))
    mummer_anis_original['reference'] = mummer_anis_original['reference'].apply(lambda x:
                                                                                x.replace('.fa', '').replace('.', '_'))

    # drop all lines that just compare a file with itself
    mummer_anis = mummer_anis_original.drop(mummer_anis_original[mummer_anis_original['query']
                                                                 == mummer_anis_original['reference']].index)

    # for each cluster, take everything where the reference contains '_bin_' - gets all relevant comparisons
    # while dropping the doubles
    mummer_anis = mummer_anis[mummer_anis['reference'].str.contains('_bin_', regex=False)]

    # cut df down to query, reference, ANI, query and reference coverage, and primary cluster
    mummer_anis = mummer_anis[['query', 'reference', 'ref_coverage', 'query_coverage', 'ani',
                               'primary_cluster']].reset_index(drop=True)

    # TODO: possible to get smarter with dropping duplicates & setting fields to NA to avoid using second file & joining?
    genomes, bins = get_bins_and_genomes(all_bins)

    # check which genomes are not in 'query' and append these
    genomes_without_bin = genomes[~genomes['genome2'].isin(mummer_anis['query'])].reset_index(drop=True)
    genomes_without_bin = genomes_without_bin.rename(columns={'genome2': 'query'})
    # get primary clusters for these:
    # inner join genomes_without_bin and unique queries & primary clusters from original Mummer dataframe
    # on query ID
    genomes_without_bin = pd.merge(genomes_without_bin,
                                   mummer_anis_original[['query', 'primary_cluster']].drop_duplicates(),
                                   on='query').reset_index(drop=True)

    mummer_anis = mummer_anis.append(genomes_without_bin, ignore_index=True)
    # check which bins are not in 'reference' and append these as well
    bins_without_ref = bins[~bins['genome2'].isin(mummer_anis['reference'])].reset_index(drop=True)
    bins_without_ref = bins_without_ref.rename(columns={"genome2": "reference"})
    # find primary cluster numbers here as well
    bins_without_ref = pd.merge(bins_without_ref,
                                mummer_anis_original[['reference', 'primary_cluster']].drop_duplicates(),
                                on='reference').reset_index(drop=True)

    mummer_anis= mummer_anis.append(bins_without_ref, ignore_index=True)

    return mummer_anis


def merge_mag_and_ref_stats(mag_stats: pd.DataFrame, ref_stats: pd.DataFrame) -> pd.DataFrame:
    """Concatenate two tables with statistics for the synthetic MAGs and the reference genomes,
    marking which of the two each is.
     Arguments:
         mag_stats: table containing statistics of synthetic MAGs
         ref_stats: table containing statistics of reference genomes
     Returns:
         A merged table with sources of genomes marked.
    """
    mag_stats['genome_type'] = 'synthetic_MAG'
    ref_stats['genome_type'] = 'reference'
    merged_stats = mag_stats.append(ref_stats, ignore_index=True)
    return merged_stats


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
    parser = ArgumentParser(description="Extract genome statistics from CheckM, dRep and bbstats output.")
    parser.add_argument("stats", action="store", help="Path to BBstats file giving stats of bins")
    parser.add_argument("genome_stats", action="store", help="Path to BBstats file giving stats of reference genomes")
    parser.add_argument("checkm", action="store", help="Path to CheckM file for bins")
    parser.add_argument("genome_checkm", action="store", help="Path to CheckM file for reference genomes")
    parser.add_argument("drep_mash", action="store", help="Path to dRep Mash file (Mdb.csv)")
    parser.add_argument("drep_mummer", action="store", help="Path to dRep Mummer file (Ndb.csv)")
    parser.add_argument("-o", "--outfile", action="store",
                        help="Name of Excel file to write to (recommended extension: .xlsx) (default: samplestats.xlsx)",
                        default="samplestats.xlsx")
    args = parser.parse_args()

    # get absolute paths
    stats = pathlib.Path(args.stats).resolve()
    original_stats = pathlib.Path(args.genome_stats).resolve()
    checkm = pathlib.Path(args.checkm).resolve()
    original_checkm = pathlib.Path(args.genome_checkm).resolve()
    mash = pathlib.Path(args.drep_mash).resolve()
    mummer = pathlib.Path(args.drep_mummer).resolve()
    outfile = pathlib.Path(args.outfile).resolve()

    drep_stats = get_drep_stats(mummer, mash)
    # Extract stats for both MAGs and original genomes
    bb_stats = get_bb_stats(stats)
    reference_stats = get_bb_stats(original_stats)
    # Identify dRep primary clusters for each
    bb_stats = pd.merge(bb_stats, drep_stats, how="left", left_on="bin_name", right_on="reference")[['bin_name',
                                                                                                     'scaf_bp',
                                                                                                     'gc_avg',
                                                                                                     'n_scaffolds',
                                                                                                     'n_contigs',
                                                                                                     'scaffold_L50',
                                                                                                     'scaffold_N50',
                                                                                                     'primary_cluster']]
    reference_stats = pd.merge(reference_stats, drep_stats, how="left", left_on="bin_name",
                               right_on="query")[['bin_name', 'scaf_bp', 'gc_avg', 'n_scaffolds', 'n_contigs',
                                                  'scaffold_L50', 'scaffold_N50', 'primary_cluster']]

    complete_stats = merge_mag_and_ref_stats(bb_stats, reference_stats)

    # Get CheckM stats for MAGs and original genomes
    checkm_stats = get_checkm_stats(checkm)
    reference_checkm = get_checkm_stats(original_checkm)
    complete_checkm = merge_mag_and_ref_stats(checkm_stats, reference_checkm)



    write_summaries(complete_stats, complete_checkm, drep_stats, outfile)