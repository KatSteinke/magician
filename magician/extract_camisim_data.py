#import os

from argparse import ArgumentParser
from pathlib import Path
from typing import Optional, List, Tuple, Union

import pandas as pd

from Bio import SeqIO


def get_metadata_from_records(records: List[str], make_abundance: Optional[bool]=False) -> Tuple[List[str],
                                                                                                List[str],
                                                                                                List[str],
                                                                                                Union[List[str], bool]]:
    """Extract metadata and genome IDs from a list of Genbank files, and convert the files
     to FASTA.
     Arguments:
         records:           Paths to the respective genbank files
         make_abundance:    flag for generating a tsv file specifying even abundance for genomes
     Returns:
         A list of genome IDs,
         a list of metadata lines,
         a list of genome ID to fasta file lines, and
         optionally a list of genome ID to abundance.
         (Fasta files are also written.)
     """
    records = [Path(record) for record in records]
    # set up metadata file with "genome_ID\tOTU\tNCBI_ID\tnovelty_category"
    metadata = ["genome_ID\tOTU\tNCBI_ID\tnovelty_category"]
    # initialize OTU count to 1 as we want all to be included
    otu_count = 1
    # initialize id_to_genome_file to blank
    id_to_genome = []
    # save record IDs separately as well
    record_ids = []
    # set id_to_abundance to False by default unless requested otherwise
    id_to_abundance = False
    # if abundance is to be generated:
    if make_abundance:
        common_abundance = 1 / len(records)
        id_to_abundance = []
    # check if fasta dir exists
    if not Path('camisim_fasta').is_dir():
        Path('camisim_fasta').mkdir()
    # for each input file:
    for genbank in records:
        record = SeqIO.read(genbank, "genbank")  # TODO: this only works with no contigs - adapt for contigs later!
        #     extract ID: version
        record_id = record.id.replace(".", "_")
        record_ids.append(record_id)
        #     extract taxonomic ID - there's no pretty way, so stateful parsing time
        taxon_id = ""
        with open(genbank, "r") as infile:
            line = infile.readline()
            while line and not taxon_id:
                if "taxon:" in line:
                    taxon_parts = line.strip().split(":")
                    taxon_id = taxon_parts[1].replace('"', "")
                line = infile.readline()
        # create and append metadata line
        metadata.append("{}\t{}\t{}\tknown_strain".format(record_id, otu_count, taxon_id))
        otu_count += 1
        # create and save FASTA file under {original_filename}.fa - not pretty, but avoids double loop
        fasta_path = Path('camisim_fasta',
                          '{}.fa'.format(genbank.stem)).resolve()  # resolves as far as possible, appends the rest
        SeqIO.convert(genbank, "genbank", fasta_path, "fasta")
        # id to file line
        id_to_genome.append("{}\t{}".format(record_id, fasta_path))
        if make_abundance:
            id_to_abundance.append("{}\t{}".format(record_id, common_abundance))

    return record_ids, metadata, id_to_genome, id_to_abundance


def write_camisim_files(metadata: List[str], id_to_genome: List[str], id_to_abundance: Union[List[str], bool],
                        metafile_name: Optional[str] = "metadata",
                        idfile_name: Optional[str]="id_to_genome_file") -> None:
    """Generate metadata and ID mapping files for CAMISIM and convert genbanks to FASTA.
        Arguments:
            metadata:           tab-separated strings representing metadata of genome files
            id_to_genome:       tab-separated strings representing genome ID and location of fasta files
            metafile_name:      optional name for metadata file to be created
            idfile_name:        optional name for id to fasta file to be created
            id_to_abundance:    if supplied, tab-separated strings mapping genome ID to uniform distributions,
                                else False
        Returns:
            Writes:
            A tab-separated metadata file containing genome ID, OTU, NCBI taxid and novelty category,
            a fasta file with the sequence,
            a tab-separated file listing the genome ID and path of the fasta file, and
            optionally, a tab-separated file listing genome ID and abundance
            for each genbank.
    """
    # check if dir exists
    if not Path('camisim_configfiles').is_dir():
        Path('camisim_configfiles').mkdir()
    # write metadata file
    with open(Path("camisim_configfiles", metafile_name), "w") as meta_file:
        meta_file.write("\n".join(metadata))
    # write id_to_genome_file
    with open(Path("camisim_configfiles", idfile_name), "w") as id_file:
        id_file.write("\n".join(id_to_genome))
    if id_to_abundance: # TODO: do we need this? Can we clean it up?
        with open(Path("camisim_configfiles", "id_to_distributions"), "w") as abundance_file:
            abundance_file.write("\n".join(id_to_abundance))

def get_camisim_per_sample(samples_file: Path, sample_col: str):
    """From a tab-separated table giving genbank files and their abundance in a given sample,
    create CAMISIM metadata, genome and abundance files.
    Arguments:
        samples_file:   Path to .tsv file
        sample_col:     column name for sample
    Returns:
        A tab-separated metadata file containing genome ID, OTU, NCBI taxid and novelty category,
        a fasta file with the sequence,
        a tab-separated file listing genome ID and abundance, and
        a tab-separated file listing the genome ID and path of the fasta file
        for each genbank listed in the table.
    """
    samples_table = pd.read_csv(samples_file, sep="\t", index_col=False)
    # create metadata and fasta files
    genome_ids, metadata, id_to_genome, id_to_abundance = get_metadata_from_records(samples_table['genomes'])
    # write metadata/id to fasta files
    write_camisim_files(metadata, id_to_genome, id_to_abundance, "metadata_{}".format(sample_col),
                        "id_to_genome_file_{}".format(sample_col))
    # extract genome IDs
    samples_table['genome_id'] = genome_ids
    with open(Path("camisim_configfiles", "id_to_distributions_{}".format(sample_col)), "w") as abundance_file:
        abundance_file.write(samples_table.to_csv(sep="\t", header=False, index=False, columns=['genome_id', sample_col]))



if __name__ == "__main__":
    # input: list of files
    parser = ArgumentParser(
        description="Prepare CAMISIM metadata and id to file mapping files and create FASTA files from gbks")
    #parser.add_argument("gbk_files", nargs="+", help="Genbank files for CAMISIM run")
    #parser.add_argument("--even_abundance", action="store_true",
    #                    help="Automatically generate abundance file with even distribution")
    # TODO: add default names for metadata etc files
    parser.add_argument("sample_file", help="Tab-separated file with abundances for samples")
    parser.add_argument("sample_column", help="Column name of sample to extract")
    args = parser.parse_args()
    sample_file = args.sample_file
    sample_column = args.sample_column
    #abundance = args.even_abundance
    #create_camisim_files(gbk_files, abundance)
    get_camisim_per_sample(sample_file, sample_column)
