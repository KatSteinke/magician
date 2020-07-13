#import os

from argparse import ArgumentParser
from pathlib import Path
from typing import List

from Bio import SeqIO


def create_camisim_files(genbanks: List[str]):
    """Generate metadata and ID mapping files for CAMISIM and convert genbanks to FASTA.
        Arguments:
            genbanks: paths of the genbank files to use for the CAMISIM run
        Returns:
            A tab-separated metadata file containing genome ID, OTU, NCBI taxid and novelty category,
            a fasta file with the sequence, and
            a tab-separated file listing the genome ID and path of the fasta file
            for each genbank.
    """
    genbanks = [Path(genbank) for genbank in genbanks]
    # set up metadata file with "genome_ID\tOTU\tNCBI_ID\tnovelty_category"
    metadata = ["genome_ID\tOTU\tNCBI_ID\tnovelty_category"]
    # initialize id_to_genome_file to blank
    id_to_genome = []
    # initialize OTU count to 1 as we want all to be included
    otu_count = 1
    # check if dirs exist
    for camisim_dir in [Path('camisim_fasta'), Path('camisim_configfiles')]:
        if not camisim_dir.is_dir():
            camisim_dir.mkdir()
    # for each input file:
    for genbank in genbanks:
        record =  SeqIO.read(genbank, "genbank")# TODO: this only works with no contigs - adapt for contigs later!
        #     extract ID: version 
        record_id = record.id.replace(".", "_")
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
        # create and save FASTA file under {original_filename}.fa
        fasta_path = Path('camisim_fasta', '{}.fa'.format(genbank.stem)).resolve() # resolves as far as possible, appends the rest
        SeqIO.convert(genbank, "genbank", fasta_path, "fasta")
        # id to file line
        id_to_genome.append("{}\t{}".format(record_id, fasta_path))

    # write metadata file - TODO: separate dir
    with open(Path("camisim_configfiles", "metadata"), "w") as meta_file:
        meta_file.write("\n".join(metadata))
    # write id_to_genome_file
    with open(Path("camisim_configfiles","id_to_genome_file"), "w") as id_file:
        id_file.write("\n".join(id_to_genome))

if __name__ == "__main__":
    # input: list of files
    parser = ArgumentParser(
        description="Prepare CAMISIM metadata and id to file mapping files and create FASTA files from gbks")
    parser.add_argument("gbk_files", nargs="+", help="Genbank files for CAMISIM run")
    # TODO: add default names for metadata etc files
    args = parser.parse_args()
    gbk_files = args.gbk_files
    create_camisim_files(gbk_files)
