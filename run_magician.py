import pathlib
import re
import subprocess

from argparse import ArgumentParser
from typing import List, Optional

def get_snake_cmd(input_file, target: str, profile_type: Optional[str] = "mbarc", profile_base: Optional[str] = "",
                  readlength: Optional[str] = "", insert_size: Optional[int] = 270, cluster_cmd: Optional[str] = "",
                  *snake_params) -> List[str]:
    """Get the Snakemake command with optional configuration parameters.
    Arguments:
        input_file:     File with paths to source genomes, sequence type (plasmid/chromosome) and desired relative
                        copy number in each community to simulate
        target:         rule or output file for Snakemake
        profile_type:   type of ART error profile to use for CAMISIM
        profile_base:   path to own ART error profile for CAMISIM
        readlength:     read length used for own ART error profile
        insert_size:    mean insert size to use for simulation (default: 270)
        cluster_cmd:    command to use for cluster mode
        snake_params:   parameters to pass to the Snakefile

    """
    # check all elements of the command
    for input_param in [target, profile_type, profile_base, readlength, insert_size]:
        bad_chars = re.search(r"""[^a-zA-Z0-9"'./_\- ]""", str(input_param))
        if bad_chars:
            raise ValueError("Arguments can only consist of alphanumeric characters, quote marks, ., /, _, - and space.")

    # sanity check parameters
    if insert_size <= 0:
        raise ValueError("Insert size needs to be above 0.")

    if profile_type == "own" and not (profile_base and readlength):
        raise ValueError("Both name of the custom error profile and read length of the error profile must be given when using own profiles.")

    if (profile_base or readlength) and profile_type != "own":
        raise ValueError("Name of the error profile and read length can only be specified when using own profiles.")

    # if custom parameters haven't been given, set to False for later processing in Snakefile
    if not profile_base:
        profile_base = "False"
    if not readlength:
        readlength = "False"

    # get path for snakefile
    snake_path = pathlib.Path(__file__).resolve().parent / "snakefiles" / "Snakefile"

    # get basic command
    snakemake_cmd = ["snakemake", target, "-s", snake_path, "--config", 'profile_type="{}"'.format(profile_type),
                     'profile_name="{}"'.format(profile_base), 'readlength="{}"'.format(readlength),
                     'insert_size={}'.format(insert_size),
                     'samples_file={}'.format(input_file),
                     *snake_params]

    # if cluster mode is specified:
    if cluster_cmd:
        snakemake_cmd.insert(4, "--cluster")
        snakemake_cmd.insert(5, cluster_cmd)

    return snakemake_cmd
    

if __name__ == "__main__":
    parser = ArgumentParser(description="Run MAGICIAN to simulate MAGs for a specified community"
                                        " or set of communities.")
    parser.add_argument("community_file", action="store",
                        help="File with paths to source genomes, sequence type (plasmid/chromosome)"
                             " and their desired relative copy number in each community to simulate "
                             "(one column with organisms' copy numbers per community)")
    parser.add_argument("--target", action="store",
                        help="Desired output file or rule "
                             "(default: MAGs, statistics and summary files for all communities)",
                        default="all_bin_summaries")
    parser.add_argument("--profile_type", action="store",
                        help="Type of ART error profile to use for CAMISIM: mbarc, hi, mi, hi150, own (default: mbarc)",
                        default="mbarc", choices=['mbarc', 'hi', 'mi', 'hi150', 'own'])
    parser.add_argument("--profile_name", action="store",
                        help="""Base name of custom error profile, if given (name of files without '[1/2].txt');\
                         required with 'own' error profile""", default="")
    parser.add_argument("--profile_readlength", action="store",
                        help="Read length of custom error profile; required with 'own' error profile", default="")
    parser.add_argument("--insert_size", action="store", type=int, default=270,
                        help="Mean insert size for read simulation (default: 270)")
    parser.add_argument("--cluster", action="store", default="",
                        help="""For use with snakemake's cluster mode; supply command for submitting jobs as you \
                        would with snakemake.""")
    parser.add_argument("--snake_flags", nargs='*', help="Flags to be passed to snakemake, enclosed in quotes")
    args = parser.parse_args()
    community_file = pathlib.Path(args.community_file)
    target_result = args.target
    profiletype = args.profile_type
    profilename = args.profile_name
    read_length = args.profile_readlength
    insert_size = args.insert_size
    cluster_cmd = args.cluster
    snake_flags = args.snake_flags[0].split()

    snake_command = get_snake_cmd(pathlib.Path(__file__).parent / "data" / "test_distribution_file.tsv", target_result,
                                  profiletype, profilename, read_length, insert_size, cluster_cmd, *snake_flags)
    subprocess.run(snake_command)

