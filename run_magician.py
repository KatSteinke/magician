import pathlib
import re
import subprocess

from argparse import ArgumentParser
from typing import Optional


def run_snakemake(target: str, profile_type: Optional[str] = "mbarc", profile_base: Optional[str] = "",
                  readlength: Optional[str] = "", insert_size: Optional[int] = 270,
                  *snake_params) -> None:
    """Runs a Snakemake command with optional configuration parameters.
    Arguments:
        target:         rule or output file for Snakemake
        profile_type:   type of ART error profile to use for CAMISIM
        profile_base:   path to own ART error profile for CAMISIM
        readlength:     read length used for own ART error profile
        insert_size:    mean insert size to use for simulation (default: 270)
    """
    # TODO: make testable
    # check all elements of the command
    for input_param in [target, profile_type, profile_base, readlength, insert_size]:
        bad_chars = re.search(r"""[^a-zA-Z0-9"'./_\- ]""", input_param)
        if bad_chars:
            raise ValueError("Invalid character(s) in one of the parameters: {}.".format(bad_chars.group(0)))

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

    subprocess.run(["snakemake", target, "-s", snake_path, "--config", 'profile_type="{}"'.format(profile_type),
                    'profile_name="{}"'.format(profile_base), 'readlength="{}"'.format(readlength),
                    'insert_size={}'.format(insert_size), *snake_params])
    # construct string
    #snakemake_cmd = 'snakemake {target} -s {snakefile} --config profile_type="{profile_type}" \
    #profile_name="{profile_base}" readlength="{read_length}"'.format(target=target, snakefile=snake_path,
    #                                                                 profile_type=profile_type,
    #                                                                 profile_base=profile_base,
    #                                                                 read_length=readlength)


if __name__ == "__main__":
    parser = ArgumentParser(description="Run MAGICIAN to create a specified file.")
    parser.add_argument("target", action="store", help="Desired output file or rule")
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
    parser.add_argument("--snake_flags", nargs='*', help="Flags to be passed to snakemake, enclosed in quotes")
    args = parser.parse_args()
    target_result = args.target
    profiletype = args.profile_type
    profilename = args.profile_name
    read_length = args.profile_readlength
    insert_size = args.insert_size
    snake_flags = args.snake_flags[0].split()

    run_snakemake(target_result, profiletype, profilename, read_length, insert_size, *snake_flags)

