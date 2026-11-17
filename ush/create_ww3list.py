"""
create_ww3list.py

PURPOSE:
    Create a WW3 model file list text file (ww3list.txt) for
    modelBuoy_collocation.py.

    The script searches model output directories based on:
        - user-defined base data directory
        - start date
        - end date
        - forecast cycles
        - filename pattern

    Existing files are written into the ww3list text file.
    Missing files are reported to screen but do not stop execution.

USAGE:
    python create_ww3list.py \
        -d /path/to/model/data \
        -s YYYYMMDD \
        -e YYYYMMDD \
        -c 00 06 12 18 \
        -p 'gfswave.t{cycle}z.bull_tar' \
        -o /path/to/output/directory \
        -f ww3list.txt

OUTPUT:
    Text file containing full paths of existing WW3 model files.

NOTE:
    The script assumes the following directory structure:

        <data-dir>/gfs.YYYYMMDD/CC/wave/station/<filename>

    where:
        YYYYMMDD = date
        CC       = cycle (00, 06, 12, 18)

    The filename pattern must contain:
        {cycle}

    Example:
        gfswave.t{cycle}z.bull_tar

AUTOR and DATE:
    05/12/2026: Ming Chen, first version

"""

import argparse
from datetime import datetime, timedelta
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create WW3 model file list for modelBuoy_collocation.py"
    )

    parser.add_argument("-d", "--data-dir", required=True,
                        help="Base model data directory, e.g. /scratch3/.../Data/gfsv16")

    parser.add_argument("-s", "--start-date", required=True,
                        help="Start date in YYYYMMDD format")

    parser.add_argument("-e", "--end-date", required=True,
                        help="End date in YYYYMMDD format, inclusive")

    parser.add_argument("-c", "--cycles", nargs="+", default=["00", "06", "12", "18"],
                        help="Forecast cycles, e.g. 00 06 12 18")

    parser.add_argument("-p", "--filename-pattern", default="gfswave.t{cycle}z.bull_tar",
                        help="Filename pattern. Use {cycle}, e.g. gfswave.t{cycle}z.bull_tar")

    parser.add_argument("-o", "--outdir", required=True,
                        help="Directory to write ww3list text file")

    parser.add_argument("-f", "--outfile", default="ww3list.txt",
                        help="Output list filename")

    return parser.parse_args()

def main():
    args = parse_args()

    data_dir = Path(args.data_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = outdir / args.outfile

    start = datetime.strptime(args.start_date, "%Y%m%d")
    end = datetime.strptime(args.end_date, "%Y%m%d")

    if end < start:
        raise ValueError("end-date must be >= start-date")

    files_written = 0
    missing_files = 0

    with outfile.open("w") as f:
        current = start
        while current <= end:
            yyyymmdd = current.strftime("%Y%m%d")

            for cycle in args.cycles:
                cycle = cycle.zfill(2)

                model_file = (
                    data_dir
                    / f"gfs.{yyyymmdd}"
                    / cycle
                    / "wave"
                    / "station"
                    / args.filename_pattern.format(cycle=cycle)
                )

                if model_file.is_file():
                    f.write(str(model_file) + "\n")
                    files_written += 1
                    print(f"FOUND   : {model_file}")
                else:
                    missing_files += 1
                    print(f"MISSING : {model_file}")

            current += timedelta(days=1)

    print("=====================================================")
    print(f"WW3 list created: {outfile}")
    print(f"Files written   : {files_written}")
    print(f"Files missing   : {missing_files}")
    print("=====================================================")


if __name__ == "__main__":
    main()
