#!/usr/bin/env python3

import argparse
import io
import re
import zipfile
from datetime import datetime
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description='Merge timestamped birandom CSV files for a given p.'
    )
    parser.add_argument('p', help='Filename prefix, for example 62')
    parser.add_argument(
        '-o',
        '--output',
        type=Path,
        help='Output ZIP path (default: data/{p}-birandom-merged-YYYYMMDD-hhmmss.csv.zip)',
    )
    args = parser.parse_args()

    data_dir = Path('data')
    pattern = re.compile(
        rf'^{re.escape(args.p)}-birandom-\d{{8}}-\d{{6}}\.csv$'
    )
    sources = sorted(
        path for path in data_dir.iterdir()
        if path.is_file() and pattern.fullmatch(path.name)
    )

    if not sources:
        parser.error(
            f'no files match data/{args.p}-birandom-YYYYMMDD-hhmmss.csv'
        )

    timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')
    output = args.output or Path(
        f'data/{args.p}-birandom-merged-{timestamp}.csv.zip'
    )
    if output.suffix != '.zip':
        output = output.with_name(output.name + '.zip')
    csv_name = output.name.removesuffix('.zip')

    records = 0
    with zipfile.ZipFile(
        output, 'w', compression=zipfile.ZIP_DEFLATED
    ) as archive:
        with archive.open(csv_name, 'w') as compressed:
            with io.TextIOWrapper(compressed, newline='') as destination:
                for source in sources:
                    with source.open(newline='') as data:
                        for line in data:
                            if (
                                line.strip()
                                and not line.lstrip().startswith('#')
                            ):
                                destination.write(
                                    line.rstrip('\r\n') + '\n'
                                )
                                records += 1

    print(f'merged {len(sources)} files and {records} records into {output}')


if __name__ == '__main__':
    main()
