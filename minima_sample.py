#!sage -python

PRIME_ALIASES = {
    62: 55340232221128654847,
    248: 2261564242916331941866620800950935700259179388000792266395655937654553313279,
}

COLUMNS = (
    ['I1', 'I2', 'I3', 'I4']
    + ['OL2', 'OL3', 'OL4']
    + ['OR2', 'OR3', 'OR4']
    + [f'inner{i}{j}' for i in range(1, 5) for j in range(1, 5)]
)


def random_ideals(p):
    import os, subprocess, threading, time

    num = os.cpu_count()

    mtx = threading.Lock()
    lines = []

    def fun():
        with subprocess.Popen( 
                    ['sage', '-python', './sample.py', str(p), '--birandom'],
                    stdout = subprocess.PIPE,
                    stderr = subprocess.DEVNULL,
                    env = os.environ | {'PYTHONOPTIMIZE': '1'},
                ) as proc:

            while True:
                try:
                    line = proc.stdout.readline().decode().strip()
                except Exception as e:
                    print(f'\x1b[31mthread crashed: {e}\x1b[0m')
                    break
                with mtx:
                    lines.append(line)

    threads = []
    while len(threads) < num:
        threads.append(threading.Thread(target=fun, daemon=True))
    for t in threads:
        t.start()

    while True:
        time.sleep(float(.1))
        with mtx:
            while lines:
                line = lines.pop(0)
                yield line

def main(p, howmany=None, output=None):
    from datetime import datetime
    from tqdm import tqdm

    if howmany is None:
        howmany = 10
    if output is None:
        timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')
        output = f'data/{p}-birandom-{timestamp}.csv'

    records = []
    interrupted = False
    try:
        with tqdm(total=howmany, desc='Collecting', unit='record') as progress:
            for record in random_ideals(p):
                records.append(record)
                progress.update()
                if len(records) >= howmany:
                    break
    except KeyboardInterrupt:
        interrupted = True
        print(f'\ninterrupted; saving {len(records)} records collected so far')

    save_records(p, output, records)
    if interrupted:
        raise SystemExit(130)


def save_records(p, output, records):
    with open(output, 'w', newline='') as f:
        prime = PRIME_ALIASES.get(p, p)
        f.write(f'# prime={prime}; columns={",".join(COLUMNS)}\n')
        for record in records:
            f.write(record)
            f.write('\n')
    print(f'saved {len(records)} records to {output}')

if __name__ == '__main__':
    import sys
    p = int(sys.argv[1])
    howmany = int(sys.argv[2]) if len(sys.argv) > 2 else None
    output = sys.argv[3] if len(sys.argv) > 3 else None
    main(p, howmany, output)
