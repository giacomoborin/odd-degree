#!/usr/bin/env python3

def random_ideals(args):
    import os, subprocess, threading, time, ast

    num = os.cpu_count()

    mtx = threading.Lock()
    lines = []

    def fun():
        with subprocess.Popen(
                    ['./sample.py'] + args,
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
                print(line, flush=True)
                yield line

if __name__ == '__main__':
    import sys
    print(f'# {sys.argv}')
    howmany = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    args = sys.argv[2:]
    results = list(zip(range(howmany), random_ideals(args)))
