#!/usr/bin/env python3
from math import sqrt

#def all_ideals(p):
#    B = BrandtModule(p)
#
#    Quat,(i,j,k) = B.quaternion_algebra().objgens()
#    print(f'{Quat = }')
#
#    O0 = B.maximal_order()
#    print(f'{O0 = }')
#    assert i in O0 and j in O0 and k in O0
#
#    yield from iter(B.right_ideals())

def random_ideal(O0):
    p = O0.discriminant()
    l = 257
    steps = ceil(log(p << 256, l))  #TODO
    I = O0*1
    for t in reversed(range(steps)):
        O = I.right_order()
        while True:
            α = sum(randrange(l)*g for g in O.basis())
            if α.reduced_norm().valuation(l) == 1:
                break
        assert α in O
        J = O*l + O*α
#        assert J.left_order() == O
        assert J.norm() == l
        I *= J
    γ = I.minimal_element()
    I *= γ.conjugate() / I.norm()
    return I

def random_ideals(p):
    import os, subprocess, threading, time, ast

    num = os.cpu_count()

    mtx = threading.Lock()
    lines = []

    def fun():
        with subprocess.Popen(
                    ['./ideal-classes-sample.py', str(p)],
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

    outf = open(f'/tmp/data-{p}.csv', 'w')

    while True:
        time.sleep(float(.1))
        with mtx:
            while lines:
                line = lines.pop(0)
                outf.write(line + '\n')
                data = ast.literal_eval(line)
                yield data

def main(p, howmany=None):

    results = zip(range(howmany), random_ideals(p))

    vals = []
    for n,minima in results:
#        val = minima[0]     # λ₁
        val = minima[-1]    # λ₄
        vals.append(val)
        print(f'{n:6} {val:9.7f}', end='')
        µ = sum(vals) / len(vals)
        print(f' | µ = {µ:9.7f}', end='')
        if len(vals) >= 2:
            σ = sqrt(sum((x - µ)**2 for x in vals) / (len(vals) - 1))
            print(f' | σ = {σ:9.7f}')
        else:
            print()

    from sage.all import histogram
    histogram(vals, bins=100).save(f'/tmp/plot.png')

if __name__ == '__main__':
    import sys
    p = int(sys.argv[1])
    howmany = int(sys.argv[2]) if len(sys.argv) > 2 else None
    main(p, howmany)

