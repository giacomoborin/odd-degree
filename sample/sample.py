#!sage -python
from sage.all import *
proof.all(False)
set_random_seed()

def serialize(O):
    s = repr([O.discriminant(), list(O.basis_matrix())])
    assert parse(s) == O
    return s

def parse(Ostring):
    p, basis = sage_eval(Ostring)
    if not is_prime(p):
        raise ValueError(f'{p:#x} is not prime')
    B = QuaternionAlgebra(p)
    O = B.maximal_order(order_basis=tuple(B(g) for g in basis))
    if p != O.discriminant():
        raise ValueError(f'{O} is not maximal')
    return O

def random_ideal(O0, deg=None):
    p = O0.discriminant()
    if deg is None:
        deg = p**2
    l = 257   # or something
    steps = ceil(log(deg, l))
    I = O0*1
    for t in reversed(range(steps)):
        O = I.right_order()
        while True:
            α = sum(randrange(l)*g for g in O.basis())
            if α.reduced_norm().valuation(l) == 1:
                break
        assert α in O
        J = O*l + O*α
        assert J.left_order() == O
        assert J.norm() == l
        I *= J
    #γ = I.minimal_element()
    #I *= γ.conjugate() / I.norm()
    return I

def random_ideals(O, p = None):
    if O is None:
        O = random_order(p)
    else:
        print(f'{O = }', file=sys.stderr)
    Obasis = (O * 1).reduced_basis()
    p = O.discriminant()
    Ominima = sorted([g.reduced_norm() / p**(2/3) for g in Obasis])[1:]
    out = ','.join(f'{float(m):7.5f}' for m in Ominima)
    print(f'minima {out = }', file=sys.stderr)

    while True:
        yield random_ideal(O)

def start_order(p):
    A = QuaternionAlgebra(p)
    O0 = A.maximal_order()
    return O0

def random_order(p):
    O0 = start_order(p)
    O = random_ideal(O0).right_order()
    return O

def minima(I, normalize=1):
    basis = I.reduced_basis()
    ORbasis = (I.right_order() * 1).reduced_basis()
    OLbasis = (I.left_order() * 1).reduced_basis()
    minima = sorted([g.reduced_norm() / I.norm() / normalize**(1/2)
                     for g in basis])
    ORminima = sorted([g.reduced_norm() / normalize**(2/3)
                      for g in ORbasis])[1:]
    OLminima = sorted([g.reduced_norm() / normalize**(2/3)
                      for g in OLbasis])[1:]
    inner_products = [[(g*h).reduced_trace() % 2 for h in basis] for g in basis]
    # see if the first 4 are all even
    # first_four_even = all([all([((g*h).reduced_trace() % 2) == 0 for h in basis[:3] ]) for g in basis[:3] ])
    return (minima, OLminima, ORminima, inner_products)


def minima_eichler(I, normalize=1):
    basis = I.reduced_basis()
    minima = sorted([g.reduced_norm() / I.norm() / normalize**(1/2) for g in basis])
    OL = I.left_order()
    E = I.left_order().intersection(I.right_order())
    Eichler_basis = (E * 1).reduced_basis()
    Eichler_minima = sorted([g.reduced_norm() / normalize**(2) for g in Eichler_basis])[1:]
    return (minima, Eichler_minima)



def sample_odd_element(I, normalize=1):
    basis = I.reduced_basis()
    odd_samples = 0.0
    basis_els = 4
    total = 1000
    # for i in range(2**basis_els):
    #     g = sum(((i >> j) & 1)*b for j,b in enumerate(basis))
    #     if g != 0 and g.reduced_trace() % 2 == 1:
    #         odd_samples += 1
    # return odd_samples/(2**basis_els)
    for _ in range(total):
        g = sum((ZZ(randrange(2)))*b for b in basis[:basis_els])
        print(g)
        if g != 0 and g.reduced_trace() % 2 == 1:
            odd_samples += 1
    return odd_samples/total

if __name__ == '__main__':
    # args are either a prime characteristic,
    # or the serialization of a maximal order
    SAVE_OUT = True
    mode = ''
    try:
        p = int(sys.argv[1])
        if p == 248:
            p = 2261564242916331941866620800950935700259179388000792266395655937654553313279
        if p == 62:
            p = 55340232221128654847
        if len(sys.argv) > 2 and sys.argv[2] == '--random':
            O = random_order(p)
            mode = 'random'
        elif len(sys.argv) > 2 and sys.argv[2] == '--birandom':
            O = None
            mode = 'birandom'
        elif len(sys.argv) > 2 and sys.argv[2] == '--eichler':
            O = None
            mode = 'eichler'
        else:
            O = start_order(p)
    except ValueError:
        O = parse(sys.argv[1])
        p = O.discriminant()

    sample_size = 0
    parities = [ [ 0.0 for _ in range(4) ] for _ in range(4) ]
    # Print forever the minima of random left O-ideas classes and
    # their right orders, scaled down by factors of p^½ and p^⅔
    # respectively
    # number_true_ch = 0
    if mode == 'eichler':
        for n,I in enumerate(random_ideals(O, p = p)):
            m, OLm  = minima_eichler(I, p)
            s = ','.join(f'{float(m):17.15f}' for m in m + OLm)
            if SAVE_OUT:
                print(s, flush=True)

    try:
        for n,I in enumerate(random_ideals(O, p = p)):
            m, OLm, ORm, inner = minima(I, p)
            s = ','.join(f'{float(m):17.15f}' for m in m + OLm + ORm)
            s += ',' + ','.join(','.join(str(i) for i in inner_i) for inner_i in inner)
            # s += f',{ch}'
            for i in range(4):
                for j in range(4):
                    parities[i][j] += inner[i][j]
            # if ch:
            #     number_true_ch += 1
            sample_size += 1
            if SAVE_OUT:
                print(s, flush=True)
    except KeyboardInterrupt:
        print(f'{sample_size} samples collected', file=sys.stderr)
        for i in range(4):
            for j in range(4):
                parities[i][j] /= sample_size
        for i in range(4):
            print(' '.join(f'{parities[i][j]:.7f}' for j in range(4)), file=sys.stderr)
        # print(f'number_true_ch = {number_true_ch}, fraction = {number_true_ch/sample_size:.7f}', file=sys.stderr)
