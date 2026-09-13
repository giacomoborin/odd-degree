#!/usr/bin/env python3
from sage.all import *
proof.all(False)
set_random_seed()

def serializeO(O):
    s = repr([O.discriminant(), list(O.basis_matrix())])
    assert parseO(s) == O
    return s

def parseO(Ostring):
    p, basis = sage_eval(Ostring)
    if not is_prime(p):
        raise ValueError(f'{p:#x} is not prime')
    B = QuaternionAlgebra(p)
    O = B.maximal_order(order_basis=tuple(B(g) for g in basis))
    if p != O.discriminant():
        raise ValueError(f'{O} is not maximal')
    return O

def random_order(B):
    return B.maximal_order().random_ideal().right_order()

def two_sided_ideal(O):
    p = O.discriminant()
    I = O.left_ideal([p, O.basis()[-1]])
    assert I.norm() == p
    return I

def random_ideals(B, O=None):
    print(f'{B = }', file=sys.stderr)
    if O is not None:
        print(f'{O = }', file=sys.stderr)

    while True:
        order = O or random_order(B)
        yield order.random_ideal()

def random_two_sided_ideals(B):
    print(f'{B = }', file=sys.stderr)

    while True:
        order = random_order(B)
        yield two_sided_ideal(order)

def sort_q(q):
    '''Either I'm missing something, or `minkowski_reduction` doesn't
    properly order the vectors
    '''
    mins = [q[i,i] for i in range(4)]
    for i in range(4):
        for j in range(i,4):
            if mins[i] > mins[j]:
                mins[i], mins[j] = mins[j], mins[i]
                q = q.swap_variables(i,j)
    return q
        
def minima(I, normalize=1, gram_mod=None):
    '''Return the 4 reduced minima of the quadratic form of the ideal.

    If `gram_mod` is an integer output the full upper-triangular gram
    matrix modulo `gram_mod`, in column-order.
    '''
    if not (isinstance(normalize, list) or isinstance(normalize, tuple)):
        normalize = [normalize] * 4
    q = sort_q(I.quadratic_form().lll().minkowski_reduction()[0])
    minima = [q[i,i] / normalize[i] for i in range(4)]
    gram = []
    if gram_mod:
        for c in range(4):
            for r in range(c+1):
                gram.append(q[r,c] % gram_mod)
    return minima, gram

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Sample random ideals')
    parser.add_argument('-p', type=ZZ, help='characteristic', default=ZZ(103))
    parser.add_argument('-O', type=parseO, help='A maximal order', default=None)
    parser.add_argument('-2', '--two-sided', action='store_true', default=False,
                        help='Sample two-sided ideals (incompatible with -O)')
    parser.add_argument('-g', '--gram', type=int, default=None)
    args = parser.parse_args()

    if args.O is not None:
        if args.two_sided:
            raise RuntimeError('Cannot use -O and -j together')
        args.p = O.discriminant()
        B = O.quaternionAlgebra()
    else:
        assert args.p > 1
        B = QuaternionAlgebra(args.p)

    # Print forever the minima of the quadratic form of random ideal
    # classes, or random two-sided ideal classes if `-2` is given.
    #
    # Scale down by [p^½, p^½, p^½, p^½] in the first case and
    # [½p^⅓, ½p^⅓, ½p^⅓, ¼p^⅔] in the second.
    #
    # If `-g` is given, also print Gram matrix modulo it
    if args.two_sided:
        ideals = random_two_sided_ideals(B)
        normalize = [args.p**(1/3)/2, args.p**(1/3)/2, args.p**(1/3)/2, args.p/4]
    else:
        ideals = random_ideals(B, args.O)
        normalize = args.p**(1/2)
    for n,I in enumerate(ideals):
        mins, gram = minima(I, normalize, args.gram)
        s = ','.join(f'{float(m):17.15f}' for m in mins)
        if args.gram:
            s += ',' + ','.join(str(g) for g in gram)
        print(s, flush=True)
