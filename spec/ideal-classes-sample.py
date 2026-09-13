#!/usr/bin/env python3
from sage.all import *
proof.all(False)

def cornacchia(M, d):
    if M < 0:
        return []
    if d == 1:
        if ZZ(M).is_prime():
            try:
                return two_squares(M)
            except ValueError:
                return []
    return []

def inert_prime(bound, d=1):
    while True:
        N = random_prime(bound)
        if kronecker(-d, N) == -1:
            return N

def represent_integer(O0, M):
    B = O0.quaternion_algebra()
    i, j, _ = B.gens()
    p = O0.discriminant()
    if i**2 != -1 or j**2 != -p:
        raise ValueError("represent_integer assumes i^2 = -1 and j^2 = -p")
    if M <= p:
        raise ValueError("represent_integer requires M > p")

    m = max(floor(sqrt(M / (2 * p))), 5)
    for _ in range(m**2):
        z = randint(-m, m)
        t = randint(-m, m)
        M_prime = M - p * (z**2 + t**2)
        xy = cornacchia(M_prime, 1)
        if xy:
            x, y = xy
            γ = x + i*y + j*(z + i*t)
            if γ in O0 and γ.reduced_norm() == M:
                return γ
    return None

def random_sqisign_v1_ideal(O0, l=2, eτ=None, Bτ=None):
    p = O0.discriminant()
    λ = ceil(log(p, 2) / 2)
    eτ = eτ or 2 * λ
    Bτ = Bτ or 2**(λ // 2)
    Nl = l**eτ

    for _ in range(10000):
        Nτ = inert_prime(Bτ)
        if Nτ * Nl > 2 * p:
            break
    else:
        raise ValueError("could not find a large enough inert prime")

    for _ in range(10000):
        γ = represent_integer(O0, Nτ * Nl)
        if γ is not None:
            Iτ = O0 * γ + O0 * Nτ
            assert Iτ.norm() == Nτ
            return Iτ

    raise ValueError("could not represent Nτ * l^eτ")

def sqisignv1_random_ideals(p):
    Quat,(i,j,k) = QuaternionAlgebra(p).objgens()
    print(f'{Quat = }', file=sys.stderr)

    O0 = Quat.maximal_order()
    print(f'{O0 = }', file=sys.stderr)
    assert i in O0 and j in O0 and k in O0

    while True:
        yield random_sqisign_v1_ideal(O0)

from sample import two_sided_ideal

def sqisignv1_random_ideals_janus(p):
    Quat,(i,j,k) = QuaternionAlgebra(p).objgens()
    print(f'{Quat = }', file=sys.stderr)

    O0 = Quat.maximal_order()
    print(f'{O0 = }', file=sys.stderr)
    assert i in O0 and j in O0 and k in O0

    while True:
        yield two_sided_ideal(random_sqisign_v1_ideal(O0).right_order())

def short_random_ideals_janus(p):
    Quat,(i,j,k) = QuaternionAlgebra(p).objgens()
    print(f'{Quat = }', file=sys.stderr)

    O0 = Quat.maximal_order()
    print(f'{O0 = }', file=sys.stderr)
    assert i in O0 and j in O0 and k in O0

    while True:
        I = short_random_ideal(O0)
        print(ZZ(I.norm()).factor())
        yield two_sided_ideal(I.right_order())

def random_ideal(O0):
    p = O0.discriminant()
    l = 257   # or something
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
        assert J.left_order() == O
        assert J.norm() == l
        I *= J
    γ = I.minimal_element()
    I *= γ.conjugate() / I.norm()
    return I


def short_random_ideal(O0):
    p = O0.discriminant()
    l = 2   # or something
    steps = ceil(log(p, l) * 0.4)  #TODO
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
    γ = I.minimal_element()
    I *= γ.conjugate() / I.norm()
    return I


def short_random_ideals(p):
    Quat,(i,j,k) = QuaternionAlgebra(p).objgens()
    print(f'{Quat = }', file=sys.stderr)

    O0 = Quat.maximal_order()
    print(f'{O0 = }', file=sys.stderr)
    assert i in O0 and j in O0 and k in O0

    while True:
        yield short_random_ideal(O0)

def random_ideals(p):
    Quat,(i,j,k) = QuaternionAlgebra(p).objgens()
    print(f'{Quat = }', file=sys.stderr)

    O0 = Quat.maximal_order()
    print(f'{O0 = }', file=sys.stderr)
    assert i in O0 and j in O0 and k in O0

    while True:
        yield random_ideal(O0)

def random_bi_ideals(p):
    Quat,(i,j,k) = QuaternionAlgebra(p).objgens()
    print(f'{Quat = }', file=sys.stderr)

    O0 = Quat.maximal_order()
    print(f'{O0 = }', file=sys.stderr)
    assert i in O0 and j in O0 and k in O0

    while True:
        yield two_sided_ideal(random_ideal(O0).right_order())


from math import log

def is_B_smooth(n, B):
    if n < 2:
        return True
    return factor(ZZ(n), limit=B)[-1][0] <= B


if __name__ == '__main__':
    p = ZZ(sys.argv[1])
    if not is_prime(p):
        raise ValueError(f'{p:#x} is not prime')
    buff = 0
    primes = 0
    smooth = 0
    sample = 0
    B = floor(p ** (1/16)) 
    print(f'# log(p) = {log(p,2):.15f}, , B = {B}', file=sys.stderr)
    try:
        # for n,I in enumerate(random_bi_ideals(p)):
        # for n,I in enumerate(sqisignv1_random_ideals_janus(p)):
        for n,I in enumerate(short_random_ideals(p)):
            # print(f'{I.reduced_basis()}', file=sys.stderr)
            # minima = [g.reduced_norm() / I.norm() / sqrt(p) for g in I.reduced_basis()]
            minima_int = [ZZ(g.reduced_norm() / I.norm()) for g in I.reduced_basis()]
            rb = I.reduced_basis()
            minima = [log(g.reduced_norm() / I.norm(),p) for g in rb]
            r = 20
            for _ in range(20):
                g = sum(choice((-1,1)) * randint(1, r) * v for v in rb[:2])
                _n = ZZ(g.reduced_norm() / I.norm())
                if is_B_smooth( _n , B):
                    smooth += 1
                    print(f'found smooth {_n.factor()}')
                    break

            s = ','.join(f'{float(m):17.15f}' for m in minima)
            print(s, flush=True)
            buff += float(minima[0])
            sample += 1
            # if minima_int[0].is_prime():
            #     primes += 1
            # if minima_int[1].is_prime():
            #     primes += 1
            # if is_B_smooth(minima_int[0], B) or is_B_smooth(minima_int[1], B): smooth += 1
            # print(ZZ(minima_int[2]).factor())
            # if int(minima_int[0]).is_smooth(B) or int(minima_int[1]).is_smooth(B): smooth += 1
    except KeyboardInterrupt:
        print(f'# {sample} samples, average log(minima) = {buff/sample:.15f}', file=sys.stderr)
        print(f'# {primes/sample*50:.2f} % primes found', file=sys.stderr)
        print(f'# {smooth/sample*100:.2f} % smooth found', file=sys.stderr)

