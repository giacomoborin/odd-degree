#!/usr/bin/env sage
# -*- coding: utf-8 -*-
"""
Balanced refinements for Wesolowski Algorithm 1:
  1. van Oorschot-Wiener under a strict 2^80-byte memory budget;
  2. multiple useful collisions;
  3. HD-isogeny speedup sensitivity.

This script uses the experimental balanced M3-C model as its baseline.

Balanced baseline
-----------------
Let

    D       = (p/2)^(1/3),
    B       = p^(1/12),
    C_bal   = 2,
    alpha   = 1/2,
    X_bal   = sqrt(C_bal * D).

The balance assumption is that, conditional on the minimal degree being
B-smooth, it has a factorization d=d1*d2 with ratio below C_bal with
probability alpha.  The success probability of one randomized curve is

    P_bal = P0(p,B) * alpha.

The balanced table-size approximation is

    M_bal = X_bal^2 * rho(log(X_bal)/log(B)).

The per-isogeny arithmetic cost is Damien's M3-C model

    C_step = log2(q) * log2(ell + 1)^2,
    q      = p^2,
    ell    = 2B/3.

No additional division by ell or batching divisor is applied.

van Oorschot-Wiener model
-------------------------
The two-sided collision-search state space is modeled as

    S = 2*M_bal.

For one useful golden collision,

    T_1 = (2.5/P_bal) * sqrt(S^3/w) * C_step,

where w is the number of distinguished-point records fitting in memory.

For K equally detectable useful collisions, the sensitivity model is

    T_K = T_1/K.

The result is bounded below by the generic birthday cost

    T_any = (sqrt(pi*S/2)/P_bal) * C_step.

Hence the implemented vOW cost is

    max(T_any, T_1/K).

The displayed cases are

    K = 1, log2(p), log2(p)^2.

The 1/K rule is a sensitivity assumption, not a theorem for Algorithm 1.

Memory records
--------------
The strict budget is 2^80 bytes.  A distinguished-point record stores

    (start state, distinguished endpoint, trail length),

using

    record_bits = 2*ceil(log2(S)) + 64.

Auxiliary memory and communication are not charged.

HD-isogeny sensitivity
----------------------
The HD table applies hypothetical h-bit speedups to the balanced full-table
M3-C cost.  An h-bit gain divides C_step by 2^h.  The default values are

    h = 0, 4, 8, 12.

These are sensitivity values, not derived HD-isogeny costs.

Gate conversion
---------------
The F_(p^2)-multiplication gate cost is derived from p with the same generic
Karatsuba-Montgomery model used by the compact baseline script.

Displayed output
----------------
Only gate-based tables are printed:

  * balanced vOW under the strict 2^80-byte budget;
  * HD sensitivity on balanced full-table M3-C.

The memory-diagnostics routine remains available for auditing.

Run:
    sage sqisign_w_refinements_vow_hd_v14_balanced.sage
"""

from functools import lru_cache

from sage.all import *


RF = RealField(160)

MEMORY_CAP_LOG2_BYTES = 80
VOW_CONSTANT = RF("2.5")
TRAIL_LENGTH_BITS = 64

BALANCED_B_EXPONENT = RF(1) / RF(12)
BALANCED_C = RF(2)
BALANCED_ALPHA = RF(1) / RF(2)

# K = log2(p)^d for these exponents d.
MULTICOLLISION_DEGREES = [0, 1, 2]

# Hypothetical constant speedups of the balanced M3-C step cost.
HD_GAIN_BITS = [0, 4, 8, 12]

# Targets contain: (label, prime p, NIST gate reference).
TARGETS = [
    ("SQIsign-I (p248.5)",       5  * 2**248 - 1, 143),
    ("Midpoint-II (p306.3)",     3  * 2**306 - 1, 146),
    ("SQIsign-III (p376.65)",   65  * 2**376 - 1, 207),
    ("Midpoint-IV (p458.3)",     3  * 2**458 - 1, 210),
    ("SQIsign-V (p500.27)",     27  * 2**500 - 1, 272),
    ("Mersenne (p521.1)",             2**521 - 1, 272),
    ("Mersenne (p607.1)",             2**607 - 1, 272),
    ("Ultra (p767.43)",          43  * 2**767 - 1, 272),
]


# ---------------------------------------------------------------------------
# Balanced table, success probability, and M3-C step cost.
# ---------------------------------------------------------------------------
def smoothness_probability(X, B):
    """Return Dickman rho(log(X)/log(B))."""
    X = RF(X)
    B = RF(B)

    if X <= 1 or B <= 1:
        return RF(1)

    return RF(dickman_rho(log(X) / log(B)))


def target_degree_bound(p):
    """Return D = (p/2)^(1/3)."""
    return (RF(p) / RF(2))**(RF(1) / RF(3))


def balanced_smoothness_bound(p):
    """Return the fixed experimental choice B = p^(1/12)."""
    return RF(p)**BALANCED_B_EXPONENT


def balanced_degree_cutoff(p):
    """Return X_bal = sqrt(C_bal * D)."""
    return sqrt(BALANCED_C * target_degree_bound(p))


def balanced_table_entries(p):
    """
    Return

        M_bal = X_bal^2 * rho(log(X_bal)/log(B)).
    """
    B = balanced_smoothness_bound(p)
    X = balanced_degree_cutoff(p)
    return X**2 * smoothness_probability(X, B)


def balanced_success_probability(p):
    """
    Return P_bal = P0(p,B) * alpha.

    alpha is conditional on the degree already being B-smooth.
    """
    B = balanced_smoothness_bound(p)
    P0 = smoothness_probability(target_degree_bound(p), B)
    return P0 * BALANCED_ALPHA


def representative_ell(p):
    """Return ell_eff = 2B/3 for the fixed balanced B."""
    B = balanced_smoothness_bound(p)
    return max(RF(2), RF(2) * B / RF(3))


def balanced_m3c_step_cost(p):
    """
    Return M3-C per generated isogeny:

        log2(p^2) * log2(ell+1)^2.
    """
    ell = representative_ell(p)
    log_q = RF(2) * log(RF(p), 2)
    return log_q * log(ell + 1, 2)**2


def balanced_full_m3c_cost(p):
    """Return (M_bal/P_bal)*C_step in F_(p^2)-multiplication units."""
    return (
        balanced_table_entries(p)
        / balanced_success_probability(p)
        * balanced_m3c_step_cost(p)
    )


# ---------------------------------------------------------------------------
# F_(p^2)-multiplication to NAND2-operation conversion.
# ---------------------------------------------------------------------------
KARATSUBA_BASE_BITS = 32


def add_nand(n):
    """NAND2 proxy for an n-bit ripple-carry addition/subtraction."""
    n = ZZ(max(1, n))
    return ZZ(9) * n - ZZ(4)


def mux_nand(n):
    """NAND2 proxy for an n-bit two-to-one multiplexer."""
    n = ZZ(max(1, n))
    return ZZ(4) * n + ZZ(1)


def schoolbook_mul_nand(n):
    """NAND2 proxy for an unsigned n-by-n schoolbook multiplier."""
    n = ZZ(max(1, n))

    if n == 1:
        return ZZ(2)

    return (
        ZZ(2) * n**2
        + ZZ(9) * n * max(ZZ(0), n - 2)
        + ZZ(5) * n
    )


@lru_cache(maxsize=None)
def karatsuba_mul_nand(n):
    """Recursive Karatsuba proxy with a 32-bit schoolbook base case."""
    n = ZZ(max(1, n))

    if n <= KARATSUBA_BASE_BITS:
        return schoolbook_mul_nand(n)

    half = (n + 1) // 2

    return (
        ZZ(2) * karatsuba_mul_nand(half)
        + karatsuba_mul_nand(half + 1)
        + ZZ(2) * add_nand(half + 1)
        + ZZ(2) * add_nand(2 * half + 2)
        + ZZ(2) * add_nand(2 * n)
    )


def fp_mod_add_nand(n):
    """NAND2 proxy for one modular addition in F_p."""
    n = ZZ(max(1, n))
    return ZZ(2) * add_nand(n + 1) + mux_nand(n + 1)


def fp_montgomery_mul_nand(n):
    """NAND2 proxy for one Montgomery multiplication in F_p."""
    n = ZZ(max(1, n))
    integer_mul = karatsuba_mul_nand(n)

    return (
        ZZ(3) * integer_mul
        + add_nand(2 * n)
        + add_nand(n)
        + mux_nand(n)
    )


def fp2_mul_nand(p):
    """
    NAND2 proxy for one F_(p^2) multiplication:

        3 F_p Montgomery multiplications + 4 F_p modular additions.
    """
    n = ZZ(ceil(log(RF(p), 2)))

    return (
        ZZ(3) * fp_montgomery_mul_nand(n)
        + ZZ(4) * fp_mod_add_nand(n)
    )


# ---------------------------------------------------------------------------
# Balanced vOW model.
# ---------------------------------------------------------------------------
def state_space_size(p):
    """Return S = 2*M_bal."""
    return RF(2) * balanced_table_entries(p)


def record_bytes(p):
    """
    Return bytes per distinguished-point record.

    The record stores two state identifiers and a 64-bit trail length.
    """
    S = state_space_size(p)
    state_bits = ZZ(ceil(log(S, 2)))
    bits = ZZ(2) * state_bits + ZZ(TRAIL_LENGTH_BITS)
    return ZZ(ceil(RF(bits) / RF(8)))


def memory_records(p):
    """Return w, the number of records fitting in 2^80 bytes."""
    cap_bytes = RF(2)**RF(MEMORY_CAP_LOG2_BYTES)
    return cap_bytes / RF(record_bytes(p))


def effective_collisions(p, degree):
    """Return K = log2(p)^degree."""
    return log(RF(p), 2)**ZZ(degree)


def balanced_vow_multicollision_cost(p, degree):
    """
    Return the balanced vOW cost for K = log2(p)^degree.

    The golden-collision estimate is divided by K and bounded below by
    the generic birthday cost.
    """
    P_bal = balanced_success_probability(p)
    S = state_space_size(p)
    w = memory_records(p)
    C_step = balanced_m3c_step_cost(p)
    K = effective_collisions(p, degree)

    golden = (
        VOW_CONSTANT
        * sqrt(S**3 / w)
        * C_step
        / (P_bal * K)
    )

    any_collision = (
        sqrt(pi * S / RF(2))
        * C_step
        / P_bal
    )

    return max(golden, any_collision)


# ---------------------------------------------------------------------------
# Compact Markdown output.
# ---------------------------------------------------------------------------
def gate_cell(log2_gates, nist, log2_B):
    """Format log2(gates) [NIST margin; log2(B)]."""
    margin = RF(log2_gates) - RF(nist)

    return "{:.2f} [{:+.2f}; {:.1f}]".format(
        float(log2_gates),
        float(margin),
        float(log2_B),
    )


def print_vow_table():
    """Print balanced vOW results for the configured K scenarios."""
    print("### Balanced vOW under a strict 2^80-byte budget")
    print("")
    print("Cells: log2(gates) [margin to NIST; log2(B)].")
    print("K is the effective number of equally detectable useful collisions.")
    print("")

    labels = [
        "K=1" if d == 0 else ("K=L" if d == 1 else "K=L^{}".format(d))
        for d in MULTICOLLISION_DEGREES
    ]

    print("| Target | NIST | " + " | ".join(labels) + " |")
    print("|:---|---:" + "|---:" * len(labels) + "|")

    for name, p, nist in TARGETS:
        fp2_gates = fp2_mul_nand(p)
        log2_B = log(balanced_smoothness_bound(p), 2)
        cells = []

        for degree in MULTICOLLISION_DEGREES:
            cost = balanced_vow_multicollision_cost(p, degree)
            log2_gates = log(cost, 2) + log(RF(fp2_gates), 2)
            cells.append(gate_cell(log2_gates, nist, log2_B))

        print("| {} | {} | {} |".format(
            name,
            nist,
            " | ".join(cells),
        ))

    print("")


def print_memory_diagnostics():
    """Print balanced record sizes and state-space diagnostics."""
    print("### Balanced vOW memory diagnostics")
    print("")
    print("| Target | bytes/record | log2(w records) | log2(S) |")
    print("|:---|---:|---:|---:|")

    for name, p, _ in TARGETS:
        print("| {} | {} | {:.2f} | {:.2f} |".format(
            name,
            record_bytes(p),
            float(log(memory_records(p), 2)),
            float(log(state_space_size(p), 2)),
        ))

    print("")


def print_hd_table():
    """Print HD speedup sensitivity on balanced full-table M3-C."""
    print("### HD-isogeny sensitivity on balanced full-table M3-C")
    print("")
    print("An h-bit gain divides the balanced M3-C step cost by 2^h.")
    print("Cells: log2(gates) [margin to NIST].")
    print("")

    labels = ["h={}".format(h) for h in HD_GAIN_BITS]
    print("| Target | NIST | " + " | ".join(labels) + " |")
    print("|:---|---:" + "|---:" * len(labels) + "|")

    for name, p, nist in TARGETS:
        fp2_gates = fp2_mul_nand(p)
        base = (
            log(balanced_full_m3c_cost(p), 2)
            + log(RF(fp2_gates), 2)
        )

        cells = []

        for gain in HD_GAIN_BITS:
            value = base - RF(gain)
            margin = value - RF(nist)

            cells.append("{:.2f} [{:+.2f}]".format(
                float(value),
                float(margin),
            ))

        print("| {} | {} | {} |".format(
            name,
            nist,
            " | ".join(cells),
        ))

    print("")


def main():
    print("## Balanced Wesolowski refinements: vOW and HD sensitivity - v14")
    print("")
    print(
        "Balanced baseline: B=p^(1/12), C=2, alpha=1/2, "
        "X=sqrt(C*(p/2)^(1/3))."
    )
    print("Arithmetic baseline: M3-C with ell=2B/3.")
    print("Every estimate includes 1/(P0*alpha).")
    print("vOW uses S=2*M_bal, constant 2.5, and dynamic records.")
    print("The 1/K multi-collision gain and HD gains are sensitivity assumptions.")
    print("")

    # Kept for auditing, but Professor L requested gate tables only.
    print_vow_table()
    print_hd_table()


if __name__ == "__main__":
    main()
