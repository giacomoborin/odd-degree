#!/usr/bin/env sage
# -*- coding: utf-8 -*-
"""
Refinement model for Wesolowski Algorithm 1:
  1. van Oorschot-Wiener under a strict 2^80-byte memory budget;
  2. multiple useful collisions;
  3. a separate HD-isogeny speedup sensitivity table.

This script is _intentionally_ separated from the main cost-model script.
The main script is reserved as the compact baseline comparison.  This script
only studies the two more speculative refinements.

Baseline used here
------------------
Only Damien's M3-C per-isogeny model is used:

    C_step = log2(q) * log2(ell + 1)^2,   q = p^2,
    ell    = 2B/3.

Every estimate includes the expected probabilistic factor 1/P0.

van Oorschot-Wiener model
-------------------------
Let M be the concrete W-paper table-size approximation.  We model the
two-sided claw-search state space as

    S = 2M.

For one useful _golden_ collision, the specialized vOW heuristic is

    T = (2.5/P0) * sqrt(S^3 / w) * C_step,

where w is the number of distinguished-point stored in memory.

For K equally detectable useful collisions, this script uses the expected value of,

    T_K = T / K.

This 1/K improvement is a sensitivity model and it can be improved/substituted by other models.
The result must always be above the birthday cost of finding any collision:

    T_any = (sqrt(pi*S/2)/P0) * C_step.

Thus the implemented cost is

    max(T_any, T/K).

The default collision scenarios are (following Benjamin's remark about the existance of polinomially many useful solutions)

    K = 1, log2(p), log2(p)^2.

These are admittedly modest polylogarithmic examples.  Change MULTICOLLISION_DEGREES to test other assumptions.

Memory records
--------------
The strict memory budget is 2^80 bytes.

A distinguished-point record is modeled as

    (start state, distinguished endpoint, trail length),

with

    record_bits = 2*ceil(log2(S)) + 64.

This is rounded up to bytes.  Auxiliary memory and communication are
not counted, so the model is still in favor of the attacker.

HD-isogeny sensitivity
----------------------
We don't assume here any concrete HD-isogeny cost formula. The HD table is
a pure placeholder: an h-bit speedup divides C_step by
2^h.  The default values are h = 0, 4, 8, 12.

The HD table applies this speedup to the full-table M3-C baseline.  The same
h-bit gain can be subtracted from a vOW result if HD isogenies also accelerate
the vOW step function by the same factor.

Gate conversion
---------------
The exact F_(p^2)-multiplication gate counts are copied from the main script's
Karatsuba-Montgomery conversion.  Keeping them here avoids duplicating the
long arithmetic-circuit code.

Displayed output
----------------
The script prints only the two gate-based result tables:

  * vOW under the strict 2^80-byte budget;
  * HD-isogeny evaluation 

The memory-diagnostics function are kept for auditing, but
main() does not print it.

Run:
    sage sqisign_w_refinements_vow_hd_v13.sage
"""

from functools import lru_cache

from sage.all import *
from scipy.optimize import minimize_scalar


RF = RealField(160)

MEMORY_CAP_LOG2_BYTES = 80
VOW_CONSTANT = RF("2.5")
TRAIL_LENGTH_BITS = 64

# K = log2(p)^d for these exponents d.
MULTICOLLISION_DEGREES = [0, 1, 2]

# Hypothetical constant speedups of the M3-C step cost.
HD_GAIN_BITS = [0, 4, 8, 12]

# Targets contain: (label, prime p, NIST gate reference).
# The F_(p^2)-multiplication gate cost is derived automatically from p.
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
# Weselowski paper table size and success probability model.
# ---------------------------------------------------------------------------
def smoothness_probability(X, B):
    """Return Dickman rho(log(X)/log(B))."""
    X = RF(X)
    B = RF(B)

    if X <= 1 or B <= 1:
        return RF(1)

    return RF(dickman_rho(log(X) / log(B)))


def degree_bound_X(p, B):
    """Return X = B^(1/2) * (p/2)^(1/6)."""
    return RF(B)**(RF(1) / 2) * (RF(p) / RF(2))**(RF(1) / 6)


def table_entries(p, B):
    """Return M(p,B) = X^2 * rho(log(X)/log(B))."""
    X = degree_bound_X(p, B)
    return X**2 * smoothness_probability(X, B)


def success_probability(p, B):
    """Return the success probability P0 of one randomized curve."""
    target_degree = (RF(p) / RF(2))**(RF(1) / 3)
    return smoothness_probability(target_degree, B)


def representative_ell(B):
    """Return ell_eff = 2B/3."""
    return max(RF(2), RF(2) * RF(B) / RF(3))


def m3c_step_cost(p, B):
    """Return Damien M3-C cost per generated isogeny."""
    ell = representative_ell(B)
    log_q = RF(2) * log(RF(p), 2)
    return log_q * log(ell + 1, 2)**2


# ---------------------------------------------------------------------------
# F_(p^2)-multiplication to NAND2-operation conversion.
# ---------------------------------------------------------------------------
KARATSUBA_BASE_BITS = 32


def add_nand(n):
    """NAND2 cost for an n-bit ripple-carry addition/subtraction."""
    n = ZZ(max(1, n))
    return ZZ(9) * n - ZZ(4)


def mux_nand(n):
    """NAND2 cost for an n-bit two-to-one multiplexer."""
    n = ZZ(max(1, n))
    return ZZ(4) * n + ZZ(1)


def schoolbook_mul_nand(n):
    """NAND2 cost for an unsigned n-by-n schoolbook multiplier."""
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
    """Recursive Karatsuba framework, with a 32-bit schoolbook base case."""
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
    """NAND2 cost for one modular addition in F_p."""
    n = ZZ(max(1, n))
    return ZZ(2) * add_nand(n + 1) + mux_nand(n + 1)


def fp_montgomery_mul_nand(n):
    """NAND2 cost for one Montgomery multiplication in F_p."""
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
    NAND2 cost for one F_(p^2) multiplication:

        3 F_p Montgomery multiplications + 4 F_p modular additions.
    """
    n = ZZ(ceil(log(RF(p), 2)))

    return (
        ZZ(3) * fp_montgomery_mul_nand(n)
        + ZZ(4) * fp_mod_add_nand(n)
    )


# ---------------------------------------------------------------------------
# Full-table and vOW costs, in total # of F_(p^2)-multiplication.
# ---------------------------------------------------------------------------
def full_m3c_cost(p, B):
    """Return (M/P0)*C_step for the explicit full-table M3-C model."""
    return (
        table_entries(p, B)
        / success_probability(p, B)
        * m3c_step_cost(p, B)
    )


def state_space_size(p, B):
    """Use S = 2M for the two sided collision search state space."""
    return RF(2) * table_entries(p, B)


def record_bytes(p, B):
    """
    Return bytes per distinguished-point record.

    We store two state identifiers and a 64-bit trail length.
    """
    S = state_space_size(p, B)
    state_bits = ZZ(ceil(log(S, 2)))
    bits = ZZ(2) * state_bits + ZZ(TRAIL_LENGTH_BITS)
    return ZZ(ceil(RF(bits) / RF(8)))


def memory_records(p, B):
    """Return w, the number of records fitting in the 2^80-byte budget."""
    cap_bytes = RF(2)**RF(MEMORY_CAP_LOG2_BYTES)
    return cap_bytes / RF(record_bytes(p, B))


def effective_collisions(p, degree):
    """Return K = log2(p)^degree."""
    return log(RF(p), 2)**ZZ(degree)


def vow_multicollision_cost(p, B, degree):
    """
    Return the budgeted vOW cost for K = log2(p)^degree useful collisions.

    The single-golden-collision formula is divided by K and then bounded 
    by the birthday cost of finding any collision.
    """
    P0 = success_probability(p, B)
    S = state_space_size(p, B)
    w = memory_records(p, B)
    C = m3c_step_cost(p, B)
    K = effective_collisions(p, degree)

    golden = VOW_CONSTANT * sqrt(S**3 / w) * C / (P0 * K)
    any_collision = sqrt(pi * S / RF(2)) * C / P0

    return max(golden, any_collision)


# ---------------------------------------------------------------------------
# One-dimensional optimization over B.
# ---------------------------------------------------------------------------
def optimize_B(p, cost_function, extra_parameter=None):
    """Minimize log2(cost) over 2 <= B <= (p/2)^(1/6)."""
    p = RF(p)
    lower = RF(1)
    upper = log((p / RF(2))**(RF(1) / 6), 2)

    def objective(log2_B):
        B = RF(2)**RF(log2_B)

        if extra_parameter is None:
            cost = cost_function(p, B)
        else:
            cost = cost_function(p, B, extra_parameter)

        return float(log(cost, 2))

    result = minimize_scalar(
        objective,
        bounds=(float(lower), float(upper)),
        method="bounded",
        options={"xatol": 1e-7},
    )

    return RF(2)**RF(result.x)


# ---------------------------------------------------------------------------
# Compact Markdown output.
# ---------------------------------------------------------------------------
def gate_cell(log2_gates, nist, log2_B):
    """Format log2(gates) [NIST margin; log2(B*)]."""
    margin = RF(log2_gates) - RF(nist)

    return "{:.2f} [{:+.2f}; {:.1f}]".format(
        float(log2_gates),
        float(margin),
        float(log2_B),
    )


def print_vow_table():
    """Print vOW results for K=1, log2(p), and log2(p)^2."""
    print("### vOW under a strict 2^80-byte budget")
    print("")
    print("Cells: log2(gates) [margin to NIST; log2(B*)].")
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
        cells = []

        for degree in MULTICOLLISION_DEGREES:
            B = optimize_B(p, vow_multicollision_cost, degree)
            cost = vow_multicollision_cost(p, B, degree)
            log2_gates = log(cost, 2) + log(RF(fp2_gates), 2)

            cells.append(gate_cell(log2_gates, nist, log(B, 2)))

        print("| {} | {} | {} |".format(
            name,
            nist,
            " | ".join(cells),
        ))

    print("")


def print_memory_diagnostics():
    """Print record sizes and memory capacity at the K=1 optimum."""
    print("### vOW memory diagnostics at the K=1 optimum")
    print("")
    print("| Target | bytes/record | log2(w records) | log2(S) |")
    print("|:---|---:|---:|---:|")

    for name, p, _ in TARGETS:
        B = optimize_B(p, vow_multicollision_cost, 0)

        print("| {} | {} | {:.2f} | {:.2f} |".format(
            name,
            record_bytes(p, B),
            float(log(memory_records(p, B), 2)),
            float(log(state_space_size(p, B), 2)),
        ))

    print("")


def print_hd_table():
    """Print hypothetical HD speedups applied to full-table M3-C."""
    print("### HD-isogeny sensitivity on full-table M3-C")
    print("")
    print("An h-bit gain divides the M3-C step cost by 2^h.")
    print("Cells: log2(gates) [margin to NIST].")
    print("")

    labels = ["h={}".format(h) for h in HD_GAIN_BITS]
    print("| Target | NIST | " + " | ".join(labels) + " |")
    print("|:---|---:" + "|---:" * len(labels) + "|")

    for name, p, nist in TARGETS:
        fp2_gates = fp2_mul_nand(p)
        B = optimize_B(p, full_m3c_cost)
        base = log(full_m3c_cost(p, B), 2) + log(RF(fp2_gates), 2)

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
    print("## Wesolowski refinements: vOW multi-collisions and HD sensitivity - v13")
    print("")
    print("Baseline: M3-C, ell=2B/3, every model includes 1/P0.")
    print("vOW: S=2M, constant 2.5, dynamic distinguished point records.")
    print("Multi-collision speedup 1/K is a sensitivity assumption.")
    print("HD gains are hypothetical constant speedups.")
    print("")

    # The memory-diagnostics routine remains implemented above for auditing
    # and debugging, but the circulated output displays only gate tables.
    print_vow_table()
    print_hd_table()


if __name__ == "__main__":
    main()
