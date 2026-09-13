#!/usr/bin/env sage
# -*- coding: utf-8 -*-
"""
Compact expected-cost model for Algorithm 1 in Wesolowski's attack.

Purpose
-------
This script compares a small set of  cost models for constructing
the meet-in-the-middle table in Algorithm 1.  

Definition of a table entry
---------------------------
In this script, one _entry_ means one newly generated non-backtracking
isogeny or, equivalently one accepted root/edge produced by expanding 
a given source curve with j-invariant.

Processing the modular polynomial Phi_ell(j, X) yields ell+1 roots.  
After excluding the backtracking edge, it produces approximately ell new entries.  
Hence:

  * M3-A-amortized and M3-B-amortized divide their per-polynomial costs by ell;
  * M3-A-raw and M3-B-raw deliberately do not amortize and are shown the fact that computing the roots
    of the modular polynomials has an non-negigible memory cost and also computational cost;
  * M3-C is Damien's per-isogeny, batching-based expression and
    must not be divided by ell again;
  * M0 and M1 are defined directly per table entry/isogeny.

Every model includes the expected probability factor 1/P0. There is no P0=1 case and no extra batching acceleration.

Common expected cost formula for all the models
-----------------------------------------------
Every model uses

    T(p, B) = (M(p, B) / P0(p, B)) * C(p, B),

where

    p         is the prime being evaluated;
    B         is the smoothness bound, optimized independently per model;
    M(p, B)   is the Wesolowski's paper concrete table size approximation;
    P0(p, B)  is the Dickman-rho success estimate for one randomized curve;
    C(p, B)   is the model-dependent cost per generated table entry.

Main models
-----------
M0 -- Wesolowski's paper expected baseline:
    C = 1.

    For gate conversion only, this unit is normalized to one complete F_(p^2)
    multiplication.

M1 -- additional cost per B-isogeny:
    C = log2(B).

M3-A-amortized:
    Per-polynomial cost: ell^2 * log2(p).
    Per-entry cost after division by ell:
    C = ell * log2(p).

M3-B-amortized:
    Per-polynomial cost: ell * log2(p).
    Per-entry cost after division by ell:
    C = log2(p).

M3-C -- D's batching-based baseline:
    C = log2(q) * log2(ell + 1)^2, where q = p^2.

    This expression is already stated per isogeny.  No further amortization
    or batching divisor is applied.

Balanced M3-C -- _experimental_ balanced-factorization refinement:
    B = p^(1/12), C = 2, and alpha_C = 1/2.

    Let D = (p/2)^(1/3).  The balance event assumes that a B-smooth degree
    d <= D admits d = d1*d2 with max(d1,d2)/min(d1,d2) < C.  Under this
    event both halves fit in a single table with cutoff

        X_bal = sqrt(C*D).

    The expected success probability is

        P_bal = P0(p,B) * alpha_C,

    where alpha_C is conditional on B-smoothness.  The per-entry arithmetic
    cost remains M3-C.  This column is experimental and is not optimized over
    B, C, or alpha_C.

Pessimistic cases
---------------------------------
M3-A-raw:
    C = ell^2 * log2(p) per entry.

M3-B-raw:
    C = ell * log2(p) per entry.

These two raw models intentionally consider the complete per-polynomial cost
to every generated entry.  They may be seen as conservative estimates for
omitted modular polynomial work, memory traffic, table management, hidden
constants, imperfect batching, etc.  

Note that M3-A-amortized is identical to M3-B-raw.  Identical
results for those two classes are expected, not a bug.

Representative ell
------------------
The M3 models use

    ell_eff = 2B/3.

This is the output-weighted prime heuristic

    ell_eff ~= (sum_{prime ell <= B} ell^2)
               / (sum_{prime ell <= B} ell)
            ~ 2B/3.

It is an approximation, by no means a precise average over Algorithm 1.
This of course can be further improved.

Optimization
------------
The original models optimize B independently over

    2 <= B <= (p/2)^(1/6).

Optimization is performed in log2(B), using SciPy.  The balanced M3-C
column instead uses the fixed experimental choice B = p^(1/12), C = 2,
and alpha_C = 1/2.  

Gate conversion
---------------
For comparison with NIST's classical gate references, one F_(p^2)
multiplication is modeled as

    3 Montgomery multiplications in F_p
    + 4 modular additions in F_p.

Each F_p Montgomery multiplication includes integer multiplication and
Montgomery reduction.  Integer multiplication uses recursive Karatsuba down
to a 32-bit schoolbook multiplier.

The resulting two-input NAND counts are circuit size models.  They are not
synthesized ASIC areas, nor timing estimates, nor a NIST-prescribed field operation conversion.

Costs not considered
---------------------
Unless represented by a model formula, the script explicitly omits:

  * modular polynomial construction or precomputation;
  * memory capacity and memory traffic;
  * sorting, hashing, collision management, and communication;
  * control, routing, fan-out, registers, and synchronization;
  * hidden constants in asymptotic root-finding algorithms;
  * Specialized prime implementation optimizations.

Output
------
The output table places raw and amortized M3-A/M3-B side by side in a
single NIST gate-comparison table.  The operation count and arithmetic
conversion table functions remain in the source for auditing, but main()
does not print them.

Operation-table cells:
    log2(T) (log2(B*)).

Gate-table cells:
    log2(NAND2-operation equivalent circuit size) [margin to achieve NIST].

The Markdown tables are printed in MD format.

Dependencies
------------
    SageMath
    SciPy

Run
---
    sage W_attack_models_v13.sage
"""

from functools import lru_cache

from sage.all import *
from scipy.optimize import minimize_scalar


RF = RealField(160)
KARATSUBA_BASE_BITS = 32

# Experimental balanced-factorization parameters.
#
# alpha_C is the conditional probability
#   Pr[balanced ratio < C | degree is B-smooth].
BALANCED_B_EXPONENT = RF(1) / RF(12)
BALANCED_C = RF(2)
BALANCED_ALPHA = RF(1) / RF(2)

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


DISPLAY_MODELS = [
    ("M0", "M0"),
    ("M1", "M1"),
    ("A raw", "M3-A-RAW"),
    ("A amort.", "M3-A-AMORT"),
    ("B raw", "M3-B-RAW"),
    ("B amort.", "M3-B-AMORT"),
    ("M3-C", "M3-C"),
    ("Bal. M3-C", "BAL-M3-C"),
]

ALL_MODELS = DISPLAY_MODELS


# ---------------------------------------------------------------------------
# Wesolowski table size and expected success probability.
# ---------------------------------------------------------------------------
def smoothness_probability(X, B):
    """
    Estimate the probability that an integer of size X is B-smooth.
    The estimate is Dickman rho(log(X)/log(B)).
    """
    X = RF(X)
    B = RF(B)

    if X <= 1 or B <= 1:
        return RF(1)

    return RF(dickman_rho(log(X) / log(B)))


def degree_bound_X(p, B):
    """
    Return X = B^(1/2) * (p/2)^(1/6), as in Algorithm 1.
    """
    p = RF(p)
    B = RF(B)
    return B**(RF(1) / 2) * (p / RF(2))**(RF(1) / 6)


def table_entries(p, B):
    """
    Return the concrete W-paper approximation to the number of table entries:
        M(p,B) = X^2 * rho(log(X)/log(B)).
    """
    X = degree_bound_X(p, B)
    return X**2 * smoothness_probability(X, B)


def success_probability(p, B):
    """
    Estimate P0, the success probability of one randomized-curve attempt.
    Every displayed model includes the expected repetition factor 1/P0.
    """
    target_degree = (RF(p) / RF(2))**(RF(1) / 3)
    return smoothness_probability(target_degree, B)


def representative_ell(B):
    """
    Return ell_eff = 2B/3, the output-weighted representative prime.
    """
    return max(RF(2), RF(2) * RF(B) / RF(3))


# ---------------------------------------------------------------------------
# Experimental balanced-factorization refinement.
# ---------------------------------------------------------------------------
def balanced_smoothness_bound(p):
    """
    Return the fixed experimental choice B = p^(1/12).
    """
    return RF(p)**BALANCED_B_EXPONENT


def balanced_degree_cutoff(p):
    """
    Return X_bal = sqrt(C*D), where D = (p/2)^(1/3).

    If a B-smooth degree d <= D factors as d=d1*d2 with d2/d1 < C,
    then both d1 and d2 are at most X_bal.
    """
    D = (RF(p) / RF(2))**(RF(1) / RF(3))
    return sqrt(BALANCED_C * D)


def balanced_table_entries(p):
    """
    Return the balanced table-size approximation

        M_bal = X_bal^2 * rho(log(X_bal)/log(B)).

    This uses the same table-size model as the original columns, but with
    the reduced cutoff X_bal = sqrt(C*D).
    """
    B = balanced_smoothness_bound(p)
    X = balanced_degree_cutoff(p)
    return X**2 * smoothness_probability(X, B)


def balanced_success_probability(p):
    """
    Return P_bal = P0(p,B) * alpha_C.

    alpha_C is interpreted as the conditional probability that a B-smooth
    minimal degree admits a factorization with ratio below C.
    """
    B = balanced_smoothness_bound(p)
    return success_probability(p, B) * BALANCED_ALPHA


def balanced_expected_cost_fp2_units(p):
    """
    Return the expected balanced M3-C work in normalized F_(p^2) units:

        T_bal = M_bal / (P0 * alpha_C) * C_M3-C(p,B).

    The M3-C factor is already a per-isogeny cost and is not amortized again.
    """
    B = balanced_smoothness_bound(p)
    probability = balanced_success_probability(p)

    if probability <= 0:
        return Infinity

    return (
        balanced_table_entries(p)
        / probability
        * entry_cost_fp2_units(p, B, "M3-C")
    )


# ---------------------------------------------------------------------------
# Per-entry costs in normalized F_(p^2) operation units.
# ---------------------------------------------------------------------------
def entry_cost_fp2_units(p, B, model):
    """
    Return the model-dependent cost C(p,B) per generated table entry.

    M3-A-RAW and M3-B-RAW charge a complete per-polynomial cost per entry.
    M3-A-AMORT and M3-B-AMORT divide those costs by approximately ell.
    M3-C is already a per-isogeny expression.
    """
    p = RF(p)
    B = RF(B)

    log_p = log(p, 2)
    log_q = RF(2) * log_p
    ell = representative_ell(B)

    if model == "M0":
        return RF(1)

    if model == "M1":
        return max(RF(1), log(B, 2))

    if model == "M3-A-RAW":
        return ell**2 * log_p

    if model == "M3-A-AMORT":
        return (ell**2 * log_p) / ell

    if model == "M3-B-RAW":
        return ell * log_p

    if model == "M3-B-AMORT":
        return (ell * log_p) / ell

    if model == "M3-C":
        return log_q * log(ell + 1, 2)**2

    raise ValueError("unknown model: {}".format(model))


def expected_cost_fp2_units(p, B, model):
    """
    Return expected work:
        T = (M(p,B) / P0(p,B)) * C(p,B).
    """
    probability = success_probability(p, B)

    if probability <= 0:
        return Infinity

    return (
        table_entries(p, B)
        / probability
        * entry_cost_fp2_units(p, B, model)
    )


# ---------------------------------------------------------------------------
# Optimization over B.
# ---------------------------------------------------------------------------
def optimize_B(p, model):
    """
    Minimize log2(expected cost) over 2 <= B <= (p/2)^(1/6).
    The search variable is log2(B).  
    """
    p = RF(p)
    lower = RF(1)  # B = 2
    upper = log((p / RF(2))**(RF(1) / 6), 2)

    def objective(log2_B):
        B = RF(2)**RF(log2_B)
        cost = expected_cost_fp2_units(p, B, model)

        if cost <= 0 or cost == Infinity:
            return float("inf")

        return float(log(cost, 2))

    result = minimize_scalar(
        objective,
        bounds=(float(lower), float(upper)),
        method="bounded",
        options={"xatol": 1e-7},
    )

    return RF(2)**RF(result.x)


# ---------------------------------------------------------------------------
# NAND2-operation-equivalent arithmetic model.
# ---------------------------------------------------------------------------
def add_nand(n):
    """
    NAND2 cost for n-bit ripple-carry addition or subtraction.
    """
    n = ZZ(max(1, n))
    return ZZ(9) * n - ZZ(4)


def mux_nand(n):
    """
    NAND2 cost for an n-bit two-to-one multiplexer.
    """
    n = ZZ(max(1, n))
    return ZZ(4) * n + ZZ(1)


def schoolbook_mul_nand(n):
    """
    NAND2 cost for an unsigned n-by-n schoolbook multiplier.
    """
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
    """
    NAND2 cost for recursive n-by-n Karatsuba integer multiplication.
    """
    n = ZZ(max(1, n))

    if n <= KARATSUBA_BASE_BITS:
        return schoolbook_mul_nand(n)

    half = (n + 1) // 2

    recursive_products = (
        ZZ(2) * karatsuba_mul_nand(half)
        + karatsuba_mul_nand(half + 1)
    )

    additions = (
        ZZ(2) * add_nand(half + 1)
        + ZZ(2) * add_nand(2 * half + 2)
        + ZZ(2) * add_nand(2 * n)
    )

    return recursive_products + additions


def fp_mod_add_nand(n):
    """
    NAND2 cost for one modular addition in F_p.
    """
    n = ZZ(max(1, n))
    return ZZ(2) * add_nand(n + 1) + mux_nand(n + 1)


def fp_montgomery_mul_nand(n):
    """
    NAND2 cost for one Montgomery multiplication in F_p.

    The model assumes three integer-multiplication terms, accumulation,
    and final conditional reduction.
    """
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
# Computation and compact Markdown output.
# ---------------------------------------------------------------------------
def operation_cell(log2_cost, log2_B):
    """
    Format a field-operation cell as log2(T) (log2(B*)).
    """
    return "{:.2f} ({:.1f})".format(float(log2_cost), float(log2_B))


def gate_cell(log2_gates, nist_reference):
    """
    Format a gate cell as log2(size) [margin to NIST].
    """
    margin = RF(log2_gates) - RF(nist_reference)
    return "{:.2f} [{:+.2f}]".format(float(log2_gates), float(margin))


def compute_rows():
    """
    Optimize B independently for every target and model.
    """
    rows = []

    for name, p, nist_reference in TARGETS:
        model_results = {}

        for _, model in ALL_MODELS:
            if model == "BAL-M3-C":
                B_star = balanced_smoothness_bound(p)
                cost = balanced_expected_cost_fp2_units(p)
            else:
                B_star = optimize_B(p, model)
                cost = expected_cost_fp2_units(p, B_star, model)

            model_results[model] = {
                "log2_B": log(B_star, 2),
                "log2_ops": log(cost, 2),
            }

        rows.append({
            "name": name,
            "p": RF(p),
            "nist": nist_reference,
            "fp2_gate_cost": RF(fp2_mul_nand(p)),
            "models": model_results,
        })

    return rows


def print_operation_table(rows, title, models):
    """
    Print a Zulip/GitHub-compatible expected-operation table.
    """
    print("### {}".format(title))
    print("")
    print("Each entry is log2(T), with log2(B*) in parentheses.")
    print("")

    headings = ["Target", "log2(p)"] + [label for label, _ in models]
    print("| " + " | ".join(headings) + " |")
    print("|:---|---:" + "|---:" * len(models) + "|")

    for row in rows:
        cells = [
            operation_cell(
                row["models"][model]["log2_ops"],
                row["models"][model]["log2_B"],
            )
            for _, model in models
        ]

        print("| {} | {:.1f} | {} |".format(
            row["name"],
            float(log(row["p"], 2)),
            " | ".join(cells),
        ))

    print("")


def print_gate_table(rows, title, models):
    """
    Print a Zulip/GitHub-compatible NIST gate-comparison table.
    """
    print("### {}".format(title))
    print("")
    print("Each entry is log2(NAND2-operation-equivalent size),")
    print("with margin to NIST in brackets.")
    print("")

    headings = ["Target", "NIST ref."] + [label for label, _ in models]
    print("| " + " | ".join(headings) + " |")
    print("|:---|---:" + "|---:" * len(models) + "|")

    for row in rows:
        gate_shift = log(row["fp2_gate_cost"], 2)

        cells = [
            gate_cell(
                row["models"][model]["log2_ops"] + gate_shift,
                row["nist"],
            )
            for _, model in models
        ]

        print("| {} | {} | {} |".format(
            row["name"],
            row["nist"],
            " | ".join(cells),
        ))

    print("")


def print_arithmetic_summary(rows):
    """
    Print the F_(p^2)-multiplication conversion used for each target.
    """
    print("### Arithmetic conversion")
    print("")
    print("| Target | F_(p^2) multiplication (NAND2 operations) | log2 |")
    print("|:---|---:|---:|")

    for row in rows:
        print("| {} | {} | {:.2f} |".format(
            row["name"],
            ZZ(row["fp2_gate_cost"]),
            float(log(row["fp2_gate_cost"], 2)),
        ))

    print("")


def main():
    print("## SQIsign / Wesolowski Algorithm 1 compact cost model - version 13")
    print("")
    print("Entry = one newly generated non-backtracking isogeny/table record.")
    print("Every model includes 1/P0; no P0=1 case is present.")
    print("M3-A/B amortized divide per-polynomial costs by ell.")
    print("M3-C is already per isogeny and is not divided again.")
    print("Raw and amortized M3-A/B are shown side by side.")
    print(
        "Balanced M3-C uses B=p^(1/12), C=2, "
        "and alpha_C=1/2 (experimental)."
    )
    print("")

    rows = compute_rows()

    # The operation-count and F_(p^2)-arithmetic tables remain implemented
    # above for auditing and debugging, but Professor L requested that the
    # circulated output display only the final gate-comparison table.
    print_gate_table(
        rows,
        "NIST classical-gate comparison",
        DISPLAY_MODELS,
    )

    print("Model summary:")
    print("- M0: W-paper one-field-operation-per-entry baseline.")
    print("- M1: log2(B) per generated isogeny.")
    print("- M3-A amort.: ell*log2(p), from ell^2*log2(p) divided by ell.")
    print("- M3-B amort.: log2(p), from ell*log2(p) divided by ell.")
    print("- M3-C: log2(q)*log2(ell+1)^2, q=p^2; already per isogeny.")
    print(
        "- Bal. M3-C: X=sqrt(C*(p/2)^(1/3)), B=p^(1/12), "
        "C=2, alpha_C=1/2."
    )
    print("- M3-A/B raw: unamortized pessimistic per-entry envelopes.")
    print("")
    print("Memory, data movement, modular-polynomial construction, and hidden")
    print("constants remain outside the model unless represented above.")


if __name__ == "__main__":
    main()
