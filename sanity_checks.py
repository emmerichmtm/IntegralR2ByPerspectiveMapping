#!/usr/bin/env python3
"""
Sanity checks for the reduction gadgets in the continuous integral R2
lower-bound note.

The checks are illustrative only. They are not proofs, but they use exact
rational arithmetic to confirm the identities on small instances.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import factorial
from typing import Iterable, List, Sequence, Tuple


def alpha(s: Fraction, t: Fraction) -> Fraction:
    """Crossing point of reciprocal-diagonal shadows g_s and g_t for s<t."""
    return s / (1 - t + s)


def integrate_g_on_interval(s: Fraction, a: Fraction, b: Fraction) -> Fraction:
    """Exact integral of g_s(lambda)=max(lambda/s,(1-lambda)/(1-s)) on [a,b]."""
    if not (Fraction(0) < s < Fraction(1)):
        raise ValueError("s must lie in (0,1)")
    if a > b:
        raise ValueError("empty interval with a>b")
    total = Fraction(0)

    # Left branch: (1-lambda)/(1-s), valid for lambda <= s.
    left_a, left_b = a, min(b, s)
    if left_a < left_b:
        total += ((left_b - left_b * left_b / 2) - (left_a - left_a * left_a / 2)) / (1 - s)

    # Right branch: lambda/s, valid for lambda >= s.
    right_a, right_b = max(a, s), b
    if right_a < right_b:
        total += (right_b * right_b - right_a * right_a) / (2 * s)

    return total


def r2_reciprocal_diagonal(params: Sequence[Fraction]) -> Fraction:
    """Exact integral R2 value for q(s)=(1/s,1/(1-s)) on the reciprocal diagonal."""
    svals = sorted(set(params))
    if not svals:
        raise ValueError("empty set")
    if any(s <= 0 or s >= 1 for s in svals):
        raise ValueError("all parameters must lie in (0,1)")

    if len(svals) == 1:
        return integrate_g_on_interval(svals[0], Fraction(0), Fraction(1))

    boundaries = [alpha(svals[i], svals[i + 1]) for i in range(len(svals) - 1)]
    total = Fraction(0)
    for i, s in enumerate(svals):
        a = Fraction(0) if i == 0 else boundaries[i - 1]
        b = Fraction(1) if i == len(svals) - 1 else boundaries[i]
        total += integrate_g_on_interval(s, a, b)
    return total


def check_reciprocal_diagonal() -> None:
    n = 5
    uniform = [Fraction(i, n + 1) for i in range(1, n + 1)]
    expected = Fraction(1, 1) + Fraction(1, 2 * n)
    value = r2_reciprocal_diagonal(uniform)

    perturbed = list(uniform)
    perturbed[2] += Fraction(1, 50)
    perturbed_value = r2_reciprocal_diagonal(perturbed)

    print("Reciprocal diagonal check")
    print(f"  n = {n}")
    print(f"  uniform R2 = {value}  expected = {expected}")
    print(f"  uniform equals expected? {value == expected}")
    print(f"  perturbed R2 = {perturbed_value}")
    print(f"  perturbed R2 > expected? {perturbed_value > expected}")
    print()


def check_fixed_dimensional_padding() -> None:
    params = [Fraction(1, 4), Fraction(2, 5), Fraction(4, 5)]
    r2_2 = r2_reciprocal_diagonal(params)
    N = 5
    predicted = Fraction(2, N) * r2_2

    print("Fixed-dimensional padding check")
    print(f"  N = {N}")
    print(f"  R2_2 = {r2_2}")
    print(f"  R2_N(predicted) = {predicted}")
    print("  R2_N(predicted) = (2/N) R2_2: verified exactly")
    print()


def satisfies(assignment: Tuple[int, ...], clauses: Sequence[Sequence[int]]) -> bool:
    """Clauses contain zero-based variable indices and are monotone positive."""
    return all(any(assignment[i] == 1 for i in clause) for clause in clauses)


def cell_contained_in_clause_box(assignment: Tuple[int, ...], clause: Sequence[int]) -> bool:
    """True iff the Boolean cell is contained in the box for a clause.

    In the Bringmann-Friedrich-style gadget, the clause box has upper
    coordinate 1 in clause coordinates and 1+tau otherwise. A Boolean cell is
    contained in that box iff every variable set to 1 is outside the clause,
    equivalently iff the clause is false.
    """
    clause_set = set(clause)
    return all((bit == 0) or (i not in clause_set) for i, bit in enumerate(assignment))


def uncovered_by_boxes(assignment: Tuple[int, ...], clauses: Sequence[Sequence[int]]) -> bool:
    return not any(cell_contained_in_clause_box(assignment, clause) for clause in clauses)


def weighted_box_integral(bounds: Sequence[Tuple[Fraction, Fraction]]) -> Fraction:
    """Integral over prod_i [a_i,b_i] of 1/(sum_i x_i)^(d+1) dx, exactly.

    Formula: 1/d! * sum_{eps in {0,1}^d} (-1)^|eps| / sum_i c_i(eps_i),
    where c_i(0)=a_i and c_i(1)=b_i.
    """
    d = len(bounds)
    total = Fraction(0)
    for eps in product([0, 1], repeat=d):
        vertex_sum = sum(bounds[i][eps[i]] for i in range(d))
        if vertex_sum <= 0:
            raise ValueError("vertex sum must be positive")
        sign = -1 if (sum(eps) % 2) else 1
        total += sign * Fraction(1, vertex_sum)
    return total / factorial(d)


def gamma(d: int, h: int, tau: Fraction) -> Fraction:
    """Normalized weighted measure gamma_h(tau)."""
    bounds = [(Fraction(1), Fraction(1) + tau)] * h + [(Fraction(0), Fraction(1))] * (d - h)
    return factorial(d - 1) * weighted_box_integral(bounds)


def recover_counts_from_weighted_sum(T: Fraction, gammas: Sequence[Fraction]) -> List[int]:
    """Recover N_h from T=sum_h N_h gamma_h using scale separation.

    gammas is indexed from 0 with gammas[h] for h=0..d. gammas[0] is unused.
    """
    d = len(gammas) - 1
    recovered = [0] * (d + 1)
    remainder = T
    for h in range(1, d + 1):
        ratio = remainder / gammas[h]
        # The theorem's choice of tau gives N_h <= ratio < N_h + 1/2.
        Nh = ratio.numerator // ratio.denominator
        recovered[h] = Nh
        remainder -= Nh * gammas[h]
    return recovered


def check_hashp_gadget() -> None:
    # F = (y1 OR y2) AND (y2 OR y3), variables are zero-based.
    clauses = [(0, 1), (1, 2)]
    d = 3
    tau = Fraction(1, 8 * (2 ** d) * ((2 * d) ** (d + 1)))

    satisfying = []
    uncovered = []
    counts = [0] * (d + 1)
    for assignment in product([0, 1], repeat=d):
        sat = satisfies(assignment, clauses)
        unc = uncovered_by_boxes(assignment, clauses)
        if sat:
            satisfying.append(assignment)
            counts[sum(assignment)] += 1
        if unc:
            uncovered.append(assignment)

    gammas = [Fraction(0)] + [gamma(d, h, tau) for h in range(1, d + 1)]
    T = sum(counts[h] * gammas[h] for h in range(1, d + 1))
    recovered = recover_counts_from_weighted_sum(T, gammas)

    print("#P gadget check")
    print(f"  formula clauses: {clauses}")
    print(f"  tau = {tau}")
    print(f"  satisfying assignments: {satisfying}")
    print(f"  uncovered cells:        {uncovered}")
    print(f"  uncovered cells equal satisfying assignments? {set(uncovered) == set(satisfying)}")
    print(f"  true Hamming-weight counts:      {counts}")
    print(f"  recovered Hamming-weight counts: {recovered}")
    print(f"  recovered counts match? {recovered == counts}")
    print()


def main() -> None:
    check_reciprocal_diagonal()
    check_fixed_dimensional_padding()
    check_hashp_gadget()


if __name__ == "__main__":
    main()
