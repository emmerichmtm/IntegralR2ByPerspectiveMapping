# Tchebycheff Shadow Perspective Mapping for Integral R2

This directory contains a standalone LaTeX report and Python verification code for computing the integral R2 indicator by the Tchebycheff shadow perspective mapping.

## Files

- `integral_r2_tchebycheff_perspective_report.tex` -- LaTeX source of the report.
- `integral_r2_tchebycheff_perspective_report.pdf` -- compiled report.
- `integral_r2_perspective.py` -- perspective-map implementation for 3-D absolute integral R2 and anchor-normalized improvement.
- `r2_subdivision3d.py` -- subdivision-based exact evaluator used for verification.
- 'sanity_checls.py` -- some sanity checks related to complexity proofs
- `verification_results.csv` -- reproducibility table comparing perspective, subdivision, and Monte Carlo values.

## Mathematical setting

The ideal point is translated to the origin. Every approximation point is assumed to be a strictly positive loss vector. Equivalently, in the original objective space, the ideal point strictly dominates every approximation point.

The main functions are:

- `perspective_r2_value(points, anchor=None)` computes the absolute unnormalised integral R2 value.
- `perspective_r2_improvement(points, anchor)` computes `R2({anchor}) - R2(points)`.
- `single_point_r2_perspective(point)` computes `R2({point})` by semi-infinite weighted boxes.

The code uses a transparent z-sweep box emitter for reproducibility. The report separates the box emitter from the weighted integration step. Replacing this emitter by an O(n log n), O(n)-box tree-free 3-D hypervolume emitter gives the O(n log n) consequence described in the paper.

## Reproduce the verification

Install dependencies:

```bash
python -m pip install numpy shapely
```

Run:

```bash
python integral_r2_perspective.py
```

This regenerates `verification_results.csv`. The CSV reports both absolute R2 values and anchor-normalized improvement values, compared against a subdivision-based exact evaluator and Monte Carlo estimates.

## Compile the report

```bash
pdflatex -interaction=nonstopmode -halt-on-error integral_r2_tchebycheff_perspective_report.tex
pdflatex -interaction=nonstopmode -halt-on-error integral_r2_tchebycheff_perspective_report.tex
```
