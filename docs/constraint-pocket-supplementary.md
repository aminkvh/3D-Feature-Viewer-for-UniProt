# Constraint-pocket residue prioritization — Methods and Supplementary Material

This document provides (i) a short Methods paragraph suitable for the main text, (ii) a
Supplementary Methods section written in the style of the existing S1–S4 routines, and
(iii) short validation addenda for the S1 (pathogenic-variant enrichment) and S3
(mutation/phenotype burden) routines. It is intended to be merged into the existing
Supplementary Methods document.

---

## Methods (main-text excerpt)

**Constraint-pocket prioritization.** To highlight residues whose evolutionary constraint
exceeds what their structural burial would predict, the viewer computes, for each residue, the
mean AlphaMissense pathogenicity over the nineteen possible substitutions and regresses it
against the residue's coordination number using a local (LOESS) fit; the positive part of the
residual is treated as "unexpected" constraint. A Getis-Ord *G*<sub>*i*</sub><sup>*</sup>
spatial-autocorrelation statistic is applied to these residuals over a Gaussian,
PAE-gated neighbourhood, and significance is assessed by a within-protein permutation null with
Benjamini-Hochberg false-discovery-rate control. The routine runs client-side in a Web Worker
and is exposed as a coloring mode with a user-adjustable sensitivity (FDR) threshold. It is an
exploratory heuristic for locating candidate buried functional sites and is not a validated
predictor of pathogenicity or of protein-protein interfaces.

---

## Supplementary Methods Sx. Constraint-pocket residue prioritization

This routine identifies residues that are more evolutionarily constrained than their structural
burial would predict and that cluster together in three dimensions — a pattern characteristic of
buried catalytic and ligand-binding sites. It decouples functional constraint from the
generic stabilizing constraint of the hydrophobic core, a known confounder in one-dimensional
variant-effect maps (the residual-against-burial idea follows established structural-evolution
analyses, e.g. Sivley et al. 2018). The specific burial normalization, spatial statistic,
permutation null, and tier thresholds below are particular to this implementation.

### Inputs and per-residue constraint

For each residue *i* of the queried protein, raw pathogenicity *P*<sub>*i*</sub> is the mean of
the AlphaMissense scores of the nineteen possible substitutions at position *i* (Cheng et al.
2023). Residues without AlphaMissense scores (including non-human entries, for which
AlphaMissense is not currently computed) are excluded; the analysis is not run when fewer than
twelve scored residues remain.

### Burial normalization (residual pathogenicity)

Structural burial is summarized by the coordination number *c*<sub>*i*</sub>, the number of
Cα atoms within 13 Å of the residue's Cα. Coordination is computed against **every** modelled
Cα in the loaded structure — including additional copies of the queried protein and the chains
of partner proteins in a complex — so that a residue buried at a subunit interface is correctly
treated as buried rather than as solvent-exposed.

The expected pathogenicity at a given burial, *E*(*c*), is estimated by locally weighted linear
regression (LOESS; Cleveland 1979) of *P* on *c* with a tricube kernel and a span of 0.3. The
residual pathogenicity is

> ΔP<sub>*i*</sub> = *P*<sub>*i*</sub> − *E*(*c*<sub>*i*</sub>).

A positive residual indicates a residue that is more intolerant to mutation than other residues
of comparable burial in the same protein. Using a per-protein, continuous LOESS fit (rather than
fixed solvent-accessibility bins) avoids the boundary artifacts and empty/degenerate-bin
failures that affect discretized burial strata in globular proteins.

### Spatial enrichment (Getis-Ord *G*<sub>*i*</sub><sup>*</sup>)

For each ordered pair of scored residues within 13 Å, a spatial weight

> *w*<sub>*ij*</sub> = exp(−*d*<sub>*ij*</sub><sup>2</sup> / 2σ<sup>2</sup>) · exp(−PAE<sub>*ij*</sub> / τ)

is assigned, with Cα–Cα distance *d*<sub>*ij*</sub>, Gaussian bandwidth σ = 5 Å, and
τ = 10 Å. The second factor is the AlphaFold Predicted Aligned Error gate: for AlphaFold models
it down-weights pairs whose relative position is uncertain, suppressing spurious clusters across
flexible linkers; for experimental structures, or when PAE is unavailable, the gate evaluates to
1 and the weighting reduces to the Gaussian distance kernel. The local statistic is the
self-inclusive Getis-Ord *G*<sub>*i*</sub><sup>*</sup> applied to the residuals,

> *G*<sub>*i*</sub><sup>*</sup> ∝ ΔP<sub>*i*</sub> + Σ<sub>*j*</sub> *w*<sub>*ij*</sub> ΔP<sub>*j*</sub>,

so that a residue scores highly when it and its spatial neighbourhood are jointly
constrained-beyond-burial.

### Significance and multiple testing

An 8–13 Å neighbourhood typically contains only a few tens of residues, so analytical *G*<sup>*</sup>
*p*-values (which assume asymptotic normality) are not used. Instead, residuals are randomly
permuted over the fixed residue positions (999 permutations for proteins up to 600 scored
residues, 599 up to 1,200, and 399 above) and the observed *G*<sub>*i*</sub><sup>*</sup> is
compared with its permutation distribution. The one-sided empirical *p*-value is
(*b*<sub>*i*</sub> + 1)/(*m* + 1), where *b*<sub>*i*</sub> is the number of permutations with a
simulated statistic greater than or equal to the observed one and *m* is the number of
permutations; the plus-one estimator prevents zero *p*-values (Phipson and Smyth 2010).
*p*-values are converted to *q*-values within the protein using the Benjamini-Hochberg
false-discovery-rate procedure (Benjamini and Hochberg 1995). A fixed pseudorandom-number-
generator seed is used so that repeated calculations return identical results.

Only residues with a positive residual are considered. For homo-oligomers each physical copy of
a residue is scored against its own environment and the copies are then collapsed to a single
call per UniProt position, retaining the most significant copy. All candidate positions are
returned with their *q*-values; the viewer displays those with *q* ≤ a user-selected threshold
(default 0.10), so the sensitivity slider re-thresholds the display without recomputation. Each
displayed residue is additionally given a geometric label — *buried (pocket)* when it is at or
above the median coordination number or its upward/downward half-sphere ratio exceeds 1.3,
otherwise *exposed* — as a visualization aid, not a validated site-type classification.

### Browser-native implementation

The coordination, spatial-weight, and permutation calculations execute in a Web Worker, leaving
the user-interface thread responsive; if a Worker cannot be created (for example under a
restrictive page content-security policy) the computation falls back to the main thread. The
required data — the AlphaMissense substitution table and the Cα coordinates already loaded for
rendering, plus the AlphaFold PAE matrix fetched on demand — are obtained from public REST
endpoints, and analysis is capped at 3,000 scored residues.

### Validation (small benchmark)

The routine was evaluated on a small set of human enzymes drawn from the Mechanism and Catalytic
Site Atlas (M-CSA; Ribeiro et al. 2018) that have AlphaMissense data, using the annotated
catalytic residues as a positive set and a within-5 Å Cα grace radius. At the default threshold
(*q* ≤ 0.10), catalytic-site recall was enriched by a median of approximately 1.7–1.8× relative
to a same-sized random-residue baseline, with roughly three-quarters of enzymes exceeding chance
and a typical flag rate near 16 % of residues. The burial normalization contributed to this
signal: the median enrichment of the residual-pathogenicity statistic exceeded that of the raw
(unnormalized) AlphaMissense mean clustered in the same way (≈1.7× versus ≈1.3×). Per-protein
enrichment varied widely, which motivates the user-adjustable sensitivity threshold. These
results indicate a modest, exploratory enrichment useful for narrowing candidate active-site
regions; the routine has not been benchmarked as a classifier, makes no diagnostic-performance
claim, and does not detect protein-protein interaction interfaces (which are frequently
mutation-tolerant and therefore carry little AlphaMissense constraint signal).

### Supplementary Table — additional parameters

| Routine | Parameter | Value | Purpose and provenance |
|---|---|---|---|
| Constraint pocket | Coordination radius | 13 Å | Cα neighbourhood for burial and spatial weights. |
| Constraint pocket | LOESS span | 0.3 | Local burial-baseline smoothing (Cleveland 1979). |
| Constraint pocket | Gaussian bandwidth σ | 5 Å | Spatial-weight decay for *G*<sub>*i*</sub><sup>*</sup>. |
| Constraint pocket | PAE decay τ | 10 Å | AlphaFold confidence gate; evaluates to 1 without PAE. |
| Constraint pocket | Permutations | 999 / 599 / 399 | Empirical null, scaled by protein size; plus-one estimator (Phipson and Smyth 2010). |
| Constraint pocket | Default display FDR | *q* ≤ 0.10 | User-adjustable visualization threshold (Benjamini and Hochberg 1995). |
| Constraint pocket | Maximum scored residues | 3,000 | Pragmatic in-browser performance guard. |

### Additional references

- Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with
  AlphaMissense. *Science* 2023;381:eadg7492.
- Cleveland WS. Robust locally weighted regression and smoothing scatterplots. *J Am Stat Assoc*
  1979;74:829–36.
- Getis A, Ord JK. The analysis of spatial association by use of distance statistics. *Geogr Anal*
  1992;24:189–206.
- Ribeiro AJM, Holliday GL, Furnham N, et al. Mechanism and Catalytic Site Atlas (M-CSA): a database
  of enzyme reaction mechanisms and active sites. *Nucleic Acids Res* 2018;46:D618–23.

(Benjamini and Hochberg 1995, Phipson and Smyth 2010, and Sivley et al. 2018 are already cited in
the main Supplementary Methods.)

---

## Validation addenda for S1 and S3

These short paragraphs can be appended to the corresponding routine descriptions. Both
benchmarks are small and exploratory and are reported as supporting evidence, not as
diagnostic-performance evaluations.

### Add to S1 (Pathogenic-variant enrichment hotspots)

> In a small benchmark, hotspot positions were enriched for high-confidence ClinVar variants —
> those reviewed by an expert panel or by multiple submitters without conflict — relative to all
> annotated pathogenic positions (approximately 1.7-fold for TP53 and 11-fold for BRCA1), and
> were correspondingly depleted of low-review-status calls. Because the routine uses only binary
> pathogenic/benign labels and not the review status, this indicates that the spatial clustering
> preferentially recovers well-supported pathogenic positions. The benchmark is exploratory (two
> proteins) and is not a diagnostic-performance evaluation.

### Add to S3 (Mutation/phenotype burden residues)

> In a small benchmark against independent deep-mutational-scanning data (TP53; Giacomelli et al.
> 2018), burden-positive positions showed substantially larger experimental functional effects
> than other variant-containing positions, with roughly half of burden-positive positions falling
> among the strongest-effect fifth of all measured positions (an approximately 2.8-fold
> enrichment). Because deep-mutational-scanning fitness is measured independently of the
> variant-count and phenotype-count inputs used here, this provides external, non-circular
> support. The benchmark is small (one to two proteins) and exploratory.

Additional reference: Giacomelli AO, Yang X, Lintner RE, et al. Mutational processes shape the
landscape of TP53 mutations in human cancer. *Nat Genet* 2018;50:1381–87.
