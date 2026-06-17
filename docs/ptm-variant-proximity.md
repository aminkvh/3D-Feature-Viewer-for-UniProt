# PTM–Variant Proximity (manuscript text)

A short descriptive paragraph for the main text or Supplementary Methods. This is a
look-up/visualization feature, not a statistical algorithm, so the text deliberately makes no
predictive claim.

---

## Main-text version (concise)

**PTM–Variant Proximity.** To help users see whether a post-translational modification (PTM) site
coincides with genetic variation in three dimensions, the viewer reports, for each annotated PTM
residue, the variants that fall on or near it in the loaded structure. Variants are grouped by
Cα–Cα distance to the PTM into three descriptive tiers: on the same residue, within 8 Å, or within
12 Å. For each PTM the viewer lists the matching variants with their distances, a count of how
many lie within 8 Å (and how many of those are pathogenic or predicted-deleterious), and the
single nearest variant. Selecting a PTM optionally draws dashed lines from it to the nearest and
pathogenic-proximal variants. This is a structural look-up that surfaces spatial co-occurrence for
hypothesis generation; it does not predict PTM–variant crosstalk, regulatory interaction, or
pathogenicity, and the tiers reflect distance only, not functional evidence.

---

## Supplementary Methods version (with the few specifics)

**PTM–Variant Proximity.** This feature summarises the spatial relationship between
post-translational-modification (PTM) sites and genetic variants in the structure currently loaded
in the viewer. For every PTM site that is resolved in the model, the Cα–Cα distance to each
variant-bearing residue is computed and the variant is assigned to one of three distance tiers:
Tier 1, a variant on the same residue as the PTM; Tier 2, a pathogenic or predicted-deleterious
variant within 8 Å; and Tier 3, any variant within 8–12 Å (or a non-pathogenic variant on the same
residue). For each PTM site the interface reports the per-tier list of variants and their
distances, the number of variants within 8 Å and the number of those that are
pathogenic/deleterious, and the closest variant. When a PTM site is selected, dashed lines are
optionally drawn from the PTM Cα to the Cα of its Tier 1 and Tier 2 variants. PTM and variant
annotations are the same UniProt/ClinVar-derived records used elsewhere in the tool; pathogenicity
labels follow the consequence categories defined in the variant view.

The 8 Å and 12 Å radii are descriptive cut-offs chosen so that Tier 2 approximates direct
side-chain contact range and Tier 3 a wider structural neighbourhood; they are not derived from a
model or fitted to data. The feature performs no statistical test and makes no prediction: it is a
distance-based look-up intended to flag PTM sites that warrant closer inspection because
disease-associated or predicted-deleterious variation occurs on or near them in 3D. Spatial
co-occurrence reported here should not be interpreted as evidence of PTM–variant crosstalk or of a
mechanistic or clinical relationship.
