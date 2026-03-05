# nf-enrich

Nextflow DSL2 module for downstream functional interpretation of ChIP peak sets.

Default inputs come from:
- `/nf-chipseeker/chipseeker_output`
- `/nf-idr/idr_output`
- `/nf-peak-consensus/peak_consensus_output`
- `/nf-diffbind/diffbind_output`

Default peak sources:
- `idr`
- `consensus`
- `diffbind`

The module performs:
- per-peak-set GO enrichment (`BP`, `MF`, `CC`)
- per-peak-set KEGG enrichment
- per-source `compareCluster` GO/KEGG comparison
- peak overlap Venn diagram when a source has exactly 2 peak sets

Outputs are written under `enrich_output/functional_enrich_results/`.
