# Software references

## MDAnalysis

This repository uses MDAnalysis for topology/trajectory loading, selections, and distance calculations. The MDAnalysis project asks researchers to cite both foundational publications when using the library in published work:

1. N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. “MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.” *Journal of Computational Chemistry* 32 (2011), 2319–2327. https://doi.org/10.1002/jcc.21787
2. R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein. “MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations.” *Proceedings of the 15th Python in Science Conference* (2016), 98–105. https://doi.org/10.25080/Majora-629e541a-00e

Official citation instructions: https://docs.mdanalysis.org/stable/documentation_pages/references.html

## Distance and periodic-boundary behavior

The trajectory analyzer calls `MDAnalysis.lib.distances.capped_distance` with frame coordinates and, when PBC is enabled and valid unit-cell dimensions exist, the box in MDAnalysis format `[lx, ly, lz, alpha, beta, gamma]`. The minimum-image behavior and accepted box format are documented in the official MDAnalysis distance-analysis API: https://docs.mdanalysis.org/stable/documentation_pages/lib/distances.html

When publishing results from this tool, report the topology/trajectory source, atom selections, cutoff, PBC policy, time-unit/offset interpretation, and frame sampling. Cite the software itself and the original scientific work that generated the analyzed simulation data.
