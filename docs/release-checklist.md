# Public-release checklist

- [ ] Choose a license and confirm that the original source scripts may be redistributed.
- [ ] Confirm permission to publish any research-derived tables or figures; this release excludes the legacy result collection.
- [ ] Run the tests on a supported Python/MDAnalysis environment.
- [ ] Test the CLI with the intended topology/trajectory pair and verify selections, box dimensions, cutoff units, frame sampling, and reported time.
- [ ] Compare a small sample of frame-level distances against an independent calculation.
- [ ] Ensure `.gitignore` excludes production trajectories and generated result directories.
