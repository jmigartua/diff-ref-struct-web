# Example: a synthetic heating run (periclase + corundum, 300 K to 1300 K)

21 powder diffraction patterns of an in-situ heating run, one every 50 K from 300 K to 1300 K, made for trying
diff-ref-struct Web without data of one's own. **The counts are simulated, not measured.** Each pattern is the
calculated pattern of a FullProf control file, computed by the diff-ref-struct core (`drs simulate`), with Poisson
noise added. No real sample, instrument or reaction lies behind them: periclase does not turn into corundum; the
numbers are chosen to give a plausible-looking series in which a second phase appears on heating.

Licence: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Credit: diff-ref-struct (J. M. Igartua, UPV/EHU).

## Files

- `example_NN_TK.xys`: scan NN (01 to 21, in temperature order) at T kelvin, FullProf's XYDATA format (Ins = 10):
  six header lines, then 2θ (degrees), counts, σ = √counts. 1201 points, 2θ 15° to 75° in steps of 0.05°.
  The header lists both phases in every file; below 1000 K corundum's scale is zero.
- `manifest.json`: the files in order, with size, temperature and seed, and the `drs` version that wrote them.

The temperature is in each file's name (`_300K`), which is where the core reads it: `drs plots propose` offers
one set, "21 scans by temperature, 300 to 1300 K", with no file left out.

## The model

Laboratory X-rays, Cu Kα1 only, λ = 1.5406 Å; no monochromator polarisation, no absorption, no zero shift.

| | periclase MgO | corundum α-Al₂O₃ |
| --- | --- | --- |
| space group | F m -3 m (225) | R -3 c (167), hexagonal axes |
| cell at 300 K | a = 4.2112 Å | a = 4.7589 Å, c = 12.991 Å |
| atoms | Mg 0 0 0 (4a); O ½ ½ ½ (4b) | Al 0 0 0.3523 (12c); O 0.3064 0 ¼ (18e) |
| B_iso (Å²), 300 K → 1300 K | Mg 0.30 → 0.45; O 0.35 → 0.50 | Al 0.30 → 0.45; O 0.35 → 0.50 |
| thermal expansion | a(T) = a₃₀₀ (1 + 1.3·10⁻⁵ K⁻¹ (T − 300 K)) | a and c alike, 8·10⁻⁶ K⁻¹ |

The structures are the standard literature values of both minerals, rounded; they were not taken from one cited
refinement and not checked against one here. The expansion coefficients are representative constants, not fits to
measured data. B_iso rises linearly with T.

**Phase fractions.** Corundum's weight fraction is zero up to 1000 K and then follows a logistic in T centred at
1100 K (width 40 K), rescaled to start at zero at 1000 K and reach 0.5 at 1300 K: 0.08 at 1050 K, 0.23 at 1100 K,
0.38 at 1150 K, 0.46 at 1200 K, 0.49 at 1250 K, 0.50 at 1300 K. Periclase holds the rest. The scale factors follow
from Hill and Howard's rule, W ∝ S·Z·M·V, with periclase's scale 0.014 at 300 K when it is the whole sample.

**Profile and background.** Pseudo-Voigt (FullProf's Npr = 5), η = 0.5, U = 0.020, V = −0.010, W = 0.020 deg²
(FWHM about 0.14° to 0.15°), no asymmetry, no preferred orientation. Background: FullProf's 6-term polynomial about
2θ = 45°: 177 counts at 15°, 160 at 30°, 150 at 45°, 147 at 60°, 150 at 75°. The strongest peak, periclase 200, has
about 5300 counts at 300 K and 2700 at 1300 K; corundum's strongest peaks reach about 800 counts at 1300 K.

**Noise.** Poisson, from `drs simulate`'s generator with seed n for scan n (1 to 21).

## Regenerating

```sh
python3 web/examples/tools/make_example.py            # from the repository's root
DRS=/path/to/drs python3 web/examples/tools/make_example.py --pcr-dir /tmp/pcr   # another drs; keep the .pcr files
```

The script (Python 3 standard library) writes one FullProf .pcr per temperature into a temporary folder, runs
`drs simulate example_NN_TK.pcr --seed NN --out example_NN_TK.xys` there, and copies the patterns and the manifest
here. Re-running it with the same `drs` version reproduces every file byte for byte; another version may change the
calculated pattern or the generator, and says its version in each file's title line.
