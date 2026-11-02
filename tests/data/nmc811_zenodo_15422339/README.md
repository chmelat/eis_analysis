# EIS data of 54 21700 cells (Zenodo 15422339)

Third-party data, **not covered by the MIT license of this project**.

- **Source:** M. Moertelmaier, M. Kasper (Keysight Technologies, Austria),
  S. Clark (SINTEF, Norway): *EIS data of 54 21700 cells*, Zenodo, 2025,
  https://doi.org/10.5281/zenodo.15422339
- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0),
  https://creativecommons.org/licenses/by/4.0/
- **Changes:** none. The 54 files `ID01.csv`..`ID54.csv` are the
  originals; their MD5 sums match those published in the Zenodo record.
  The record's metadata files (`metadata.jsonld`, `metadata.html`, the
  CSVW schema) are not included.

## Content

As stated by the authors: galvanostatic EIS of a batch of 54 cylindrical
21700 NMC811 lithium-ion cells at 30 % state of charge, measured with a
Keysight SL1007A Scienlab Battery Test System. The manufacturer and the
anode are not stated.

Each file: columns `Frequency_Hz, Real_Ohm, Imag_Ohm` (imaginary part
signed), 44 points from 10 kHz to 50 mHz, |Z| about 14-20 mOhm.

What the spectra show (our reading, not the authors'):

- high frequencies: series inductance ~ 90-96 nH and Re(Z) rising with
  frequency (leads and wound electrode, skin effect);
- middle band: contact resistance, SEI and charge transfer;
- low frequencies: a ~45 degree tail (solid-state diffusion) that
  continues below 50 mHz.

The cells are nominally identical, so differences between their analyses
point at the method rather than the cells. Used by `tests/test_real_nmc811.py`.

## Equivalent circuit (our fits)

Proposed circuit, fitted with DE (eis_analysis v0.64.2, seed 1) to cells
ID11, ID45 and ID48 (cold-agent task T16, `doc/COLD_AGENT.md`, round 8):

    R0 - (R1 | L) - (R2 | Q) - W

```bash
eis ID11.csv --circuit "R(0.0132) - (R(0.039) | L(1e-7)) - (R(0.0045) | Q(0.63, 0.77)) - W(0.00102)"
```

| element | ID11 | ID45 | ID48 | reading |
|---|---|---|---|---|
| R0 [mOhm] | 13.21 | 13.08 | 13.25 | ohmic: electrolyte, contacts, current collectors |
| R1 [mOhm] | 39.4 | 37.9 | 39.0 | R1 \| L: inductance of the cell and leads; R1 makes Re(Z) rise toward f_max (skin effect, eddy currents); corner R1/(2 pi L) ~ 61 kHz, above the window |
| L [nH] | 99.4 | 100.2 | 100.9 | |
| R2 [mOhm] | 5.44 | 4.25 | 4.13 | one depressed arc, tau 0.63 / 0.46 / 0.49 ms (250-350 Hz): SEI and charge transfer merged |
| Q [F s^(n-1)] | 0.71 | 0.59 | 0.59 | |
| n | 0.75 | 0.78 | 0.79 | |
| W sigma [Ohm s^-1/2] | 1.04e-3 | 1.00e-3 | 1.01e-3 | semi-infinite diffusion tail below ~3 Hz |

Fit error 0.32-0.36 %. Against `R0 - L - (R2 | Q) - W` (series L, the
circuit the cold agent chose) BIC is lower by 167-185 and the fit error
halves (0.65-0.68 %).

Limits - approximate values, a start for a fit, not reference values:

- The residuals stay non-random (runs test p ~ 1e-6): the data resolve
  structure the circuit lacks.
- A second (R|Q) splits the middle arc into ~0.28 ms and ~1.5 ms (fit
  error 0.17-0.18 %, BIC lower) on ID11 and ID45, but is not identifiable
  on ID48 (one arc runs to its bound). R2 + R3 is well determined
  (5.23 / 4.12 mOhm), the split is not (correlation -1.00).
- `Ws`, `Wo` and an anomalous Warburg were not compared; at 50 mHz the
  tail does not show yet whether it bends toward a finite length.
- ID11 differs from the other two in the arc (R2 ~30 % higher), in both
  circuits; whether that is the cell or the spread of the batch, three
  cells do not tell.
