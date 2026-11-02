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
