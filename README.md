# bSPOD – Band-Ensemble Spectral Proper Orthogonal Decomposition

MATLAB implementation of **Band-Ensemble Spectral Proper Orthogonal Decomposition (bSPOD)** with frequency attribution.

Described in:  
**Band-Ensemble Spectral Proper Orthogonal Decomposition with Frequency Attribution**  
Jakob G.R. von Saldern, Oliver T. Schmidt, Philipp Godbersen, J. Moritz Reumschüssel, Tim Colonius

**Contact:** Jakob G. R. von Saldern, j.vonsaldern@tu-berlin.de

---

## Description

bSPOD estimates SPOD modes from consecutive Fourier coefficients obtained
from a single Fourier transform of the full time record, avoiding time
segmentation. Compared to Welch-based SPOD, bSPOD reduces spectral
leakage, permits increased frequency resolution, and retains frequency
information of tonal components at comparable computational cost.

These properties enable reduced estimator variance while maintaining low
bias for tonal components, making bSPOD particularly effective for
broadband–tonal flows.

---
## Function

```matlab
[GAINS, MODES, FREQS] = bSPOD(data, Nf)
[GAINS, MODES, FREQS, EXPANSION] = bSPOD(data, Nf, ...)
```

## Inputs

```text
data    : (Nt x Ndof) time-series data (time in rows, DOFs in columns)
Nf      : integer, number of consecutive Fourier bins per frequency band / number of POD modes

Options (name–value pairs):
dt      : time step (default: 1)
ell     : frequency attribution exponents (default: [inf 1 0])
weights : spatial weights vector (Ndof x 1)
ovlp    : number of overlapping Fourier bins between bands (default: 0)
window  : apply Hann window before FFT (0 or 1, default: 0)

        -> the number of frequnecy bands is 
          Nwin = floor((Nt − Nf)/HOP) + 1 with the band hop size HOP = Nf − ovlp
```

## Outputs

```text
GAINS     : (Nwin x Nf) band-wise modal gains  (number of frequnecy bands x number of modes)
MODES     : (Nwin x Nf x Ndof) spatial mode shapes
FREQS     : (Nwin x Nf x (numel(ell))) frequencies for each exponent in ell
EXPANSION : (Nwin x Nf x Nf) band-wise eigenvectors (optional)
```

## Citation

```bibtex
@article{vonSaldern2026bspod,
  title   = {Band-Ensemble Spectral Proper Orthogonal Decomposition with Frequency Attribution},
  author  = {von Saldern, Jakob G. R. and Schmidt, Oliver T. and Godbersen, Philipp and Reumsch{\"u}ssel, J. Moritz and Colonius, Tim},
  journal = {},
  year    = {2026}
}
```

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).
See the `LICENSE` file for details.
