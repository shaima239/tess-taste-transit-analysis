# Transit Light Curve Analysis of WASP‑12b with TASTE and TESS Observations

This repository contains the data reduction, transit modeling, and Bayesian inference code used for the analysis of the hot Jupiter WASP‑12b. The study combines ground‑based photometry from the TASTE project (Asiago 1.82 m telescope) with space‑based observations from NASA’s TESS mission (Sectors 43, 44, and 45). The goal is to refine the planetary parameters, particularly the radius, orbital period, and transit depth, and to compare the results with published values.

## Data Sources

- **TASTE** (Asiago Search for Transit Timing Variation of Exoplanets)  
  Ground‑based data from the Asiago 1.82 m telescope, observed on 2019‑12‑29.  
  *Filter:* Sloan r.

- **TESS** (Transiting Exoplanet Survey Satellite)  
  Sectors 43, 44, and 45 (2021‑09‑16 to 2021‑11‑07), 2‑minute cadence.  
  *Products used:* Target Pixel Files (TPF) and Light Curve Files (LCF).  
  *Detrending:* Manual filtering with `wotan` (Huber spline, biweight, window sizes) to remove instrumental and systematic signals.

## Methods

- **TASTE data reduction**  
  - Bias and flat‑field correction with error propagation  
  - Barycentric time correction (JD → BJD_TDB)  
  - Centroid refinement and aperture photometry (95% enclosed flux)  
  - Differential photometry using two reference stars  

- **TESS data preprocessing**  
  - Removal of NaN/infinite values and bad quality flags (pessimistic selection)  
  - Detrending with `wotan` (masked transit regions)  
  - Selection of optimal flux type (PDCSAP for Sector 43, SAP for Sectors 44 & 45)  

- **Transit modeling**  
  - `batman` for Mandel‑Agol quadratic limb‑darkened transit model  
  - Limb‑darkening coefficients from LDTk (Monte Carlo sampling)  
  - Circular orbit assumed  

- **Bayesian inference**  
  - Initial optimization with `pyde` (Differential Evolution)  
  - MCMC sampling with `emcee` (28 walkers, 20 000 steps)  
  - Priors: Gaussian for limb‑darkening, uniform for other parameters  
  - Posterior analysis with corner plots  

## Key Results

| Parameter                        | This Work                     | Literature (Ozturk & Erdem 2019) |
|----------------------------------|-------------------------------|----------------------------------|
| \(T_c\) (BJD)                    | \(2454509.000 \pm 0.006\)     | \(2454508.97824\)                |
| \(P\) (days)                     | \(1.091418 \pm 0.000001\)     | \(1.0914199\)                    |
| \(R_p/R_\ast\)                   | \(0.117112 \pm 0.000649\)     | \(0.119\)                        |
| \(a/R_\ast\)                     | \(2.97994 \pm 0.04167\)       | \(2.933\)                        |
| \(i\) (deg)                      | \(82.145 \pm 0.688\)          | \(81.92\)                        |
| Transit depth                    | \(0.013715 \pm 0.000152\)     | –                                |
| Planet radius (\(R_J\))          | \(1.88 \pm 0.016\)            | \(1.90\) (approx.)               |
| Planetary density (g/cm³)        | \(0.263 \pm 0.069\)           | \(0.311\)                        |

The differences in \(a/R_\ast\), \(i\), and density may reflect different modeling choices, detrending methods, or physical evolution of the planet (mass loss, Roche‑lobe overflow).


