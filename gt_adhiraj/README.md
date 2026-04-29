# Divergent Neural Dynamical Regimes in Meditators and Controls

## Overview

This project analyzes EEG data from a meditator (013AR) and a matched control (064PK) using Higuchi Fractal Dimension (HFD) as a measure of neural complexity.

The aim is to examine how neural dynamics differ under visual stimulation and meditation, beyond standard spectral analysis.

## Method

EEG signals are analyzed in broadband (1–90 Hz). For each trial:

ΔHFD = HFD_stimulus − HFD_baseline

Baseline window: −500 to 0 ms  
Stimulus window: 250 to 750 ms  

HFD is computed with kmax = 10. Statistical comparisons are performed using the Wilcoxon signed-rank test. Analysis is focused on occipito-parietal electrodes.

## Experimental Conditions

EO1 → EC1 → G1 → M1 → G2 → EO2 → EC2 → M2

EO/EC: eyes open / closed  
G1, G2: grating stimuli  
M1, M2: early and late meditation blocks  

## Results

Controls show consistent reductions in HFD under stimulus conditions.  
The meditator shows stability or increases in HFD.

For grating stimuli, controls exhibit monotonic suppression with increasing spatial frequency, while the meditator shows scale-dependent modulation.

Differences between meditator and control are minimal in M1 but become pronounced in M2, indicating that effects emerge over time rather than being present initially.

## Interpretation

The results suggest that meditation is associated with a change in neural dynamical regime, rather than a simple change in activity levels.

## Limitations

Data is epoched and short, limiting reliable estimation of long-range temporal correlations (Hurst exponent).  
Surrogate analysis (IAAFT) has not yet been implemented.  
Analysis is based on a single subject pair.

## Data

Data is derived from the ProjectDhyaanBK1 EEG dataset.

