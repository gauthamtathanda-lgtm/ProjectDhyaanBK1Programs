# Fractal Dimension Analysis of EEG During Meditation

This project investigates the fractal properties of EEG signals in meditators and non-meditators using nonlinear complexity measures.

We compute two primary measures:

- Higuchi Fractal Dimension (HFD): quantifies local temporal irregularity (signal roughness)
- Hurst Exponent (H) using Detrended Fluctuation Analysis (DFA): quantifies long-range temporal correlations

Our objective is to compare neural complexity:

- Between meditators and non-meditators
- Across different stimulus conditions
- Across time within a meditation session especially within specific frequency bands (alpha and theta)

## Research Questions

1. Do meditators exhibit different broadband EEG fractal dimension compared to controls?
2. Does EEG complexity change during a meditation session?
3. Are alpha and theta bands associated with altered long-range temporal correlations?
4. Do visual grating stimuli modulate fractal complexity?


## Planned Methodology

1. Preprocessing
   - Apply bandpass filtering to reduce slow drifts and high-frequency noise
   - Examine data for clearly noisy segments or channels
   - Re-reference signals to give voltage measurements are consistently defined
     
2.Windowing

   - Divide EEG into sliding time windows

3. Feature Computation (per channel and per band)
   - Compute Higuchi Fractal Dimension
   - Compute Hurst exponent using DFA
   - Repeat analysis for:
        - Broadband EEG
        - Alpha Band (8-12 Hz)
        - Theta Band (4-8 Hz)

4. Comparative Analysis
   - Meditators vs controls
   - Within-session dynamics
   - Stimulus-based differences

5. Visualization
   - Channel-wise complexity
   - Time evolution plots
   - Band-specific analysis


## Theoretical Background (Brief)

EEG signals are nonlinear and may exhibit scale-invariant properties.

Higuchi Fractal Dimension estimates the fractal dimension of a time series directly in the time domain. Higher values indicate greater short-scale complexity.

The Hurst Exponent measures long-range temporal correlations:
- H ≈ 0.5 indicates random dynamics
- H > 0.5 indicates persistent correlations
- H < 0.5 indicates anti-persistent behavior

For one-dimensional signals, fractal dimension relates to H as:
D = 2 − H

Alpha and theta oscillations have been reported to change in opposite directions during meditation. Rather than examining only band power, we use the Hurst exponent to study whether these bands themselves differ in their long-range temporal correlation structure. Specifically, we aim to test whether theta activity becomes more persistent (higher H) while alpha becomes less persistent, indicating a shift in large-scale temporal organization within bands of different frequencies during meditation.
