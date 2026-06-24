# Wavelet Transform Tutorial

A hands-on MATLAB tutorial introducing the fundamentals and practical application of wavelet transforms for signal and image analysis.

This repository accompanies a series of tutorial exercises that guide students through:

- Understanding the theory behind wavelet transforms
- Performing multi-resolution signal analysis
- Applying discrete wavelet transforms (DWT)
- Visualizing time-frequency representations
- Exploring wavelet-based denoising and feature extraction
- Working with real-world data examples

## Repository Structure

```text
waveletTransform_tutorial/
│
├── code/                    # Tutorial scripts and example solutions
├── data/                    # Example datasets used throughout the exercises
├── functions/               # Helper functions for wavelet analysis
│
├── Tutorial_IIA_instructions.pdf
├── Tutorial_IIB_instructions.pdf
├── Tutorial_IIC_instructions.pdf
│
└── README.md
```

## Prerequisites

### MATLAB

The tutorials were developed using MATLAB.

Recommended MATLAB toolboxes:

- Wavelet Toolbox
- Signal Processing Toolbox


## Getting Started

1. Clone the repository:

```bash
git clone https://github.com/stefan-ehrlich/waveletTransform_tutorial.git
cd waveletTransform_tutorial
```

2. Open MATLAB.

3. Add the repository to your MATLAB path:

```matlab
addpath(genpath(pwd));
```

4. Open the tutorial instructions:

- `Tutorial_IIA_instructions.pdf`
- `Tutorial_IIB_instructions.pdf`
- `Tutorial_IIC_instructions.pdf`

5. Follow the exercises in sequence.

## Tutorial Overview

### Tutorial IIA – Fundamentals of Wavelets

Topics include:

- Fourier transform limitations
- Time-frequency analysis
- Wavelet concepts
- Scaling and translation
- Mother wavelets

### Tutorial IIB – Discrete Wavelet Transform

Topics include:

- Filter banks
- Multi-level decomposition
- Approximation and detail coefficients
- Signal reconstruction

### Tutorial IIC – Applications

Topics include:

- Signal denoising
- Compression
- Feature extraction
- Biomedical or engineering signal analysis

## Example Workflow

```matlab
% Load example data
load('data/example_signal.mat')

% Perform wavelet decomposition
[c,l] = wavedec(signal,4,'db4');

% Visualize coefficients
plot(c)

% Reconstruct approximation
approx = wrcoef('a',c,l,'db4',4);
```

## Learning Objectives

After completing the tutorials, participants should be able to:

- Explain the motivation behind wavelet transforms
- Distinguish between Fourier and wavelet analysis
- Perform wavelet decomposition and reconstruction
- Interpret wavelet coefficients
- Apply wavelet methods to practical data analysis problems

## References

Recommended reading:

1. Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.
2. Mallat, S. *A Wavelet Tour of Signal Processing*.
3. Daubechies, I. *Ten Lectures on Wavelets*.
4. MATLAB Wavelet Toolbox Documentation.

## License

Please refer to the repository license information before redistribution or modification.

## Author

Stefan Ehrlich
GitHub: https://github.com/stefan-ehrlich
