# American Option Pricing Using ADI Schemes and Predictor-Corrector FDM

This project implements numerical methods to price American-style options, focusing on American put options. It leverages Finite Difference Methods (FDM) with a predictor-corrector scheme and explores Finite Element Methods (FEM) for comparison. The primary goal is to accurately compute the early exercise boundary and option price while performing convergence analysis.

## Features

- **American Put Option Pricing**: Computes the option price and optimal exercise boundary \(S_f(t)\) using a predictor-corrector FDM.
- **Numerical Schemes**: Implements explicit Euler predictor and Crank-Nicolson corrector steps.
- **Convergence Analysis**: Quantitative assessment for both time and asset price discretizations.
- **Flexible Parameters**: Modify strike price, volatility, risk-free rate, dividend yield, and grid sizes.
- **Efficient Computation**: Uses sparse matrices for tridiagonal system solutions.

## Mathematical Model

- Black-Scholes PDE for American put options:

\[
\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + (r - D_0) S \frac{\partial V}{\partial S} - r V = 0
\]

- Terminal condition at expiry:

\[
V(S, T) = \max(K - S, 0)
\]

- Free boundary \(S_f(t)\) determines the optimal early exercise region.

## Installation

1. Install [MATLAB](https://www.mathworks.com/products/matlab.html) (R2021a or later recommended).
2. Clone this repository:

```bash
git clone https://github.com/yourusername/american-option-pricing.git
cd american-option-pricing
```

## Usage
Run the main function to compute the early exercise boundary:
```matlab
% Example: Compute American put option early exercise boundary
N = 800; % Number of time steps
M = 100; % Number of asset price steps
S_f = american_option(N, M);
disp(S_f(end)); % Display final early exercise boundary
```
## Convergence Tests
- Time Convergence:
```matlab
N_list = [200, 400, 800, 1600, 3200];
S_f_values = zeros(size(N_list));
for i = 1:length(N_list)
    N = N_list(i);
    S_f = american_option(N, 100);
    S_f_values(i) = S_f(end) * 100; 
end
disp(S_f_values);
```
- Asset Price Convergence: Similar approach with varying ```M``` and fixed ```N```.

## Parameters
| Parameter            | Value  |
| -------------------- | ------ |
| Time to maturity (T) | 1 year |
| Volatility (σ)       | 0.3    |
| Risk-free rate (r)   | 0.1    |
| Dividend yield (D₀)  | 0      |
| Strike price (X)     | 100    |
## Key References
- Zhu, S.-P., & Zhang, J. (2011). A new predictor-corrector scheme for valuing American puts. Applied Mathematics and Computation, 217(9), 4439–4452.
