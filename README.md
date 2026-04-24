# Stochastic Geometry-Based Cellular Network Simulation (MATLAB)

## Title  
Stochastic Geometry-Based Analysis of Coverage Probability in SISO and MISO Cellular Networks

## Introduction  

This repository contains MATLAB code used for simulating and analysing the performance of a single-tier downlink cellular network using stochastic geometry. The project was created in MATLAB R2024b. 

The results from the Monte Carlo simulations are used to validate the tractable expressions derived in the report using stochastic geometry.

## Installation Instructions  

### Required Toolboxes
- Statistics and Machine Learning Toolbox 

### How to run the software
1. Clone the repository
2. Open MATLAB
3. Set the project folder as the current working directory
4. Ensure all .m files are accessible in the MATLAB path

No external libraries are required.

## Software Overview
The table below summarises the scripts used to generates the figures in the report. Each script corresponds to system model and varied parameter.

| Parameter | SISO Monte Carlo | SISO Expression | MISO Monte Carlo | MISO Expression |
|----------|------------------|------------------|------------------|------------------|
| SINR Threshold | MonteCarloFinal | TractableExpression4 | MonteCarloBF3 | TractableExpressionBF5 |
| Base Station Density | MonteCarloBSDensity | TractableExpression4BSDensity | MonteCarloBF3BSDensity | TractableExpressionBF5BSDensity |
| SNR | MonteCarloSNR1 | TractableExpression4SNR | MonteCarloBF3SNR | TractableExpressionBF5SNR |

Beamforming beam pattern: Beamforming_simulation3

## Technical Details
The following tractable expression derived for the SISO system was used to obtain data for SISO expression curves: 	P_c(T,\lambda,\alpha,\sigma) 	=\int_{0}^{\infty} e^{-\mu \sigma^2 r^\alpha}\,e^{-\pi\lambda r^2\beta(T,\alpha)}\,2\pi\lambda r\,dr where 	\beta(T,\alpha) = \frac{2(\mu T)^{2/\alpha}}{\alpha} \mathbb{E}_G\left[g^{2/\alpha} \left(\Gamma\left(-\frac{2}{\alpha},\mu T g\right)-\Gamma\left(-\frac{2}{\alpha}\right)\right)\right].

For the MISO system curves, the following expression was used: 	P_c(T_{BF},\lambda,\alpha,\sigma) 	= \pi\lambda\int_{0}^{\infty}\exp\big(-\mu T_{BF}\sigma^2 v^{\frac{\alpha}{2}}\big)\,\exp\big(-\pi\lambda v\,\beta(T_{BF},\alpha)\big)\,dv where \beta(T_{BF},\alpha)= \frac{2(\mu T_{BF})^{\frac{2}{\alpha}}}{\alpha}\,\mathbb{E}_G\!\left[g^{\frac{2}{\alpha}}\big(\Gamma(-\tfrac{2}{\alpha},\mu T_{BF} g)-\Gamma(-\tfrac{2}{\alpha})\big)\right].

## Known Limitations
- High computational requirement for Monte Carlo simulation: Long simulation time for Monte Carlo simulations with high base station density, large simulation area, or high iteration counts
- Monte Carlo model discrepancy: Inaccuracies in low noise, low path loss cases
