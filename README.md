## Introduction

A more detailed exposition can be found in the presentation and the pdf of the paper. 
Investments are a balance of risk and rewards. Markowitz started a revolution in portfolio management that continues to this day (factor models, multi-period models, etc.). Popular KPIs within these frameworks are generally in the form of ratios:
* Sharpe ratio
* Information ratio
* Maximum diversification
* Minimum concentration

Here we:
* Develop a generic framework for KPIs in this framework
* Recast optimisation problem as a traditional ML 
* Design efficient algorithms to solve

We find that all KPIs can be reduced to:

$$\underset{w \in C}{\mathop{\text{max}}} \frac{w^T[\delta\mu\mu^T + (1-\delta)\sigma\sigma^T]w}{w^T[\gamma\Sigma + (1-\gamma)diag(\sigma^2)]w}$$

which can be rewritten with new altered weights as a Rayleigh ratio

$$\underset{w \in C}{\mathop{\text{max}}}\frac{w^T A w}{w^T B w}$$

This can be solved as a bi-convex problem and we show out of sample performace for 2 KPIs on 3 stock indices. 


