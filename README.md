# Adaptive Spline Fitting with Automatic Knot Selection

## Project Overview
This project implements an adaptive method for B-spline curve and surface fitting. The core idea is to approximate data points using as few knots as possible while maintaining high accuracy. The algorithm, based on the work of Kang et al., operates in two stages:
1.  **Sparse Fitting:** Solves a convex optimization problem with $l_1$-norm (or mixed $l_{\infty,1}$-norm) regularization to identify "active" knots.
2.  **Knot Adjustment:** Refines the position of the selected knots to minimize the least-squares error.


## Key Features
* **Sparse Optimization:** Uses `cvx` to minimize the control polygon's derivative jumps, effectively selecting significant knots automatically.
* **Extensions Implemented:**
    * **Parametric Curves:** Support for both open (clamped) and closed (periodic) curves in 2D/3D space.
    * **Noise Robustness:** Tested on datasets with Gaussian noise to evaluate stability.
    * **Discontinuities:** Detection and fitting of discontinuous signals using the approach by Storath and Weimann.
    * **Surface Fitting:** Extension to tensor-product B-spline surfaces using mixed $l_{\infty,1}$ regularization.

    
* **Alternative Algorithms:** Implementation of a **Differential Evolution (DE)** algorithm for global heuristic knot placement (in presesence of noise).

## Technologies & Methods
* **Language:** MATLAB
* **Packages:** CVX (for convex optimization), MATLAB DE package (for an alternative approach).
* **Math Concepts:** B-Splines, Convex Optimization, Least Squares, Sparse Regularization, Numerical Analysis.

## Results
The method successfully reduces the number of parameters (knots) needed to represent complex geometries.
* **Cubic vs. Quintic Splines:** Compared performance between degrees $p=3$ and $p=5$.
* **Surface Reconstruction:** Successfully reconstructed 3D surfaces starting from a dense $40 \times 40$ knot grid, reducing it to essential active knots only.

## Execution
To reproduce the numerical experiments, please refer to the **`CODE_GUIDE`** file included in this repository.


## 👥 Authors
* **Chiara Garuglieri** 
* Ilaria Girotti
* Agnese Pizzi

## 📚 References

The algorithms and methodologies implemented in this project are based on the following academic papers:

1. H. Kang, F. Chen, Y. Li, J. Deng, and Z. Yang.
    *Computer-Aided Design*, 58, 179-188, (2015).

2.  J. Luo, H. Kang, and Z. Yang.
    *Computer Aided Geometric Design*, 73, 54-69, (2019).

3.  J. Luo, H. Kang, and Z. Yang.
    *Journal of Computational Mathematics*, 40, 589-606, (2021).

4.  V. Goepp, O. Bouaziz, and G. Nuel.
    *Computational Statistics and Data Analysis*, 202, 108043, (2025). 

5.  M. Buehren.
    *MATLAB Central File Exchange*, (2025). Available at: (https://it.mathworks.com/matlabcentral/fileexchange/18593-differential-evolution).

6.  M. Storath and A. Weinmann.
    *Journal of Computational and Graphical Statistics*, 33:2, 651-664, (2024). 

7.  S. Volz, M. Storath, and A. Weinmann.
    *Information and Inference: A Journal of the IMA*, 1-32, (2025). 
