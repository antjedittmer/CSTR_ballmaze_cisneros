# Code for Paper "Nonlinear Model Predictive Control for models in quasi-LPV form" #

This folder contains the files to reproduce the results and figures (in MATLAB version R2023a) presented in

P. S. Gonzalez Cisneros & H. Werner, "Nonlinear Model Predictive Control for models in quasi-LPV form", International Journal of Robust and Nonlinear ControlVolume 30, Issue 10, Jul 2020, Pages 3801-4163

The paper can be found here:

https://onlinelibrary.wiley.com/doi/epdf/10.1002/rnc.4973

The following code can be used to recreate the figures from the paper:

* Figure 1: cstr/terminal_ingredients/compareellipsoids.m 
* Figure 2: cstr/Plotellipsoids.m 
* Figure 4: unicycle_ballmaze/Ballmaze_main.m



# Pablo's comments as a Readme #
Unicycle:

*	The original simulations used the qpOASES solver, which offered both speed and robustness. For the shared version of this code, quadprog is used instead for broader compatibility. However, quadprog is significantly slower and less reliable. To handle solver failures, the code is designed to fall back to the previous solution. If demonstrating algorithm speed is important, consider replacing quadprog with a more efficient alternative such as the Lemke algorithm, qpOASES, OSQP, or MATLAB’s mpcActiveSet (note: the latter is unavailable in older MATLAB versions and was therefore not used here).
*	At the end of the script, there is a variable animate_plot (set to false by default). When set to true, the simulation includes a pause requiring a key press to proceed. This was added to help inspect potential infeasibilities in the trajectory. The resulting animation illustrates how the prediction avoids obstacles and may be useful for presentations.

CSTR:
* In the terminal_ingredients folder, there are three files that compute terminal ingredients using the methods discussed in the paper: a constant terminal set, a BMI-based approach using the iterative method from earlier work, and the LMI-based method proposed in this paper. Precomputed .mat files are included for all three methods, as the BMI-based method (F_YZ) may take considerable time to compute.
* In the main file, the variable use_constant_PW is set to false by default. This controls whether the simulation uses constant (true) or parameter-dependent (false) matrices P and W. The parameter-dependent approach is recommended, as it yields larger ellipsoids. These were computed using the method introduced in this paper.
* The file plot_ellipsoids.m generates Figure 2 from the paper, which compares the performance of constant versus parameter-dependent terminal ingredients.


Synthesis tools:
This folder's contents have to be added to the path to calculate the terminal ingredients and to solve the cstr problem, as it is an SOCP and can be solved using quadprog. It contains yalmip and sdpt3. It is included because this is the version that works with this code.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or  FITNESS FOR A PARTICULAR PURPOSE.