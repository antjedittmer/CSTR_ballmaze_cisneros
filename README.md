# Code for Paper "Nonlinear Model Predictive Control for models in quasi-LPV form" #

This folder contains the files to reproduce the results and figures presented in

P. S. Gonzalez Cisneros & H. Werner, "Nonlinear Model Predictive Control for models in quasi-LPV form", International Journal of Robust and Nonlinear ControlVolume 30, Issue 10, Jul 2020, Pages3801-4163

The paper can be found here:

https://onlinelibrary.wiley.com/doi/epdf/10.1002/rnc.4973



# Pablo's comment as a Readme #
Unicycle:

*	for simulation, qpOASES solver  was used originally, as it was fast and robust. For sharing, this code uses quadprog, which is very slow and brittle (I added the condition to give the old solution whenever it fails). If you want to highlight the speed of the algorithms, consider changing to something else if you have it (Lemke algorithm, qpoases, osqp, or maybe matlab's built-in mpcativeset, I didn't use it because I'm using a very old version of MATLAB, which doesn't have it).
*	At the end there is a variable animate_plot set to false by default, if you set it to true, just be aware that there is a pause so you have to press a key to keep running (did this to inspect possible infeasibilities in the path), the animation actually looks pretty cool how the prediction goes around the obstacles so you can maybe use it effectively in a presentation.


CSTR:
*	In the folder terminal_ingredients, there are 3 files, each calculating the terminal ingredients with one of the methods mentioned in the paper (constant, BMI with the old papers' iterative method, and the LMI method suggested in this paper). I also included mat files with the results of each of the methods already saved, especially the F_YZ iterative method can take a while.
* In the main file, there is a variable use_constant_PW set to false by default, with which you can select to use constant P and W (true) or parameter dependent (false). The parameter-dependent is better as the ellipsoids are larger. These were calculated with the method suggested in the paper.
* The file plot_ellipsoids generates Figure 2 from the paper comparing the prediction using constant and PD terminal ingredients.

Synthesis tools:
This folder's contents have to be added to the path to calculate the terminal ingredients and to solve the cstr problem, as it is an SOCP and can be solved using quadprog. It contains yalmip and sdpt3. I included it because this is the old version, known to work with my code; maybe newer versions won't work.

