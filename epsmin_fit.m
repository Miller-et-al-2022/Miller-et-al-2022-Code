% This is the main code for fitting epsmin
% Epsmin is the deformation threshold after which membrane tension begins to build
% The model tested here uses the WT initial model assumptions (no wall deformation)

x0 = 0.16;
[x,fmin] = fminsearch(@epsmin_fit_funct,x0);
disp("The minimum is")
disp(fmin)
disp("And at that point, epsmin =")
disp(x)
disp('')