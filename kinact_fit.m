% This is the main code for fitting kinact
% kinact is the rate at which MSL8 inactivated
% This is fitted to the WT data because that is the only simulation that has MSL8


x0 = 0.03;
[x,fmin] = fminsearch(@kinact_fit_funct,x0);
disp("The minimum is")
disp(fmin)
disp("And at that point, kclose =")
disp(x)
disp('')