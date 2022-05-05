
% Use this code to fit the non-linear elastic wall parameters (i.e. pc & kp)
% Note that there is a function code you need as well (fit_NLE.m).
% See the note at the bottom about how to modify that code properly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calling fminsearch
% Y0 is the initial search values of (pc,kp)
Y0 = [2e06;8e-09];
[Y,fmin] = fminsearch(@fit_NLE,Y0,[]);
disp('The minimum value is')
disp(fmin)
disp('and is located at (pc,kp) = ')
disp(Y)
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a surface plot - this is an option if you want to view the data
% The color coding is indicative of the error value
% Modification to the linear space may be necessary

% pc = linspace(0,1e10)
% kp = linspace(0,1)
% 
% for j = 1:length(pc)
%     for i = 1:length(kp)
%         f(i,j) = fit_NLE([pc(j);kp(i)]);
%     end
% end
% 
% surf(pc,kp,f)
% %set(gca,'Fontsize',24)
% colorbar;
% title('Plastic Deformation Parameter Fitting Results');
% xlabel('pc');
% ylabel('kp');
% zlabel('Error (MSE)');
% set(gca,'Cameraposition',[-100.2588 -10.804 220.895])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE FOR FUTURE USERS:
% The fitting function (fit_NLE.m) is written in such a way that you have to change
%two things in order to switch between fitting the WT and msl8 curves.
% The two things are:
%1) sigmac value - either 5 for WT or 500000+ for msl8
%2) the reference data file name - WT_Data vs. msl8_Data (towards the end of the script)


%To switch back and forth between basic, time inactivation and tension inactivation models, 
% (un)comment out the indicated code

%Last note: use the function fit_NLE_cellwall.mlx to fit the cell wall strengthening model.