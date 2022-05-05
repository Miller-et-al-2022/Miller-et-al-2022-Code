% This is code for the basic model.
% These are the basic assumptions 
% (non-linear elastic cell wall, membrane unfolding, traditional osmotic safety valve)

% Simulation length and time steps
tmax = 150; dt = 0.554; nt = tmax/dt;
% Define parameters - physical properties of the cell
r0 = 1.05d-05;
k = 87.5; 
km = 0.23;
c0 = 1.81; 
cext = 0;
lp = 8.89d-14; 
epsmin = 0.1644;

% Define parameters - thermodynamic value constants
rt = 8.314*1000*298;

% Define parameters - channel properties
sigmac_WT = 0.005;
sigmac_msl8 = 5000000000000;
kflux = 0.00236;

% Define parameters - non-linear elastic cell wall
pc = 2.56e06;
p = 0;
kp = 8.26e-9;

%Making matrices to print the output data into
outvol_WT = [0 0];
outmem_WT = [0 0];
outvol_msl8 = [0 0];
outmem_msl8 = [0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run WT simulation first
r = r0;

%Start the time loop from 1 to nt
for it = 1:nt 
    %Remove any initial stress
    sigma = 0;
    
    %Then the stress becomes nonzero if the strain (r/r0-1) exceeds epsmin
    if (r/r0-1) > epsmin
        sigma = 1.*km*((r/r0-1) - epsmin);
    end 

    %Pressure calculation - this is for plastic deformation
    p = k*(r/r0-1)*2/r;
    
    %Scale the MS leakage rate by an inverse of the cell radius
    kleak = kflux*(r0/r);
    %Set the MS channel leakage rate to zero if stress is less than the opening stress
    if sigma < sigmac_WT
        kleak = 0;
    end
    dc0dt=-kleak * c0;
    
   
    %Integrate the first-order differential equation
    c0 = c0 + dc0dt * dt;
    
    %The actual concentration is taken to be proportional to 1/r^3 (i.e. volume)
    c = c0*(r0/r)^3;

    %Calculate the rate of change of radius as a difference between the osmotic pressure term 
    % and the elastic force from the cell wall (membrane is ignored - too small)
    %The external osmotic concentration (cext) can be modified here.
    drdt = lp*(rt*c-rt*cext-2.*(k/r)*(r-r0)/r0);
    
    %Non-linear elastic deformation of the cell wall
    if p > pc 
        drdt=drdt+kp*(p-pc)*r;
    end
    
    %Integrate the first order ODE for radius
    r = r+drdt*dt;
    %Give total time in terms of time step and time increment
    t = it*dt;
    %Time and volume change are written down at each time step
    vol = (r/r0)^3-1.;
    newvol = [t vol];
    outvol_WT = [outvol_WT;newvol];
    %Time and membrane stress are written down at each time step
    newmem = [t sigma];
    outmem_WT = [outmem_WT;newmem];
   
end
%Final volume change is written at the end of the time loop
disp('The final relative volume change for WT is:' ) 
disp((r/r0)^3-1.)


%%%% Now run the msl8 simulation
%Start the cell radius at its initial value
r = r0
% Restate c0, p, and the tension-related parameters since they changed during the WT simulation
c0 = 1.81;
p = 0;
dsdt = 0;
s_prev = 0;

%Start the time loop from 1 to nt
for it = 1:nt 
    %Remove any initial stress
    sigma = 0;
    %Then the stress becomes nonzero if the strain (r/r0-1) exceeds epsmin
    if (r/r0-1) > epsmin
        sigma = 1.*km*((r/r0-1) - epsmin);
    end 
    %Calculate strain (epsilon)
    eps = (r/r0-1);
    
    %Added pressure calculation - this is for plastic deformation
    %p=2*sigma/r;
    p = k*(r/r0-1)*2/r;
    
    %Scale the MS leakage rate by an inverse of the cell radius
    kleak = kflux*(r0/r);
    %Set the MS channel leakage rate to zero if stress is less than the opening stress
    if sigma < sigmac_msl8
        kleak = 0;
    end
    
    dc0dt=-kleak * c0;
    %Integrate the first-order differential equation
    c0 = c0 + dc0dt * dt;
    
    %The actual concentration is taken to be proportional to 1/r^3 (i.e. volume)
    c = c0*(r0/r)^3;

    %Calculate the rate of change of radius as a difference between the osmotic pressure term 
    % and the elastic force from the cell wall (membrane is ignored - too small)
    %The external osmotic concentration (cext) can be modified here.
    drdt = lp*(rt*c-rt*cext-2.*(k/r)*(r-r0)/r0);
    
    %Non-linear elastic deformation of the cell wall
    if p > pc 
        drdt=drdt+kp*(p-pc)*r;
    end
    
    %Integrate the first order ODE for radius
    r = r+drdt*dt;
    %Give total time in terms of time step and time increment
    t = it*dt;
    %Time and volume change are written down at each time step
    vol = (r/r0)^3-1.;
    newvol = [t vol];
    outvol_msl8 = [outvol_msl8;newvol];
    %Time and membrane stress are written down at each time step
    newmem = [t sigma];
    outmem_msl8 = [outmem_msl8;newmem];
   
end
%Final volume change is written at the end of the time loop
disp('The final relative volume change for msl8 is:' ) 
disp((r/r0)^3-1.)

%For display purposes - Plot the volume change over time of both the data and simulation
plot(outvol_WT(:,1),outvol_WT(:,2), 'linewidth', 2)
hold on
plot(outvol_msl8(:,1),outvol_msl8(:,2), 'linewidth', 2)
load WT_Data.csv;
plot(WT_Data(:,1),WT_Data(:,2), 'linewidth', 2)
xlabel('Time (sec)')
ylabel('Relative Volume Change')
ylim([0,1.0])
xlim([0,150])

hold on
load msl8_Data.csv;
plot(msl8_Data(:,1),msl8_Data(:,2), 'linewidth', 2, 'color', 'green')


legend('WT model','msl8 model', 'WT', 'msl8-5', 'Location','southeast')

%For display purposes - plot the membrane tension
figure;
plot(outmem_WT(:,1),outmem_WT(:,2), 'linewidth', 2)
hold on
plot(outmem_msl8(:,1),outmem_msl8(:,2), 'linewidth', 2)
xlabel('Time (sec)')
ylabel('Membrane Tension (N/m)')