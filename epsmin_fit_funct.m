function f = epsmin_fit_funct(epsmin)

    % Part 1 - run the model:

    tmax = 167.4;
    dt = 0.554;nt = tmax/dt;

    r0 = 1.05d-05; lp = 8.89d-14; rt = 8.314*1000*298;
    k = 87.5; c0 = 1.81; km = 0.23;
    sigmac = 0.005; kflux = 0.00236;
    
    % These parameters are fitted in another function
%     pc = 2.56e06;
%     p = 0;
%     kp = 8.26e-9;
   
    outvol = [0 0];
    
    r=r0;
    
    %Start the time loop from 1 to nt
    for it = 1:nt 
        sigma = 0;
        if (r/r0-1) > epsmin
            sigma = 1.*km*((r/r0-1) - epsmin);
        end 
        eps = (r/r0-1);
      
        p = k*(r/r0-1)*2/r;
        
        kleak = kflux*(r0/r);
        if sigma < sigmac
            kleak = 0;
        else
            kleak = kleak;
        end
        
        dc0dt=-kleak * c0;
        c0 = c0 + dc0dt * dt;
        c = c0*(r0/r)^3;
    
        drdt = lp*(rt*c-2.*(k/r)*(r-r0)/r0);
    
%         if p > pc 
%             drdt=drdt+kp*(p-pc)*r;
%         end
        
        r = r+drdt*dt;
        t = it*dt;
        %time and volume change are written down at each time step - outvol is what is compared to the data
        vol = (r/r0)^3-1.;
        newvol = [t vol];
        outvol = [outvol;newvol];
    
    end

% Part 2 - Compare the model results to the experimental data
    % Load in the experimental data - pull out the "reference y"
    % Prepare the estimated data that we'll compare it to
    load WT_Data.csv;
    yref = WT_Data(:,2);
    yest = outvol(:,2);
    
    % Calculate the goodness of fit, or error norm, between the measured and estimated outputs. 
    % Used the Systems Identification Toolbox Add-On to access the "goodnessOfFit()" function
    % Specify the cost function: 
        % normalized root mean squared error (NRMSE)
        % mean squared error (MSE)
        % normalized mean squared error (NMSE)

    %This is MSE - can change the cost function if you want.
    f = goodnessOfFit(yest,yref,'MSE'); 
    % This fit is the output - should be a single value for each iteration
   
       
end