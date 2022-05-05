% This is the corresponding function for fitting the non-linear elastic deformation parameters (Pc and kp)
%This code uses the cell wall strengthening model.

function f = fit_NLE_cellwall(Y)

        % Part 1 - run the model:
        
        %Index out the parameters we are fitting here
        pc = Y(1);
        kp = Y(2);
        
        
        tmax = 167.4; 
        dt = 0.554;nt = tmax/dt;
        
        r0 = 1.05d-05; lp = 1.043e-13; rt = 8.314*1000*298;
        k = 87.5; c0 = 1.55; epsmin = 0.126; km = 0.23;
        
        %The pressure value needs to start at 0
        p = 0;
        
        r = r0;
        
        outvol = [0 0];
        
        %Start the time loop from 1 to nt
        for it = 1:nt 
            sigma = 0;
            if (r/r0-1) > epsmin
                sigma = 1.*km*((r/r0-1) - epsmin);
            end 
            eps = (r/r0-1);
          
            p = k*(r/r0-1)*2/r;
            
            kleak = 0;
            dc0dt=-kleak * c0;
            c0 = c0 + dc0dt * dt;
            c = c0*(r0/r)^3;
        
            drdt = lp*(rt*c-2.*(k/r)*(r-r0)/r0);
        
            if p > pc 
                drdt=drdt+kp*(p-pc)*r;
            end
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
        load msl8_Data.csv;
        yref = msl8_Data(:,2); %For selecting just the beginning of the curve, index out 1-111
        yest = outvol(:,2);
        
        % Calculate the goodness of fit, or error norm, between the measured and estimated outputs. 
        % Used the Systems Identification Toolbox Add-On to access the "goodnessOfFit()" function
        % Specify the cost function: 
            % normalized root mean squared error (NRMSE)
            % mean squared errer (MSE)
            % normalized mean squared error (NMSE)
      % This fit is the output - should be a single value for each iteration
      %This is SSE
        %f = sum(((yest(:)-yref(:))).^2);
        %f = sum(((outvol(:,2)-msl8_Data(:,2)).^2));
      %This is MSE - can also get RNMSE or NMSE here if you want
        f = goodnessOfFit(yest,yref,'MSE'); 
        
        
       
           
end