% This is the corresponding function for fitting the non-linear elastic deformation parameters (Pc and kp)
%This code uses the basic, tension-based inactivation, or time-based inactivation models.

function f = fit_NLE(Y)

        % Part 1 - run the model:
        
        %Index out the parameters we are fitting here
        pc = Y(1);
        kp = Y(2);
        
        
        tmax = 167.4; 
        dt = 0.554;nt = tmax/dt;

        r0 = 1.05d-05; lp = 8.89e-14; 
        rt = 8.314*1000*298;
        k = 87.5; 
        c0 = 1.81; 
        epsmin = 0.1644; km = 0.23;
        kflux = 0.00236;
        
        sigmac = 5000000000000;
        
       %FOR TIME-BASED; kclose was fitted previously
%        kclose = 0.0322
        
        %FOR TENSION-BASED; Change in membrane tension
%         dsdt = 0;
%         s_prev = 0;
        
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
            
            kleak = kflux*(r0/r);
            if sigma < sigmac
                kleak = 0;
             else
                kleak = kleak;
            end
            
            %FOR TIME-BASED:
%             if sigma > sigmac 
%                 kleak0=kleak0*(1.-dt*kclose);
%             end

            %FOR TENSION-BASED
%             dsdt = (sigma-s_prev)/dt;
%             if dsdt < 0
%                 kleak=0;
%             end
            
            dc0dt=-kleak * c0;
            c0 = c0 + dc0dt * dt;
            c = c0*(r0/r)^3;
        
            drdt = lp*(rt*c-2.*(k/r)*(r-r0)/r0);
        
            if p > pc 
                drdt=drdt+kp*(p-pc)*r;
            end
            
            %Save the current stress value
            s_prev = sigma;
            
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
        yref = msl8_Data(:,2); 
        yest = outvol(:,2);
        
        % Calculate the goodness of fit, or error norm, between the measured and estimated outputs. 
        % Used the Systems Identification Toolbox Add-On to access the "goodnessOfFit()" function
        % Specify the cost function: 
            % normalized root mean squared error (NRMSE)
            % mean squared errer (MSE)
            % normalized mean squared error (NMSE)
      % This fit is the output - should be a single value for each iteration
      %This is MSE - can also get RNMSE or NMSE here if you want
        f = goodnessOfFit(yest,yref,'MSE'); 
        
        
       
           
end