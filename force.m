% force.m file
function [Ft,zt,ut] = force(r,h,V,x0,debug)     
    if (nargin < 4)
        [~,~,~,sphereradius,l1,l2] = equilibrium_distance_for_radius(r,V);
        A = 0;
        C = sphereradius;
        mu = 2/(A+C);
        m = (C^2-A^2)/2;
        n = (C^2+A^2)/2;
        t1=-sign(l1)*acos((r^2-n)/m)/mu;
        t2=sign(l2)*acos((1-n)/m)/mu;
        x0 = [A,C,t1,t2]; % start optimizing from equilibrium shape
    end    
    if nargin < 5        
        debug = 0;
    end
    % set options for the fsolve algorithm
    foptions = optimset('TolFun',1e-6,'TolX',1e-8,'Display','off','MaxIter',1000);   
    % fsolve finds the solution to the system of equations m = 0
    % argm(1) = du/dz at z = 0, argm(2) = pressure
    lb = [-Inf,0,-Inf,-Inf];    
    if debug
        f = @debug_plot;
    else
        f = @(x)0;
    end
    [optm,~,exitflag] = fmincon(f,x0,[0 0 1 -1],0,[],[],lb,[],@mycon,foptions);    
    
    % fsolve didn't converge, warn and return 0 instead of returning
    % potentially erroneous values
    if (exitflag ~= 1)                    
        warning('Did not converge, exitflag = %d',exitflag);                
    end    
    
    [zt,ut] = findshape(optm);
    
    % Force
    H = 1/(optm(1)+optm(2));
    B = (optm(2)-optm(1))*H;
    Ft = -pi * (B^2-1) / (2*H);          

    function ret = debug_plot(x)
        [z,u] = findshape(x);
        plot(u,z);
        drawnow;
        ret = 0;
    end
    
    function [c,ceq] = mycon(arg)                
        [z,~,r1,r2,Vcurrent] = findshape(arg);   
        %c = [arg(1)-arg(2);-arg(1)-arg(2)];
        c = [];
        ceq = [V - Vcurrent;z(end)-h;r1-r;r2-1];
    end
end