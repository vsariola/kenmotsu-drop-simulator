% force.m file
function [Ft,zt,ut] = force(r,h,V,x0)   
    if (nargin < 4)
        x0 = [1,1,-1.5,0.5];        
    end
    % set options for the fsolve algorithm
    foptions = optimset('TolFun',1e-6,'TolX',1e-8,'Display','off','MaxIter',1000);   
    % fsolve finds the solution to the system of equations m = 0
    % argm(1) = du/dz at z = 0, argm(2) = pressure
    lb = [-Inf,0,-Inf,-Inf];    
    [optm,Fval,exitflag] = fmincon(@(x)0,x0,[0 0 1 -1],0,[],[],lb,[],@mycon,foptions)
    [c,ceq] = mycon(optm)     
    
    % fsolve didn't converge, warn and return 0 instead of returning
    % potentially erroneous values
    if (exitflag ~= 1)                    
        warning('Did not converge, exitflag = %d',exitflag);                
    end    
    
    [zt,ut] = findshape(optm);
    
    % Force
    Ft = -pi * (optm(1)^2-1) / (2*optm(2));          

    function [c,ceq] = mycon(arg)                
        [z,rr,r1,~,r2,~,Vcurrent,f_a,f_b] = findshape(arg);                
        c = [-f_a,-f_b,1e-6-min(rr)];        
        ceq = [V - Vcurrent;z(end)-h;r1-r;r2-1];
    end
end