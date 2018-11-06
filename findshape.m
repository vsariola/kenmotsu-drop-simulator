function [z_ret,r_ret,r_a,r_b,V,anglea,angleb] = findshape(A,C,a,b)
    if length(A) == 4 && nargin == 1
        b = A(4);
        a = A(3);
        C = A(2);
        A = A(1);
    end
    
    mu = 2/(A+C);
    m = (C^2-A^2)/2;
    n = (C^2+A^2)/2;
    
    r = @(s) sqrt(m*cos(mu*s)+n);    
    dz = @(s) 1/(A+C)*(r(s)+A*C./(r(s)));
    [z_ret,s] = cumquad(dz,a,b,1e-12);
    r_ret = r(s);        
    
    if nargout >= 3
        r_a = r(a);
    end
    if nargout >= 4
        r_b = r(b);
    end
    if nargout >= 5
        area = @(s) pi*r(s).^2.*dz(s);
        V = integral(area,a,b);
    end
    drdz = -sin(2*a/(A+C))*(A-C)/((A-C)*cos(2*a/(A+C))-A-C);
    anglea = 180-sign(drdz)*asind(1/sqrt(1+drdz^2));
    drdz = -sin(2*b/(A+C))*(A-C)/((A-C)*cos(2*b/(A+C))-A-C);
    angleb = 180+sign(drdz)*asind(1/sqrt(1+drdz^2));    
end
        
        
    
    
