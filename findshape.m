function [z_ret,r_ret,r_a,angle_a,r_b,angle_b,V,f_a,f_b] = findshape(B,H,a,b)
    if length(B) == 4 && nargin == 1
        b = B(4);
        a = B(3);
        H = B(2);
        B = B(1);
    end
        
    r = @(s) sqrt(1+B^2+2*B*cos(2*H*s))/(2*H);    
    f = @(s) (1+B*cos(2*H*s))./(2*H*r(s));
    [z_ret,s] = cumquad(f,a,b,1e-12);
    r_ret = r(s);    
    if nargout >= 3
        r_a = r(a);
    end
    if nargout >= 4
        drdz_a = 2*B*cos(H*a)*sin(H*a)/(2*B*cos(H*a)^2-B+1);
        angle_a = acotd(-drdz_a);
    end    
    if nargout >= 5
        r_b = r(b);
    end
    if nargout >= 6
        drdz_b = 2*B*cos(H*b)*sin(H*b)/(2*B*cos(H*b)^2-B+1);
        angle_b = acotd(-drdz_b);
    end
    if nargout >= 7
        A = @(s) sqrt(4*B*cos(H*s).^2+(B-1)^2).*(2*B*cos(H*s).^2-B+1)/(4*H^2);
        V = pi * integral(A,a,b);
    end
    if nargout >= 8
        f_a = f(a);
    end
    if nargout >= 9
        f_b = f(b);
    end
end
        
        
    
    
