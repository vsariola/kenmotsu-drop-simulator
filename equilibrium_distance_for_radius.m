function [ret,angle,angle2,sphereradius,l1,l2] = equilibrium_distance_for_radius(radius1,volume,radius2)
    if (nargin < 3)
        radius2 = 1;
    end        

    % Volume of a spherical segment
    V = @(h) pi*h/6.*(3*radius1.^2+3*radius2.^2+h.^2);
    
    if radius1 > radius2        
        % Parameterize in terms of l1
        l2 = @(x) sqrt(radius1.^2+x.^2-radius2.^2);                                
        l1 = fzero(@(x) V(x+l2(x))-volume,1);
        l2 = l2(l1);        
    else
        % Parameterize in terms of l2
        l1 = @(x) sqrt(radius2.^2+x.^2-radius1.^2);        
        l2 = fzero(@(x) V(x+l1(x))-volume,1);        
        l1 = l1(l2);                      
    end            
    ret = l1+l2;    
    sphereradius = sqrt(radius1^2+l1^2);
    angle = mod(180-atand(radius1/l1),180);
    angle2 = mod(180-atand(radius2/l2),180);   
end