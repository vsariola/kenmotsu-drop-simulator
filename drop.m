classdef drop < handle
    
    properties (SetAccess = immutable)
        A
        C
        s1    
        s2
        tol
    end
    
    properties (Transient = true,Access = private)
        r_cache
        z_cache
        s_cache
        volume_cache
    end
    
    methods
        function obj = drop(A,C,s1,s2,tol)
            if nargin < 5
                tol = 1e-12;
            end
            if nargin == 1
                if length(A) >= 5
                    tol = A(5);
                end
                if length(A) >= 4
                    s2 = A(4);
                    s1 = A(3);
                    C = A(2);
                    A = A(1);
                else                
                    error('At least four params expected');
                end
            end            
            obj.A = A;
            obj.C = C;
            obj.s1 = s1;
            obj.s2 = s2;
            obj.tol = tol;
        end
        
        function params = params(obj)
            params = [obj.A obj.C obj.s1 obj.s2];
        end
        
        function h = height(obj)
            obj.calczs_();
            h = obj.z_cache(end);
        end
        
        function r = r(obj)
           obj.calcr_();
           r = obj.r_cache;
        end
        
        function z = z(obj)
            obj.calczs_();
            z = obj.z_cache();
        end
        
        function s = s(obj)
            obj.calczs_();
            s = obj.s_cache();
        end
        
        function r1 = radius1(obj)
            r1 = obj.r_(obj.s1);
        end
        
        function r2 = radius2(obj)
            r2 = obj.r_(obj.s2);
        end
        
        function a1 = angle1(obj)
            a1 = atan2d(obj.dr_(obj.s1),obj.dz_(obj.s1))+90;
        end
        
        function a2 = angle2(obj)
            a2 = -atan2d(obj.dr_(obj.s2),obj.dz_(obj.s2))+90;  
        end
        
        function V = volume(obj)
            obj.calcvolume_();
            V = obj.volume_cache();        
        end
        
        function F = force(obj)
            H = 1/(obj.A+obj.C);
            B = (obj.C-obj.A)*H;
            F = -pi * (B^2-1) / (2*H); 
        end
        
        function show(obj)       
            obj.calcr_();
            obj.calczs_();
            rr = [obj.r_cache fliplr(-obj.r_cache) obj.r_cache(1)];
            zz = [obj.z_cache fliplr(obj.z_cache) obj.z_cache(1)];
            plot(rr,zz,'b-'); 
            axis equal;
        end
        
        function obj_flipped = flip(obj)
            obj_flipped = drop(obj.A,obj.C,-obj.s2,-obj.s1);
        end
    end
    
	methods(Static)
        function obj = create_rr(height,volume,radius1,radius2)
        	d0 = drop.segment_rr(volume,radius1,radius2);
            constraint = @(d) [height volume radius1 radius2] -...
                [d.height d.volume d.radius1 d.radius2];
            obj = d0.optimize_(constraint);
        end
        
        function obj = create_ar(height,volume,angle1,radius2)
        	d0 = drop.segment_ar(volume,angle1,radius2);
            constraint = @(d) [height volume angle1/100 radius2] -...
                [d.height d.volume d.angle1/100 d.radius2];
            obj = d0.optimize_(constraint);
        end
        
        function obj = create_aa(height,volume,angle1,angle2)
        	d0 = drop.segment_aa(volume,angle1,angle2);
            constraint = @(d) [height volume angle1/100 angle2/100] -...
                [d.height d.volume d.angle1/100 d.angle2/100];
            obj = d0.optimize_(constraint);
        end
        
        function obj = create_ra(height,volume,radius1,angle2)
            obj = drop.create_ar(height,volume,angle2,radius1).flip;
        end
        
        function obj = segment_rr(volume,radius1,radius2)             
            if radius1 > radius2
                obj = drop.segment_rr(volume,radius2,radius1).flip;               
                return;
            end
            V = @(h) pi*h/6.*(3*radius1.^2+3*radius2.^2+h.^2);
            l1 = @(x) sqrt(radius2.^2+x.^2-radius1.^2);        
            l2 = fzero(@(x) V(x+l1(x))-volume,1);        
            l1 = l1(l2);                                             
            R = sqrt(radius1^2+l1^2);            
            s1 = -sign(l1)*acos(2*radius1^2/R^2-1)*R/2;
            s2 = sign(l2)*acos(2*radius2^2/R^2-1)*R/2;
            obj = drop([0 R s1 s2]);
        end
        
        function obj = segment_ar(volume,angle1,radius2)
            % Volume of a segment segment
            V = @(h,r1) pi*h/6.*(3*r1.^2+3*radius2.^2+h.^2);    

            % x is distance from the sphere center to the plane 2
            l1 = @(l2) -sqrt(l2^2+radius2^2)*cosd(angle1);                        
            l2 = fzero(@(l2) V(l2+l1(l2),sqrt(l2^2+radius2^2)*sind(angle1))-volume,1);       
            l1 = l1(l2);
            R = sqrt(radius2^2+l2^2);  
            radius1 = R*sind(angle1); 
            s1 = -sign(l1)*acos(2*radius1^2/R^2-1)*R/2;
            s2 = sign(l2)*acos(2*radius2^2/R^2-1)*R/2;
            obj = drop([0 R s1 s2]);
        end
        
        function obj = segment_aa(volume,angle1,angle2)
            if angle1+angle2 < 180
                error('The sum of contact angles should be larger than 180 (was: %g)',angle1+angle2);
            end
            V = @(l1,l2) pi*(l1+l2)/6.*(3*(l1*tand(angle1)).^2+3*(l2*tand(angle2)).^2+(l1+l2).^2);
            l2 = fzero(@(l2) V(l2*cosd(angle1)/cosd(angle2),l2)-volume,1);
            l1 = l2*cosd(angle1)/cosd(angle2);
            R = -l2/cosd(angle2);  
            radius1 = -l1*tand(angle1); 
            radius2 = -l2*tand(angle2);
            s1 = -sign(l1)*acos(2*radius1^2/R^2-1)*R/2;
            s2 = sign(l2)*acos(2*radius2^2/R^2-1)*R/2;
            obj = drop([0 R s1 s2]);
        end
            
        function obj = segment_ra(volume,radius1,angle2)
            obj = drop.segment_ar(volume,angle2,radius1).flip;
        end
        
        function obj = cap_r(volume,radius1)
            obj = drop.segment_rr(volume,radius1,0);
        end
        
        function obj = cap_a(volume,angle1)
            obj = drop.segment_ar(volume,angle1,0);
        end            
    end
    
    methods (Access = private)
        function obj_optimized = optimize_(obj,constraints)
            foptions = optimset('TolFun',obj.tol,'TolX',obj.tol,'Display','off','MaxIter',1000);   
            lb = [-Inf,0,-Inf,-Inf];    
            [params_opt,~,exitflag] = fmincon(@(x)0,obj.params,[],[],[],[],lb,[],@mycon,foptions);
            
            obj_optimized = drop(params_opt);
           
            if (exitflag ~= 1)                    
                warning('Did not converge, exitflag = %d',exitflag);
                disp(constraints);
                constraints(obj_optimized)
            end 
            
            function [c,ceq] = mycon(arg)                
                c = [];
                ceq = constraints(drop(arg));
            end
        end
        
        function r = r_(obj,s)
            H = 1/(obj.A+obj.C);
            B = (obj.C-obj.A)*H;
            r = sqrt(1 + B^2 + 2*B*cos(2*H*s)) / (2*H);    
        end
         
        function dr = dr_(obj,s)
            mu = 2/(obj.A+obj.C);
            m = (obj.C^2-obj.A^2)/2;
            n = (obj.C^2+obj.A^2)/2;  
            dr = 0.5 * 1./sqrt(m*cos(mu*s)+n) * m .* -sin(mu*s) * mu;
        end   
        
        function dz = dz_(obj,s)
            H = 1/(obj.A+obj.C);
            B = (obj.C-obj.A)*H;
            dz = (1+B*2*cos(H*s).^2-B) ./ sqrt((B-1)^2 + 4*B*cos(H*s).^2);
        end              
        
        function drdz = drdz_(obj,s)
           drdz = -sin(2*s/(obj.A+obj.C))*(obj.A-obj.C) / ...
               ((obj.A-obj.C)*cos(2*s/(obj.A+obj.C))-obj.A-obj.C);         
        end
        
        function calczs_(obj)  
            if isempty(obj.z_cache) || isempty(obj.s_cache)
                [z,s] = cumquad(@obj.dz_,obj.s1,obj.s2,obj.tol);            
                obj.z_cache = [0 z];
                obj.s_cache = [obj.s1 s];
            end            
        end
        
        function calcr_(obj)
            if isempty(obj.r_cache)
                obj.calczs_();
                obj.r_cache = obj.r_(obj.s_cache);            
            end
        end

        function calcvolume_(obj)
            if isempty(obj.volume_cache)
                area = @(s) pi*obj.r_(s).^2.*obj.dz_(s);
                obj.volume_cache = integral(area,obj.s1,obj.s2,'AbsTol',obj.tol);
            end
        end       
    end
end