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
            drdz_at_s1 = obj.drdz_(obj.s1);
            a1 = 180-sign(drdz_at_s1)*asind(1/sqrt(1+drdz_at_s1^2));
        end
        
        function a2 = angle2(obj)
            drdz_at_s2 = obj.drdz_(obj.s2);
            a2 = 180+sign(drdz_at_s2)*asind(1/sqrt(1+drdz_at_s2^2));    
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
            plot(obj.r,obj.z,'b-',-obj.r,obj.z,'b-'); 
            axis equal;
        end
    end
    
	methods(Static)
        function obj = spherical_rr(volume,radius1,radius2)             
            if radius1 > radius2
                d = drop.spherical_rr(volume,radius2,radius1);
                obj = drop(d.A,d.C,-d.s2,-d.s1);
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
        
        function obj = spherical_ar(volume,angle1,radius2)
            obj = []; % TBW
        end
        
        function obj = spherical_aa(volume,angle1,angle2)
            obj = []; % TBW
        end
            
        function obj = spherical_ra(volume,radius1,angle2)
            d = drop.spherical_ar(volume,angle2,radius1);
            obj = drop(d.A,d.C,-d.s2,-d.s1);
        end
            
    end
    
    methods (Access = private)
        function r = r_(obj,s)
            mu = 2/(obj.A+obj.C);
            m = (obj.C^2-obj.A^2)/2;
            n = (obj.C^2+obj.A^2)/2;    
            r = sqrt(m*cos(mu*s)+n);    
        end
         
        function dz = dz_(obj,s)
            dz = 1/(obj.A+obj.C)*(obj.r_(s)+obj.A*obj.C ./(obj.r_(s)));
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
                obj.volume_cache = integral(area,obj.s1,obj.s2);
            end
        end
    end
end