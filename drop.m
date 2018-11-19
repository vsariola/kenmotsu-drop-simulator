classdef drop < handle
    
    properties (SetAccess = immutable)
        B
        R
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
        function obj = drop(B,R,s1,s2,tol)
            if nargin < 5
                tol = 1e-9;
            end
            if nargin == 1
                if length(B) >= 5
                    tol = B(5);
                end
                if length(B) >= 4
                    s2 = B(4);
                    s1 = B(3);
                    R = B(2);
                    B = B(1);
                else                
                    error('At least four params expected');
                end
            end            
            obj.B = B;
            obj.R = R;
            obj.s1 = s1;
            obj.s2 = s2;
            obj.tol = tol;
        end
        
        function params = params(obj)
            params = [obj.B obj.R obj.s1 obj.s2];
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
            a1 = obj.angle_(obj.s1);
        end
        
        function a2 = angle2(obj)
            a2 = 180-obj.angle_(obj.s2);  
        end
        
        function V = volume(obj)
            obj.calcvolume_();
            V = obj.volume_cache;        
        end
        
        function F = force(obj)
            F = -pi * (obj.B^2-1) * obj.R / 2; 
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
            obj_flipped = drop(obj.B,obj.R,-obj.s2,-obj.s1);
        end				
		
        function d = solve(d0,f,target,varargin)
            if ischar(f)
                f = @(x)x.(f);
            end            
            fun = @(x)(f(x) - target).^2/max(1e-6,abs(target));            
			d = d0.optimize_(fun,varargin{:});
		end
        
		function d = minimize(d0,f,varargin)
            if ischar(f)
                f = @(x)x.(f);
            end                  
			d = d0.optimize_(f,varargin{:});
		end
		
		function d = maximize(d0,f,varargin)
            if ischar(f)
                f = @(x)-x.(f);
            end
			d = d0.optimize_(f,varargin{:});
		end
    end
    
	methods(Static)	
        function d = segment_solve(volume,type1,value1,type2,value2,f,value)
        	d0 = drop.segment(volume,type1,value1,type2,value2);            
            d = d0.solve(f,value,'volume',[type1 '1'],[type2 '2']);
        end
    
        function d = segment_maximize(volume,type1,value1,type2,value2,f)
        	d0 = drop.segment(volume,type1,value1,type2,value2);            
            d = d0.maximize(f,'volume',[type1 '1'],[type2 '2']);
        end
        
        function d = segment_minimize(volume,type1,value1,type2,value2,f)
        	d0 = drop.segment(volume,type1,value1,type2,value2);            
            d = d0.minimize(f,'volume',[type1 '1'],[type2 '2']);
        end
            
        function obj = segment(volume,type1,value1,type2,value2)
			validatestring(type1,{'angle','radius'});
			validatestring(type2,{'angle','radius'});
			if strcmp(type1,'angle')
				if strcmp(type2,'angle')
					obj = drop.segment_aa_(volume,value1,value2);
				else
					obj = drop.segment_ar_(volume,value1,value2);
				end
			else
				if strcmp(type2,'angle')
					obj = drop.segment_ar_(volume,value2,value1).flip;
				else
					obj = drop.segment_rr_(volume,value1,value2);
				end
			end
        end
        
        function obj = cap(volume,type1,value1)
            obj = drop.segment(volume,type1,value1,'radius',0);
        end
    end
    
    methods (Access = private)
        function obj_optimized = optimize_(obj,fun,varargin)        
            constraints = @(x)obj.constraint_(x,varargin{:});          
            
            foptions = optimset('TolFun',obj.tol,'TolX',obj.tol,'Display','off','MaxIter',1000);   
            lb = [-Inf,0,-Inf,-Inf];                
			
			d_cache = [];
			d_cache_arg = [];
            [params_opt,~,exitflag] = fmincon(@myfun,obj.params,[],[],[],[],lb,[],@mycon,foptions);
            
            obj_optimized = drop(params_opt);
           
            if (exitflag < 1)                    
                warning('drop:didnotconverge','Did not converge, exitflag = %d',exitflag);
                disp(varargin);                
                for i = 1:length(varargin)
                    fprintf('%s is %g, should be %g.\n',varargin{i},obj_optimized.(varargin{i}),obj.(varargin{i}));
                end
            end 
			
			% Small optimization: share the computed drop with the constraint function
			% if it was already computed in myfun
			function check_cache(arg)
				if isempty(d_cache_arg) || any(d_cache_arg ~= arg)
					d_cache_arg = arg;
					d_cache = drop(arg);
				end	
			end
			
			function value = myfun(arg)
				check_cache(arg);
                value = fun(d_cache) ;                
			end
            
            function [c,ceq] = mycon(arg)     
				check_cache(arg);                
                c = [];
                ceq = constraints(d_cache);
            end                    
            
        end
        
        function r = r_(obj,s)
            r = sqrt(1 + obj.B^2 + 2*obj.B*cos(2*s)) * obj.R / 2;    
        end
         
        function dr = dr_(obj,s)
			dr = -(obj.B * sin(2*s))./sqrt(4*obj.B*cos(s)^2+(obj.B-1)^2) * obj.R;
        end   
        
        function dz = dz_(obj,s)
            dz = (1+obj.B*2*cos(s).^2-obj.B) ./ sqrt((obj.B-1)^2 + 4*obj.B*cos(s).^2) * obj.R;
        end         
        
        function a = angle_(obj,s)
            a = atan2d(obj.B*cos(2*s)+1,obj.B * sin(2*s));
        end
                   
        function drdz = drdz_(obj,s)
           drdz = -obj.B*sin(2*s) / (obj.B*cos(2*s)+1);
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
                area = @(s) pi*sqrt(4*obj.B*cos(s).^2+(obj.B-1)^2) .*...
                    (2*obj.B*cos(s).^2-obj.B+1)/4*obj.R^3;
                obj.volume_cache = integral(area,obj.s1,obj.s2,'AbsTol',obj.tol);
            end
        end     

        function ceq = constraint_(drop,drop2,varargin)
			n = length(varargin);
			ceq = zeros(1,n);
			for i = 1:n
                fieldname = varargin{i};				
				target = drop.(fieldname);
                value = drop2.(fieldname);				
				if strcmp(fieldname,'angle1') || strcmp(fieldname,'angle2')
					[target,value] = deal(0.1 * target - value);
                elseif strcmp(fieldname,'volume')
                    [target,value] = deal(nthroot(target,3),nthroot(value,3));
                end
                ceq(i) = target - value;				
			end
        end		
    end
    
    methods (Access = private,Static)
		function obj = segment_rr_(volume,radius1,radius2)                  
            % Volume of a segment is pi*h/6*(3*r1.^2+3*r2^2+h^2)
            % Solve h using https://en.wikipedia.org/wiki/Cubic_function
            % noting that b = 0
            a = pi/6;
            c = pi*(radius1^2+radius2^2)/2;
            d = -volume;
            delta0 = -3*a*c;
            delta1 = 27*a^2*d;
            C = nthroot((delta1+sign(delta1)*sqrt(delta1^2-4*delta0^3))/2,3);
            h = -1/(3*a)*(C+delta0/C);
            l1 = (radius2^2-radius1^2+h^2)/(2*h);
            l2 = h-l1;            
            R = sqrt(radius1^2+l1^2);                        
			obj = drop([1 R -atan2(l1,radius1) atan2(l2,radius2)]);
        end
        
        function obj = segment_ar_(volume,angle1,radius2)
            % Volume of a segment segment
            V = @(h,r1) pi*h/6.*(3*r1.^2+3*radius2.^2+h.^2);    

            % x is distance from the sphere center to the plane 2
            l1 = @(l2) -sqrt(l2^2+radius2^2)*cosd(angle1);                        
            l2 = fzero(@(l2) V(l2+l1(l2),sqrt(l2^2+radius2^2)*sind(angle1))-volume,1);
            R = sqrt(radius2^2+l2^2);  
			obj = drop([1 R (90-angle1)/180*pi atan2(l2,radius2)]);
        end
        
        function obj = segment_aa_(volume,angle1,angle2)
            if angle1+angle2 <= 180
                error('drop:SumOfContactAnglesLessThan180','The sum of contact angles should be larger than 180 (was: %g)',angle1+angle2);
            end
            R = nthroot(3*volume/pi/(4-(2-cosd(angle1))*(1+cosd(angle1))^2-...
                (2-cosd(angle2))*(1+cosd(angle2))^2),3);
            obj = drop([1 R (90-angle1)/180*pi (angle2-90)/180*pi]); 
        end
    end
end