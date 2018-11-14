classdef TestCaps < matlab.unittest.TestCase

	properties (TestParameter)
		radius = struct('small',1e-1,'medium',1,'large',10);
		volume = struct('small',1e-3,'medium',1,'large',1000);
		angle = struct('small',1,'medium',90,'large',179);
    end
    
    methods(TestMethodSetup)
        function setup(~)
            addpath('..');
        end
    end
 
    methods(TestMethodTeardown)
        function teardown(~)
            rmpath('..');
        end
    end
	
	methods (Test)
		function testCapRadius(testCase,radius,volume)
			tolerance = 1e-4;	
			V = volume / radius^3;
			k = nthroot(3*V/pi+sqrt(1+(3*V/pi).^2),3);
			height = (k - 1/k) * radius;
			d = drop.cap(volume,'radius',radius);
			testCase.verifyLessThan(abs(d.radius1 - radius) / radius,tolerance);
			testCase.verifyLessThan(abs(d.volume - volume) / volume,tolerance);
			testCase.verifyLessThan(abs(d.height - height) / height,tolerance);
			testCase.verifyLessThan(abs(d.radius2),tolerance);
			testCase.verifyLessThan(abs(d.angle2 - 180) / 180,tolerance);			
		end
		
		function testCapAngle(testCase,angle,volume)
			tolerance = 1e-4;	
			% Based on http://mathworld.wolfram.com/SphericalCap.html
			alpha = 90 - angle; 
			R = nthroot(3*volume/pi/(2-3*sind(alpha)+sind(alpha)^3),3);
			height = R*(1-sind(alpha));
			d = drop.cap(volume,'angle',angle);
			testCase.verifyLessThan(abs(d.angle1 - angle) / angle,tolerance);
			testCase.verifyLessThan(abs(d.volume - volume) / volume,tolerance);
			testCase.verifyLessThan(abs(d.height - height) / height,tolerance);
			testCase.verifyLessThan(abs(d.radius2),tolerance);
			testCase.verifyLessThan(abs(d.angle2 - 180) / 180,tolerance);	
		end
	end
end