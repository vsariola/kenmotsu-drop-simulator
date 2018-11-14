classdef TestForces < matlab.unittest.TestCase

	properties (TestParameter)
		volume = {1 1 5 5 5 25 25 25 100 100 100};
		radius = {1e-2 1e-1 1e-2 1e-1 1 1e-2 1e-1 1 1e-2 1e-1 1};
		angle = {179 150 179 150 90 179 150 90 179 150 90};
		rsnap_force = {0.000069707852128   0.016842232394507   0.000056351128791   0.010553096479367  3.330514148760286   0.000028161535132   0.004517989635302   1.038936912840385 0.000014855782323   0.002282756713972   0.431098734844533};
		rmax_force = {0.062298868251724   0.583300078069553   0.062279674177058   0.575962817222772  3.501752416101278   0.062488688623505   0.594895821705963   3.480196833091670 0.062614238355868   0.606901134600560   4.386198535465340};
		asnap_force = {0.000251757881186   0.529625524682726   0.000199875623673   0.379399772424148 3.128633121188764   0.000286566962231   0.494867130601875   3.752316799622716  0.000418285878518   0.677432469208828   4.632517820620887};
		amax_force = {NaN   0.529625014279598  NaN   0.466038093064407 3.141687497155645   NaN   0.744683191498650   3.949955692960168 NaN   1.173278962359322   4.743947039728787};
        % NaNs should be replaced with correct values, once we know them.
    end
    
    methods(TestMethodSetup)
        function setup(~)
            addpath('..');
            warning off backtrace;
        end
    end
 
    methods(TestMethodTeardown)
        function teardown(~)
            rmpath('..');
        end
    end
	
    methods (Static)
        function [d,height] = snapin(volume,type,value)
            k = nthroot(3*volume/pi + sqrt(1 + (3*volume/pi)^2),3);
            height = k - 1 / k;            
			d = drop.segment_solve(volume,type,value,'radius',1,'height',height);
        end       
    end
    
	methods (Test, ParameterCombination='sequential')
		function testRadiusSnapin(testCase,volume,radius,rsnap_force)
			[d,height] = testCase.snapin(volume,'radius',radius);
			tolerance = 1e-3;			
			testCase.verifyLessThan(abs(d.height - height) / height,tolerance);
			testCase.verifyLessThan(abs(d.volume - volume) / volume,tolerance);
			testCase.verifyLessThan(abs(d.radius1 - radius) / radius,tolerance);
			testCase.verifyLessThan(abs(d.radius2 - 1),tolerance);
			testCase.verifyLessThan(abs(d.force - rsnap_force) / rsnap_force,tolerance);
		end
		
		function testRadiusMaxforce(testCase,volume,radius,rmax_force)
			tolerance = 1e-3;
			d = drop.segment_maximize(volume,'radius',radius,'radius',1,'force');
			testCase.verifyLessThan(abs(d.volume - volume) / volume,tolerance);
			testCase.verifyLessThan(abs(d.radius1 - radius) / radius,tolerance);
			testCase.verifyLessThan(abs(d.radius2 - 1),tolerance);
			testCase.verifyLessThan(abs(d.force - rmax_force) / rmax_force,tolerance);
		end
		
		function testAngleSnapin(testCase,volume,angle,asnap_force)
			tolerance = 1e-3;
			[d,height] = testCase.snapin(volume,'angle',angle);
			testCase.verifyLessThan(abs(d.height - height) / height,tolerance);
			testCase.verifyLessThan(abs(d.volume - volume) / volume,tolerance);
			testCase.verifyLessThan(abs(d.angle1 - angle) / angle,tolerance);
			testCase.verifyLessThan(abs(d.radius2 - 1),tolerance);
			testCase.verifyLessThan(abs(d.force - asnap_force) / asnap_force,tolerance);
		end
		
		function testAngleMaxforce(testCase,volume,angle,amax_force)
			tolerance = 1e-3;
			d = drop.segment_maximize(volume,'angle',angle,'radius',1,'force');
			testCase.verifyLessThan(abs(d.volume - volume) / volume,tolerance);
			testCase.verifyLessThan(abs(d.angle1 - angle) / angle,tolerance);
			testCase.verifyLessThan(abs(d.radius2 - 1),tolerance);
			testCase.verifyLessThan(abs(d.force - amax_force) / amax_force,tolerance);         
		end
	end
end
