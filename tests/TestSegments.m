classdef TestSegments < matlab.unittest.TestCase
    
    properties (TestParameter)        
        volume = struct('small',1e-2,'medium',10,'large',1e4);
        angle1 = struct('small', 1, 'medium', 90, 'large', 179, 'max', 180);        
        angle2 = struct('small', 1, 'medium',90, 'large', 179, 'max', 180);        
        radius1 = struct('min',0,'small', 0.1, 'medium', 1, 'large', 10);        
        radius2 = struct('min',0,'small', 0.1, 'medium', 1, 'large', 10);        
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
         function testAngleAngle(testCase, volume, angle1, angle2)
             tolerance = 1e-4;
             fun = @()drop.segment(volume,'angle',angle1,'angle',angle2);                
             if (angle1+angle2<=180)                
                 testCase.verifyError(fun,'drop:SumOfContactAnglesLessThan180')
             else
                 d = fun();
                 testCase.assertLessThan(abs(d.volume - volume)/volume,tolerance);
                 testCase.assertLessThan(abs(d.angle1 - angle1),tolerance);
                 testCase.assertLessThan(abs(d.angle2 - angle2),tolerance);
             end
         end
        
        function testAngleRadius(testCase, volume, angle1, radius2)
            tolerance = 1e-4;
            d = drop.segment(volume,'angle',angle1,'radius',radius2);                            
            testCase.assertLessThan(abs(d.volume - volume)/volume,tolerance);
            testCase.assertLessThan(abs(d.angle1 - angle1),tolerance);
            testCase.assertLessThan(abs(d.radius2 - radius2),tolerance);           
        end
        
        function testRadiusAngle(testCase, volume, radius1, angle2)
            tolerance = 1e-4;
            d = drop.segment(volume,'radius',radius1,'angle',angle2);                            
            testCase.assertLessThan(abs(d.volume - volume)/volume,tolerance);
            testCase.assertLessThan(abs(d.radius1 - radius1),tolerance);
            testCase.assertLessThan(abs(d.angle2 - angle2),tolerance);           
        end
        
        function testRadiusRadius(testCase, volume, radius1, radius2)
            tolerance = 1e-4;
            d = drop.segment(volume,'radius',radius1,'radius',radius2);                            
            testCase.assertLessThan(abs(d.volume - volume)/volume,tolerance);
            testCase.assertLessThan(abs(d.radius1 - radius1),tolerance);
            testCase.assertLessThan(abs(d.radius2 - radius2),tolerance);           
        end        
    end
end  