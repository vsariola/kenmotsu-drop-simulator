T = [1e-3 1 0.5737 0.000016604854900;
     1e-2 1 0.5737 0.000095451982678;
     1e-1 1 0.5737 0.016904638685679;
     1e-3 5 1.6596 0.000038458082818;
     1e-2 5 1.6596 0.000112224642805];
     
abs_tol = 1e-6;
rel_tol = 1e-4;
     
warn = 0;
for i = 1:size(T,1)
    fprintf('Testing r=%g V=%g h=%g. Expecting: %g. ',T(i,1),T(i,2),T(i,3),T(i,4));
    F = drop.create_rr(T(i,3),T(i,2),T(i,1),1).force;
    fprintf('Got: %g. ',F);
    abs_err = abs(F-T(i,4));
    rel_err = abs((F-T(i,4))/T(i,4));    
    if abs_err > 1e-5 || rel_err > 1e-3
        fprintf('FAIL. Absolute error: %g, relative error: %g.\n',abs_err,rel_err);
        warn = 1;
    else
        fprintf('OK\n');
    end
end

if ~warn
    disp('All tests completed succesfully.');
else
    disp('Some tests failed.');    
end