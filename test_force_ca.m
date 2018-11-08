T = [170 1 0.5737 0.045501318772901;
     140 1 0.5737 0.969456779969718;
     120 1 0.5737 2.148105551373183;
     170 5 1.6596 0.032850151322825;
     140 5 1.6596 0.701119900098340;
     120 5 1.6596 1.575656991488880;
     120 25 3.3528 1.988470279487033]; 
    
abs_tol = 1e-6;
rel_tol = 1e-4;
     
warn = 0;
for i = 1:size(T,1)
    fprintf('Testing r=%g V=%g h=%g. Expecting: %g. ',T(i,1),T(i,2),T(i,3),T(i,4));
    F = force_ca(T(i,1),T(i,3),T(i,2));
    fprintf('Got: %g. ',F);
    abs_err = abs(F-T(i,4));
    rel_err = abs((F-T(i,4))/T(i,4));    
    if abs_err > 1e-6 || rel_err > 1e-4
        fprintf('FAIL. Absolute error: %g, relative error: %g.\n',T(i,1),T(i,2),T(i,3),abs_err,rel_err);
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