T = [170 1 0.045501318772901;
     140 1 0.969456779969718;
     120 1 2.148105551373183;
     170 5 0.032850151322825;
     140 5 0.701119900098340;
     120 5 1.575656991488880;
     120 25 1.988470279487033]; 
    
abs_tol = 1e-6;
rel_tol = 1e-4;
     
warn = 0;
for i = 1:size(T,1)
    h = drop.cap_r(T(i,2),1).height;
    fprintf('Testing r=%g V=%g h=%g. Expecting: %g. ',T(i,1),T(i,2),h,T(i,3));
    d = drop.create_ar(h,T(i,2),T(i,1),1);
    F = d.force;
    fprintf('Got: %g. ',F);
    abs_err = abs(F-T(i,3));
    rel_err = abs((F-T(i,3))/T(i,3));    
    if abs_err > 1e-4 || rel_err > 1e-3
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