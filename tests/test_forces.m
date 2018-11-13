function test_forces

start_test

fun = @(x) drop.create(x(1),x(2),'radius',x(3),'radius',1).force;
cases = ...
    [0.5737 1 1e-3 0.000016604854900;
     0.5737 1 1e-2 0.000095451982678;
     0.5737 1 1e-1 0.016904638685679;
     1.6596 5 1e-3 0.000038458082818;
     1.6596 5 1e-2 0.000112224642805];
test('Snap-in force for a known radius',fun,cases);

fun = @(x) drop.create(x(1),x(2),'angle',x(3),'radius',1).force;
cases = ...
    [0.5737 1 170 0.045501318772901;
     0.5737 1 140 0.969456779969718;
     0.5737 1 120 2.148105551373183;
     1.6596 5 170 0.032850151322825;
     1.6596 5 140 0.701119900098340;
     1.6596 5 120 1.575656991488880;
     3.3528 25 120 1.988470279487033]; 
test('Snap-in force for a known angle',fun,cases);

end_test
