function test_caps

start_test

k = @(V)nthroot(3*V/pi+sqrt(1+(3*V/pi).^2),3);
V = @(k) k - 1./k;

fun = @(V) drop.cap(V,'radius',1).height;
params = logspace(-5,5,11)';
cases = [params V(k(params))];
test('Cap height',fun,cases,2e-5);

end_test