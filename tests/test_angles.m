function test_segments

start_test

fun = @(x) drop.segment_ar(5,x(1),1).angle1;
params = linspace(1,179,10)';
cases = [params params];
test('Segment ar angle',fun,cases);

fun = @(x) drop.segment_ar(5,x(1),1).radius2;
params = linspace(1,179,10)';
cases = [params ones(size(params))];
test('Segment ar radius',fun,cases);

fun = @(x) drop.segment_aa(5,x(1),x(2)).angle1;
fun2 = @(x) drop.segment_aa(5,x(1),x(2)).angle2;
params = linspace(2,179,10)';
test('Segment aa angle',fun,[params 181-params params]);
test('Segment aa angle',fun2,[params 181-params 181-params]);

end_test