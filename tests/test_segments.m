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


params = rand(20,3);
fun_V = @(x) drop.segment_rr(x(1),x(2),x(3)).volume;
fun_r1 = @(x) drop.segment_rr(x(1),x(2),x(3)).radius1;
fun_r2 = @(x) drop.segment_rr(x(1),x(2),x(3)).radius2;
test('Segment rr volume',fun_V,[params params(:,1)]);
test('Segment rr radius1',fun_r1,[params params(:,2)]);
test('Segment rr radius2',fun_r2,[params params(:,3)]);

end_test