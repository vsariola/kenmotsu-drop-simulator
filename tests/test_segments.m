function test_segments

start_test

fun = @(x) drop.segment(5,'angle',x(1),'radius',1).angle1;
fun2 = @(x) drop.segment(5,'angle',x(1),'radius',1).radius2;
params = linspace(1,179,10)';
test('Segment ar angle',fun,[params params]);
test('Segment ar radius',fun2,[params ones(size(params))]);

fun = @(x) drop.segment(5,'radius',1,'angle',x(1)).radius1;
fun2 = @(x) drop.segment(5,'radius',1,'angle',x(1)).angle2;
params = linspace(1,179,10)';
test('Segment ar angle',fun,[params ones(size(params))]);
test('Segment ar radius',fun2,[params params]);

fun = @(x) drop.segment(5,'angle',x(1),'angle',x(2)).angle1;
fun2 = @(x) drop.segment(5,'angle',x(1),'angle',x(2)).angle2;
params = linspace(2,179,10)';
test('Segment aa angle',fun,[params 181-params params]);
test('Segment aa angle',fun2,[params 181-params 181-params]);

params = rand(20,3);
fun_V = @(x) drop.segment(x(1),'radius',x(2),'radius',x(3)).volume;
fun_r1 = @(x) drop.segment(x(1),'radius',x(2),'radius',x(3)).radius1;
fun_r2 = @(x) drop.segment(x(1),'radius',x(2),'radius',x(3)).radius2;
test('Segment rr volume',fun_V,[params params(:,1)]);
test('Segment rr radius1',fun_r1,[params params(:,2)]);
test('Segment rr radius2',fun_r2,[params params(:,3)]);

end_test