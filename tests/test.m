function test(name,fun,cases,tol)
    if nargin < 4
        tol = 1e-6;
    end
    fprintf('Test: %s\n',name);
    for i = 1:size(cases,1)
        fprintf('  ');
        params = cases(i,1:(end-1));
        expected = cases(i,end);
        for j = 1:length(params)
            if j > 1
                fprintf(',');
            end
            fprintf('%g',params(j));          
        end
        fprintf(' Expect: %g',expected);
        got = fun(params);
        error = abs(got - expected);
        fprintf(' Got: %g ',got);
        if error > tol
            warning('Test failed, error was %g (tol: %g).',error,tol);
        else
            fprintf('OK.\n');
        end
    end
end
    