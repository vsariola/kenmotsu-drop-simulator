function end_test

    if isempty(lastwarn)
        disp('All tests completed succesfully.');
    else
        disp('Some tests failed.');    
    end

    rmpath ..