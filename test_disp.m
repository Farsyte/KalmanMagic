function test_disp(name, etag, value)

    if (value > eps*4)
        disp(sprintf('%-40s %9.4e %s FAIL', name, value, etag));
    else
        disp(sprintf('%-40s %9.4e %s pass', name, value, etag));
    end

