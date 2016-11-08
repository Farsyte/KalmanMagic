function test_disp_xe(name,xchk,vchk,xref,vref,z,H,z_,x_)

    n2err = n2([xchk-xref vchk-vref]);
    if (n2err < 1e-4)
        disp(sprintf('%-40s %9.4e %s', name, n2err, 'pass'));
    else
        disp(sprintf('%-40s %9.4e %s', name, n2err, 'FAIL'));
        plot_zxe(['results from ' name],z,H,xchk,vchk,z_,x_);
    end
