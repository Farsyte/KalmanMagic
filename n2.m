function s = n2(a)

    s = real(dot( a(:) , conj(a(:)) ));
