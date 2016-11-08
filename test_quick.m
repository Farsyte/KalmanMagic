%% Check error magnitudes for service functions

disp(' ');
disp('checking precision of worker functions');
disp(' ');

test_disp('hh_tri',                'ssqe', test_hh_tri(7,7,4));
test_disp('hh_tri_quad',           'ssqe', test_hh_tri_quad(7,7,4));
test_disp('hh_tri_hex',            'ssqe', test_hh_tri_hex(7,6,1));
test_disp('udut_fact',             'ssqe', test_udut_fact(7));
test_disp('udut_caat',             'ssqe', test_udut_caat(7));
test_disp('mwgso',                 'ssqe', test_mwgso(5,7));
test_disp('mwgso_mf',              'ssqe', test_mwgso_mf(5,7));

disp(' ');

test_disp('hh_tri_mex',            'ssqe', test_hh_tri_mex(7,7,4));
test_disp('hh_tri_quad_mex',       'ssqe', test_hh_tri_quad_mex(7,7,4));
test_disp('hh_tri_hex',            'ssqe', test_hh_tri_hex_mex(7,6,1));
test_disp('udut_fact_mex',         'ssqe', test_udut_fact_mex(7));

