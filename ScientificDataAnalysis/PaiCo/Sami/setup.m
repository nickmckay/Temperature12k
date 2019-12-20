function setup

mex common/fastrandperm.cpp -output common/fastrandperm
mex common/inplaceadd.cpp -output common/inplaceadd
mex paico/private/paico_mex.cpp -output paico/private/paico_mex
disp('Done.');