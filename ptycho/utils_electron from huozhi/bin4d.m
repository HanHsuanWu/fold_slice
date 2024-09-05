% This function bins a 4D matrix x by b1 (rows) and b2 (columns), after the matrix has been loaded in the workspace
% output matrix: xbin
% input x: nv, nv, ny, nx

function xbin = bin4d(x, b1, b2)
[m1,m2,m3,m4] = size(x);
m11 = b1*floor(m1/b1);
m22 = b2*floor(m2/b2);
xbin = reshape(x(1:m11,1:m22,:,:), b1, m11/b1, b2, m22/b2, m3, m4);
xbin = sum(sum(xbin,3),1);
xbin = squeeze(xbin);
