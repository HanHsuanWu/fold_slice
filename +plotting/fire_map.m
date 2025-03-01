%FIRE_MAP Hot blue-purlple-yellow colormap
%
% cmap = plotting.fire_map
%
% *optional*
% ** N         Size of output colormap, default is 64
%
% returns
% ++ cmap      Nx3 array with the colormap
%
%
% see also: plotting.imagesc3D, plotting.create_colormap
%
% EXAMPLES:
%   colormap(plotting.fire_map)

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function cmap = fire_map(N)
if ~exist('N','var')
    N = 64;
end

cmap = plotting.create_colormap([0 0.0625 0.25 0.38 0.5 0.8 1],[0 0 0; 0 0 0.3765; 0.5 0 1; 0.78 0  0.35; 1 0.25 0; 1 1 0; 1 1 1],N);

end
