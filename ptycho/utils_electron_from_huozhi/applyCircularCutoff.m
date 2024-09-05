function [dp] = applyCircularCutoff(dp, cutoff)
%Apply circular cutoff to diffraction patterns
%   Yi Jiang (yj245@cornell.edu)

N_dp = size(dp,1);
x = (-fix(N_dp/2):N_dp-fix(N_dp/2)-1);
[Y,X] = meshgrid(x,x);
R = sqrt(X.^2+Y.^2);

for i=1:size(dp,3)
    for j=1:size(dp,4)
        temp = dp(:,:,i,j);
        temp(R>cutoff)=0;
        dp(:,:,i,j) = temp;
    end
end

end