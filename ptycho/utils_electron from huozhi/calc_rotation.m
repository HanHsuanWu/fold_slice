function [rot_input, trans_needed] = calc_rotation(dp)
%% load 4Ddata
% datafile = "F:\Collaboration\PTO_avi\TEAMI\2023.10.7\scan47\nv256_ny300_nx300.npy";
% load(datafile);
% dp = readNPY(datafile);
% dp = double(dp);

%%% The shape of the 4D dataset dp should be ndp, ndp, nscan, nscan (the first two dims are diffraction space)
%%% otherwise dp should be permuted or reshaped.

% dp = permute(dp, [3 4 1 2]); % depend on your data
% dp = reshape(dp, [48, 48, 230, 230]); % nky, nkx, ny, nx
% dp = m(:,:,:,:);

%% calculate DPC and its curl with different transformation and rotation angle
trans_method = ["no trans.", "flip dp ud", "flip dp lr", "ky,kx -> kx, ky"];
[nky, nkx, ny, nx] = size(dp);
nr = 2400;
rot_angle = linspace(-120, 120, nr);
curl_abs = zeros(2, nr);

[comx, comy] = COM(dp);
curl_abs(1, :) = rot_search(comx, comy, rot_angle);

dp_flip = flip(dp, 1);  % flip up-down
[comx, comy] = COM(dp_flip);
curl_abs(2, :) = rot_search(comx, comy, rot_angle);
clearvars dp_flip

% dp_flip = flip(dp, 2);  % flip left-right
% [comx, comy] = COM(dp_flip);
% curl_abs(3, :) = rot_search(comx, comy, rot_angle);
% clearvars dp_flip
% 
% dp_flip = permute(dp, [2 1 3 4]);  % flip left-right
% [comx, comy] = COM(dp_flip);
% curl_abs(4, :) = rot_search(comx, comy, rot_angle);
% clearvars dp_flip

[~, linearIndex] = min(curl_abs(:));
[row, col] = ind2sub(size(curl_abs), linearIndex);
rot_input = -rot_angle(col);
trans_needed = trans_method(row);

%% plot
figure; hold on
color = ["k", "b", "g", "r"];
for i = 1:size(curl_abs, 1)
    ll = plot(-rot_angle, curl_abs(i, :), 'Color', color(i));
end
legend(trans_method)
xlim([-120, 120])
xlabel('Rotation (degree)')
ylabel('Curl abs.')
text = sprintf('Rotation angle: %.1f degree\nTransformation needed: %s', rot_input, trans_needed);
title(text);

end