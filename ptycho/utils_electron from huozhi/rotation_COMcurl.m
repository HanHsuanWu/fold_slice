clear;
%% load 4Ddata
% load(datafile);

data_dir = 'D:\Wuhanhsuan\20240806BTO Trial1 bulk\Spectrum Image EELS Image.npy'; %change this
data_dir = strrep(data_dir, '\', '/'); % needed this to avoid having \U

if endsWith(data_dir, '.npy')
    % np file
    dp = readNPY(data_dir);
    dp = permute(dp, [3 4 1 2]); % depend on your data
elseif endsWith(data_dir, '.mat')
    %mat file
    cbed_strut=load(strcat(data_dir,'Spectrum Image EELS Image.mat'));
    fields = fieldnames(cbed_strut);
    % Initialize an empty array to hold the concatenated result
    dp = [];
    % Concatenate all fields along the 4th dimension
    for i = 1:numel(fields)
        dp = cat(4, dp, cbed_strut.(fields{i}));
    end
    clear cbed_strut;
end

dp = double(dp);

% calculate DPC and its curl with different transformation and rotation angle
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

% plot
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

%% functions
function [COMx, COMy] = COM(diffpat)
    % calculate COM
    [nky, nkx, ny, nx] = size(diffpat);
    kx = linspace(-floor(nkx/2), ceil(nkx/2)-1, nkx);
    ky = linspace(-floor(nky/2), ceil(nky/2)-1, nky);
    [kX, kY] = meshgrid(kx, ky);
    
    COMx = zeros(ny, nx);
    COMy = zeros(ny, nx);
    
    for i = 1:1:ny
        for j = 1:1:nx
            COMx(i, j) = mean(kX .* diffpat(:, :, i, j), 'all');
            COMy(i, j) = mean(kY .* diffpat(:, :, i, j), 'all');
        end
    end
end

function [curl_abs] = MeanCurl(COMx, COMy)
    % calculate curl abs. averaged over the whole image
    [~, comx_y] = gradient(COMx);
    [comy_x, ~] = gradient(COMy);
    curl = comy_x - comx_y;
    curl_abs = mean(abs(curl), 'all');

end

function [curl_abs] = rot_search(COMx, COMy, rot)
    % calculate curl abs. for different rotation angles
    nr = length(rot);
    curl_abs = zeros(1, nr);
    for i = 1:1:nr
        r = deg2rad(rot(i));
        comx_rot = COMx * cos(r) - COMy * sin(r);
        comy_rot = COMx * sin(r) + COMy * cos(r);
    
        curl_abs(i) = MeanCurl(comx_rot, comy_rot);
    end

end