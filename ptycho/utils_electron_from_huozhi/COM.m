function [COMx, COMy] = COM(diffpat)
    % calculate COM
    % created by Haozhi Sha, 2023.10.10
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