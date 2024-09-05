function [curl_abs] = MeanCurl(COMx, COMy)
    % calculate curl abs. averaged over the whole image
    % created by Haozhi Sha, 2023.10.10
    [~, comx_y] = gradient(COMx);
    [comy_x, ~] = gradient(COMy);
    curl = comy_x - comx_y;
    curl_abs = mean(abs(curl), 'all');

end

