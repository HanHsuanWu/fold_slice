function [curl_abs] = rot_search(COMx, COMy, rot)
    % calculate curl abs. for different rotation angles
    % created by Haozhi Sha, 2023.10.10
    nr = length(rot);
    curl_abs = zeros(1, nr);
    for i = 1:1:nr
        r = deg2rad(rot(i));
        comx_rot = COMx * cos(r) - COMy * sin(r);
        comy_rot = COMx * sin(r) + COMy * cos(r);
    
        curl_abs(i) = MeanCurl(comx_rot, comy_rot);
    end

end

