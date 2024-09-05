function [dp, pty] = io_TEAM(pty, resx, resy, probex, probey)
%% this function is part of ePIE codes developed by Miao's group

if contains(pty, '.mat')
    load(pty, 'dp');
elseif contains(pty, '.emd')
    dp = h5read(pty, '/data/data/data');
    dp = reshape(dp, [resx resy probex probey]);
elseif contains(pty, '.h5')
    data = h5read(pty, '/electron_events/frames');
    resx = h5readatt(pty, '/electron_events/frames', 'Nx');
    resy = h5readatt(pty, '/electron_events/frames', 'Ny');
    probey = h5readatt(pty, '/electron_events/scan_positions', 'Ny');
    probex = h5readatt(pty, '/electron_events/scan_positions', 'Nx')-1;
    dp = ones(resx, resy, probex, probey, 'uint8');
    rx = 1;
    ry = 1;
    for i=1:length(data)
        frame = data{i};
        dy = mod(frame, probex) + 1;
        dx = floor(frame / probex) + 1;
        dx(dx > probex) = probex;
        dp(rx, ry, dx, dy) = dp(rx, ry, dx, dy) + 1;
        rx = rx + 1;
        if rx > resx
            ry = ry + 1;
            rx = 1;
        end
    end
end

dp = single(dp);
dp(dp < 0) = 0;

pty = split(pty, '/');
pty = split(pty{end}, '.');
pty = pty{1};

end
