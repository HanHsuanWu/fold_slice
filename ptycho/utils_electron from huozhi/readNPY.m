

function dp = readNPY(filename)
% Function to read NPY files into matlab.
% *** Only reads a subset of all possible NPY files, specifically N-D arrays of certain data types.
% See https://github.com/kwikteam/npy-matlab/blob/master/tests/npy.ipynb for
% more.
% this function is part of ePIE codes developed by Miao's group

[shape, dataType, fortranOrder, littleEndian, totalHeaderLength, ~] = readNPYheader(filename);

if littleEndian
    fid = fopen(filename, 'r', 'l');
else
    fid = fopen(filename, 'r', 'b');
end

try

    [~] = fread(fid, totalHeaderLength, 'uint8');

    % read the data
    dp = fread(fid, prod(shape), [dataType '=>' dataType]);

    if length(shape)>1 && ~fortranOrder
        dp = reshape(dp, shape(end:-1:1));
        dp = permute(dp, [length(shape):-1:1]);
    elseif length(shape)>1
        dp = reshape(dp, shape);
    end

    fclose(fid);

catch me
    fclose(fid);
    rethrow(me);
end
