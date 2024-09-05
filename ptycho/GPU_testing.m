% Get the number of GPUs
numGPUs = gpuDeviceCount;

% Initialize variables to track the best GPU
maxMemory = 0;
bestGPU = 1;

% Loop through each GPU to find the one with the most available memory
for i = 1:numGPUs
    gpuInfo = gpuDevice(i);
    % Update the best GPU based on available memory
    if gpuInfo.AvailableMemory > maxMemory
        maxMemory = gpuInfo.AvailableMemory;
        bestGPU = i;
    end
end
