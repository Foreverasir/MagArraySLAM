function [output] = movingMeanFilt(data,windowSize)
% A moving mean filter for convenience.More args should be considered.
% data: the data to be filted. windowSize: length of window.
a = 1;
b = (1 / windowSize) * ones(1, windowSize);
output = filter(b, a, data, [], 1);
output = output(windowSize:size(data,1), :);
end

