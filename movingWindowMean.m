function [ result ] = movingWindowMean( data, windowSize )
% Calculate the mean of moving window datas 
% data: the raw data input, windowSize: the length of the moving window
[n, width]= size(data);
result = zeros(n, width);
for i = 1:n
    if i >= windowSize
        result(i,:) = sum(data(i-windowSize+1:i,:)) / windowSize;
    else
        result(i,:) = sum(data(1:i,:)) / i;
    end
end
end