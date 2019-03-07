function s = get01Sequence(oxInfo,length)
% Generate 01 sequence.
% Input: 
%   oxInfo: OX matrix
%   length: the corresponding vector.

if size(oxInfo,1) ~= size(length,1)
    return;
end

s = '';
stepCount = oxInfo(1,6);
for j = 1 : size(oxInfo,1)   
    if oxInfo(j,end) == 0
        a = '0';
    else
        a = '1';
    end
    if oxInfo(j,6)~=stepCount
        s = [s, '#'];
        stepCount = oxInfo(j,6);
    end
    for k = 1:ceil(length(j))
        s = [s, a];
    end    
    s = [s,'_'];
end
s = [s,'#'];
end

