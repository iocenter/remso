function [ out ] = ifEmptyZero(x,dim)
% if the element is empty, return a zero matrix of the given dimension

if isempty(x)
    out = zeros(dim);
else
    out = x;
end

end

