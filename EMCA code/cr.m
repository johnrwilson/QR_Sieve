function [vector] = cr(vector)
[eigvw, eigvl] = size(vector);
if eigvw>1 && eigvl>1
    error('input is a matrix')
end

if eigvl == 1
    vector = vector.';
end


