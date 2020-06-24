function y = repeat_HL(x, times)

[x_r, x_c] = size(x);

y = x(:,1);
for j = [2:x_c]
    
    for j_t = [1:times]
        
        y = [y, x(:,j)];
    
    end
    
end
