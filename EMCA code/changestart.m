function [ start_new ] = changestart( start,lower,upper )
L_arg = length(start);
for j_arg = [1:L_arg]
    if start(j_arg) <= lower(j_arg) || start(j_arg) >= upper(j_arg)
         start_new(j_arg) = 1/2*(lower(j_arg)+upper(j_arg));
    else
        start_new(j_arg) = start(j_arg);
    end
end


end


