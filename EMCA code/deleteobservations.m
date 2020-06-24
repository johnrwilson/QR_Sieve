function [ output ] = deleteobservations( input, deleteindex )
deleteindex = sort(deleteindex);
L_delete = length(deleteindex);
output = input;
for j_delete = [L_delete: (-1) : 1]
    i_delete = deleteindex(j_delete);
    output = [output([1:i_delete-1],:);output([i_delete+1 :  end],:)];
end



end


