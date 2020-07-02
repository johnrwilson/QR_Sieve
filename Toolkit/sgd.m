% This is a rough outline of the Stochastic Gradient Descent function
function [opt_x] = sgd(f, start, y, X, n_batches, n_epochs, learning_rate, verbose)

% f must be a function of x, y, and X and must return [level, gradient]
sample_size = round(size(y,1) / n_batches);

x = start;

if verbose
    disp(sprintf("|Epoch |Loss Function|"));
end

for i = 1:n_epochs
    
    shuffle_indices = randperm(size(y,1));
    y_shuffle = y(shuffle_indices);
    X_shuffle = X(shuffle_indices,:);
    
    for j = 1:n_batches
        begin_sample = 1 + (j-1) * sample_size;
        if j == n_batches
            y_batch = y_shuffle(begin_sample:end);
            X_batch = X_shuffle(begin_sample:end,:);
        else
            y_batch = y_shuffle(begin_sample:begin_sample + sample_size);
            X_batch = X_shuffle(begin_sample:begin_sample + sample_size,:);
        end
        
        [~, gradient] = f(x, y_batch, X_batch);
        x = x - learning_rate * gradient
        
    end
    
    if verbose && mod(i, 10) == 0
        loss = f(x, y, X);
        disp(sprintf("| % 5.2d | % 5.2d |", [i, loss]));
    end

end
opt_x = x;
