% This is a rough outline of the Stochastic Gradient Descent function
function [opt_x] = sgd(f, start, y, X, n_batches, n_epochs, learning_rate, decay, verbose)

% f must be a function of x, y, and X and must return [level, gradient]
sample_size = round(size(y,1) / n_batches);

% NOTE/TODO: If you want to specify a tolerance instead of a number of
% epochs, you should measure norm(learning_rate * gradient) / norm(x) as
% the object that must meet that criteria. This will decrease given the use
% of momentum parameter, though it will not be monotonic. For the OLS
% example, 10^-6 seemed to be a good tolerance level.

x = start;

% num_increase_loss_epochs = 0;\
start_loss = f(start, y, X);
best_loss = start_loss;
best_parameters = start;

if verbose
    disp("-----------------------------------------");
    disp(sprintf("| Epoch | Loss Function | Gradient Norm |"));
    disp("-----------------------------------------");
    s1 = sprintf("| % 5.2d | % 5.6d |", [0, start_loss]);
    s2 = sprintf("       % s     |", "N/A");
    disp(strcat(s1, s2));
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
        x = x - learning_rate * gradient;
        
    end

    learning_rate = learning_rate * decay;
    
    if mod(i, 2) == 0
        loss = f(x, y, X);
        if loss < best_loss
            best_loss = loss;
            best_parameters = x;
        else
            x = best_parameters;
            
        end
        if verbose && mod(i, 50) == 0
            disp(sprintf("| % 5.2d | % 5.6d |   % 5.2d   |", [i, best_loss, norm(gradient * learning_rate)]));
        end
    end
    
end
opt_x = x;
