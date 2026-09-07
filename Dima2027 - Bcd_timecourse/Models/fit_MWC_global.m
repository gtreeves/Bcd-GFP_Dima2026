function results = fit_MWC_global(a, Y_exp, n_vals, L_vals, kn_vals, ko_vals)
% FIT_MWC_GLOBAL fitting MWC model to data by systematic exploration
% 
%
% This function fits an MWC model to experimental dose/response using a
% global grid search over n, L, kn, and ko. For each combination of these
% four parameters, the capacity (ltot) is determined analytically.
%
% INPUTS
%   a       - Free concentration values
%   Y_exp   - Experimental correlated concentration values
%   n_vals  - Values of n to test
%   L_vals  - Values of L to test
%   kn_vals - Values of kn to test
%   ko_vals - Values of ko to test
%
% OUTPUT
%   results - Structure containing:
%       best_params  : [n, L, kn, ko, ltot]
%       best_error   : minimum sum of squared errors
%       best_Y_model : best-fit model evaluated at experimental a values
%       all_params   : all parameter combinations tested
%       all_errors   : SSE for each parameter combination
%
% MWC MODEL
%
%   Y = ltot * (a/ko) *
%       [ (1+a/ko)^(n-1)
%         + L*(ko/kn)*(1+a/kn)^(n-1) ]
%       ---------------------------------------------------------------
%       [ (1+a/ko)^n + L*(1+a/kn)^n ]
%
% For fixed n, L, kn, and ko, ltot is optimized analytically.

    %% Input checks

    a = a(:);
    Y_exp = Y_exp(:);

    if numel(a) ~= numel(Y_exp)
        error('a and Y_exp must have the same number of elements.');
    end

    %% Preallocate

    n_total = length(n_vals) * length(L_vals) * ...
              length(kn_vals) * length(ko_vals);

    % Columns: [n, L, kn, ko, ltot]
    all_params = zeros(n_total, 5);
    all_errors = zeros(n_total, 1);

    %% Initialize best fit

    best_error = Inf;
    best_params = [];
    best_Y_model = [];

    cnt = 0;

    %% Global grid search

    for koi = 1:length(ko_vals)

        ko = ko_vals(koi);

        for ni = 1:length(n_vals)

            n = n_vals(ni);

            for kni = 1:length(kn_vals)

                kn = kn_vals(kni);

                for li = 1:length(L_vals)

                    L = L_vals(li);

                    cnt = cnt + 1;

                    % Calculate MWC model shape with ltot = 1
                    f = MWC_shape(a, n, L, kn, ko);

                    % Analytical least-squares solution for ltot
                    denom = sum(f.^2);

                    if denom == 0
                        ltot = 0;
                    else
                        ltot = sum(f .* Y_exp) / denom;
                    end

                    % Model prediction
                    Y_model = ltot .* f;

                    % Sum of squared errors
                    err = sum((Y_model - Y_exp).^2);

                    % Store parameters and error
                    all_params(cnt,:) = [n, L, kn, ko, ltot];
                    all_errors(cnt) = err;

                    % Update global best fit
                    if err < best_error

                        best_error = err;
                        best_params = [n, L, kn, ko, ltot];
                        best_Y_model = Y_model;

                    end

                end
            end
        end
    end

    %% Store results

    results.best_params = best_params;
    results.best_error = best_error;
    results.best_Y_model = best_Y_model;

    results.all_params = all_params;
    results.all_errors = all_errors;

end


%% Local MWC model function

function f = MWC_shape(a, n, L, kn, ko)
%MWC_SHAPE model 

    f = (a ./ ko) .* ...
        ((1 + a ./ ko).^(n - 1) + ...
        L .* (ko ./ kn) .* (1 + a ./ kn).^(n - 1)) ./ ...
        ((1 + a ./ ko).^n + ...
        L .* (1 + a ./ kn).^n);

end