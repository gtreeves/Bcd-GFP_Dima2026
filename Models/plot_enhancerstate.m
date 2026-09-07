function [best_params, best_SSE, y_best, y_low, y_high] = ...
    plot_enhancerstate(x, B, fraction_active, y1, y2, ...
                       params_filtered_ratio, n1, n2, state)

% PLOT_ENHANCERSTATE
%
% Uses the MWC parameters to calculate the fraction active nuclei assuming
% it is equal to the probability that the enhancer is in the active state. 
% At each AP position, plots the minimum and maximum values predicted
% across all parameter sets, along with the best-fit prediction.
%
% INPUTS
%   x                     AP coordinates (x/L)
%   B                     Bicoid concentration at each AP position
%   fraction_active       Experimental fraction of active nuclei data
%   y1                    Lower experimental error
%   y2                    Upper experimental error
%   params_filtered_ratio MWC parameter sets
%   n1                    Number of binding sites in subdomain 1
%   n2                    Number of binding sites in subdomain 2
%   state                 Enhancer state:
%                         'fullybound'
%                         'anybound'
%                         'atleastonebound'
%
% OUTPUTS
%   best_params           Best-fit parameter set
%   best_SSE              SSE of best-fit parameter set
%   y_best                Best-fit prediction
%   y_low                 Lower prediction envelope
%   y_high                Upper prediction envelope


%% Check inputs

if nargin < 9
    error(['Nine inputs are required: x, B, fraction_active, y1, y2, ' ...
           'params_filtered_ratio, n1, n2, and state.'])
end

% Convert state to lowercase for case-insensitive selection
state = lower(state);

valid_states = {'fullybound', 'anybound', 'atleastonebound'};

if ~ismember(state, valid_states)
    error(['Unknown enhancer state: %s\n' ...
           'Valid options are: fullybound, anybound, atleastonebound'], ...
           state)
end


%% Number of parameter sets

num_sets = size(params_filtered_ratio, 1);

% Store predictions and SSE for every accepted parameter set
Y_pred_all = zeros(length(x), num_sets);
SSE_all = zeros(num_sets, 1);


%% Calculate model predictions and SSE

for m = 1:num_sets

    % MWC parameters
    L  = params_filtered_ratio(m,2);
    KN = params_filtered_ratio(m,3);
    KO = params_filtered_ratio(m,4);

    % Dimensionless Bicoid concentrations
    bn = B ./ KN;
    bo = B ./ KO;

    % MWC denominator
    denom = ...
        L .* (1 + bn).^n1 .* (1 + bo).^n2 + ...
        L .* (1 + bn).^n2 .* (1 + bo).^n1 + ...
        L^2 .* (1 + bn).^(n1+n2) + ...
        (1 + bo).^(n1+n2);


    %% Enhancer-state definition

    switch state

        case 'anybound'

            f = ((1 + bo).^(n1+n2) - 1) ./ denom;

        case 'atleastonebound'

            f = ((1 + bo).^n1 - 1) .* ...
                ((1 + bo).^n2 - 1) ./ denom;

        case 'fullybound'

            f = bo.^(n1+n2) ./ denom;

    end


    %% Store prediction

    Y_pred_all(:,m) = f;

    % SSE relative to experimental fraction active data
    SSE_all(m) = sum((f - fraction_active).^2);

end


%% Find best-fit parameter set

[best_SSE, best_idx] = min(SSE_all);

best_params = params_filtered_ratio(best_idx,:);

% Best-fit prediction
y_best = Y_pred_all(:,best_idx);


%% Prediction envelope across all accepted parameter sets

y_low  = min(Y_pred_all, [], 2);
y_high = max(Y_pred_all, [], 2);



%% Plot

figure
hold on

% Range of predictions from all accepted parameter sets
x2 = [x; flipud(x)];

fill(x2, ...
     [y_low; flipud(y_high)], ...
     'g', ...
     'FaceAlpha', .25, ...
     'EdgeColor', 'none')

% Best-fit prediction
plot(x, y_best, ...
     'k--', ...
     'LineWidth', 2)

% Experimental data
errorbar(x, fraction_active, y1, y2, ...
         'ko', ...
         'MarkerFaceColor', 'k', ...
         'LineStyle', 'none', ...
         'LineWidth', 2)

xlabel('AP coordinate (x/L)')
ylabel('Fraction active')

set(gca, 'FontSize', 18)
set(gcf, 'Renderer', 'painters')

box on
hold off

end