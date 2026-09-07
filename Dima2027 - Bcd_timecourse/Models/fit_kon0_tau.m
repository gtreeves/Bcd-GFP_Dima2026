function [kon0, tau, t0, kon0_SSE, tau_SSE, t_on, P_on] = ...
    fit_kon0_tau(x_exp, fraction_active_exp, ...
                 AP_coord, mean_onsets, ...
                 L1, Kn, Ko, co, l, options)


% FIT_KON0_TAU
%
% Fits:
%   1. kon0 from steady-state fraction active data
%   2. tau and t0 from transcriptional onset time data
%
% p_full assumptions:
%
%   'anybound'
%       p_full = ((1 + bo).^9 - 1) ./ denom
%
%   'atleastonebound'
%       p_full = (((1 + bo).^5 - 1) .* ...
%                 ((1 + bo).^4 - 1)) ./ denom
%
%   'fullybound'
%       p_full = bo^9 ./ denom
%
% Inputs:
%
%   x_exp                Experimental AP positions for fraction active
%   fraction_active_exp  Experimental fraction active
%   AP_coord             AP positions for onset times
%   mean_onsets          Experimental onset times
%   L1                   Equilibrium constant between  unbound open andfany
%                        nucleosomal states
%   Kn                   Bcd affinity for nucleosomal DNA
%   Ko                   Bcd affinity for open DNA
%   co                   Free gradient amplitude
%   l                    Free gradient decay length
%   options              Fitting/model options
%
% Outputs:
%
%   kon0                 Best-fit kon0
%   tau                  Best-fit tau
%   t0                   Best-fit t0
%   kon0_SSE             SSE for kon0 fit
%   tau_SSE              SSE for tau/t0 fit
%   t_on                 Predicted onset times
%   P_on                 Probabilities of occupying the respective number of states
%




%% CONSTANTS


c5 = [1 5 10 10 5 1];

c4 = [1 4 6 4 1];

c9 = [1 9 36 84 126 126 84 36 9 1];

N = 4;



% FITTING SETTINGS


kon0_range = options.kon0_range;

tau_range = options.tau_range;

t0_range = options.t0_range;

p_full_type = options.p_full_type;

nx = options.nx;

nt = options.nt;


%% CALCULATE OCCUPANCY AT EXPERIMENTAL POSITIONS

B = co * exp(-x_exp/l);

bo = B ./ Ko;

bn = B ./ Kn;


W1 = zeros(length(B),10);

W2 = zeros(length(B),10);


for i = 0:5

    for j = 0:4

        W1(:,i+j+1) = ...
            W1(:,i+j+1) + ...
            L1 .* ...
            (c5(i+1).*bn.^i) .* ...
            (c4(j+1).*bo.^j);


        W2(:,i+j+1) = ...
            W2(:,i+j+1) + ...
            L1 .* ...
            (c4(j+1).*bn.^j) .* ...
            (c5(i+1).*bo.^i);

    end

end


W3 = (L1 * L1) .* ...
     (bn.^(0:9)) .* ...
     c9;


W4 = (bo.^(0:9)) .* ...
     c9;


b_states = [W1 W2 W3 W4];

denom = sum(b_states,2);



% SELECT p_full ASSUMPTION FOR kon0 FIT


switch p_full_type

    case 'anybound'

        % At least one Bcd molecule bound across the 9 sites

        p_full_fit = ...
            ((1 + bo).^9 - 1) ./ denom;


    case 'atleastonebound'

        % At least one Bcd molecule bound on each enhancer

        p_full_fit = ...
            (((1 + bo).^5 - 1) .* ...
             ((1 + bo).^4 - 1)) ./ denom;


    case 'fullybound'

    % All 9 Bcd binding sites occupied

    p_full_fit = W4(:,10) ./ denom;


    otherwise

        error(['Invalid p_full_type. ' ...
               'Use ''anybound'', ''atleastonebound'', ' ...
               'or ''fullybound''.']);

end



%% FIT kon0


v = [zeros(N-1,1); 1];


bestSSE = Inf;

kon0 = NaN;


for i = 1:length(kon0_range)

    kon0_test = kon0_range(i);


    P_on_model = ...
        zeros(length(p_full_fit),1);


    for k = 1:length(p_full_fit)

        kon = ...
            kon0_test * p_full_fit(k);


        A = spdiags( ...
            repmat([kon -(kon+1) 1],N,1), ...
            [-1 0 1], ...
            N, ...
            N);


        A(1,1) = -kon;

        A(end,end) = -1;

        A(end,:) = 1;


        P = A \ v;


        P_on_model(k) = P(end);

    end


    sse = sum( ...
        (P_on_model - fraction_active_exp).^2);


    if sse < bestSSE

        bestSSE = sse;

        kon0 = kon0_test;

    end

end


kon0_SSE = bestSSE;



%% CREATE SPATIAL MODEL GRID


x = linspace(0,1,nx)';

Bx = co * exp(-x/l);

bo = Bx ./ Ko;

bn = Bx ./ Kn;



%% CALCULATE OCCUPANCY ON SPATIAL GRID


W1 = zeros(nx,10);

W2 = zeros(nx,10);


for i = 0:5

    for j = 0:4

        W1(:,i+j+1) = ...
            W1(:,i+j+1) + ...
            L1 .* ...
            (c5(i+1).*bn.^i) .* ...
            (c4(j+1).*bo.^j);


        W2(:,i+j+1) = ...
            W2(:,i+j+1) + ...
            L1 .* ...
            (c4(j+1).*bn.^j) .* ...
            (c5(i+1).*bo.^i);

    end

end


W3 = (L1 * L1) .* ...
     (bn.^(0:9)) .* ...
     c9;


W4 = (bo.^(0:9)) .* ...
     c9;


b_states = [W1 W2 W3 W4];

denom = sum(b_states,2);



% SELECT SAME p_full ASSUMPTION FOR PROMOTER DYNAMICS


switch p_full_type

    case 'anybound'

        % At least one Bcd molecule bound across the 9 sites

        p_full = ...
            ((1 + bo).^9 - 1) ./ denom;


    case 'atleastonebound'

        % At least one Bcd molecule bound to each subdomain

        p_full = ...
            (((1 + bo).^5 - 1) .* ...
             ((1 + bo).^4 - 1)) ./ denom;

    case 'fullybound'

    % All 9 Bcd binding sites occupied

    p_full = ...
        W4(:,10) ./ denom;
end



%% STEADY-STATE 


P_on = zeros(nx,1);


for k = 1:nx

    kon = ...
        kon0 * p_full(k);


    A = spdiags( ...
        repmat([kon -(kon+1) 1],N,1), ...
        [-1 0 1], ...
        N, ...
        N);


    A(1,1) = -kon;

    A(end,end) = -1;

    A(end,:) = 1;


    P = A \ v;


    P_on(k) = P(end);

end



%% TIME-DEPENDENT 


t_base = linspace(0,15,nt)';


P0 = zeros(N,1);

P0(1) = 1;


P_base = zeros(nt,nx);


for k = 1:nx

    konhat = ...
        kon0 * p_full(k);


    A = spdiags( ...
        repmat([konhat -(konhat+1) 1],N,1), ...
        [-1 0 1], ...
        N, ...
        N);


    A(1,1) = -konhat;

    A(end,end) = -1;


    fhandle = @(t,P) A * P;


    [~,P_res] = ...
        ode15s( ...
            fhandle, ...
            t_base, ...
            P0);


    P_base(:,k) = P_res(:,end);

end


P_base = cummax(P_base,1);



%% FIT tau AND t0


t_cross_base = nan(nx,1);


for k = 1:nx

    theta = ...
        0.5 * P_on(k);


    p_vec = P_base(:,k);


    idx = ...
        find(p_vec >= theta,1,'first');


    if ~isempty(idx) && idx > 1

        y1v = p_vec(idx-1);

        y2v = p_vec(idx);


        frac = ...
            (theta-y1v) / ...
            (y2v-y1v);


        t_cross_base(k) = ...
            t_base(idx-1) + ...
            frac * ...
            (t_base(idx)-t_base(idx-1));

    end

end



exp_x = AP_coord ./ 100;


min_error = Inf;

tau = NaN;

t0 = NaN;


for tau_test = tau_range

    t_scaled = ...
        tau_test .* t_cross_base;


    for t0_test = t0_range

        model_t = ...
            t_scaled + t0_test;


        model_exp = ...
            interp1( ...
                x, ...
                model_t, ...
                exp_x, ...
                'linear', ...
                NaN);


        valid = ...
            ~isnan(model_exp) & ...
            ~isnan(mean_onsets);


        if any(valid)

            err = ...
                sum( ...
                    (model_exp(valid) - ...
                     mean_onsets(valid)).^2);


            if err < min_error

                min_error = err;

                tau = tau_test;

                t0 = t0_test;

            end

        end

    end

end


tau_SSE = min_error;



%% TRANSCRIPTIONAL ONSET TIMES


t_on = ...
    tau .* t_cross_base + t0;


end