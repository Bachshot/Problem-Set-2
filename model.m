%% File Info.

%{

    model.m
    -------
    This code sets up the model.

%}

%% Model class.

classdef model
    methods(Static)
        %% Set up structure array for model parameters and set the simulation parameters.
        
function par = setup(firm_type)
            %% Structure array for model parameters.
            par = struct();
            
            %% Economic parameters
            par.beta = 0.96;            % Discount factor
            par.alpha = 0.33;           % Capital's share (Cobb-Douglas)
            par.delta = 0.08;           % Depreciation rate
            par.gamma = 0.10;           % Adjustment cost parameter
            par.r = 0.05;               % Interest rate on debt
            
            %% Financial constraints (firm-specific)
            if strcmp(firm_type, 'small')
                par.Bmax = 10;          % Tight borrowing limit for small firms
            elseif strcmp(firm_type, 'large')
                par.Bmax = 100;         % Loose borrowing limit for large firms
            else
                error('Specify firm type as either "small" or "large".');
            end

            %% Productivity shocks
            par.sigma_eps = 0.07;       % Std. deviation of productivity shocks
            par.rho = 0.85;             % Persistence of AR(1)
            par.mu = 0;                 % Mean of productivity

            %% Simulation parameters
            par.seed = 2025;
            par.T = 1000;

            %% Generate grids
            par = model.gen_grids(par);
        end

        %% Generate state grids
        function par = gen_grids(par)
            %% Capital grid
            par.klen = 300;
            par.kmax = 30.0;
            par.kmin = 1e-4;
            par.kgrid = linspace(par.kmin, par.kmax, par.klen)';
            
            %% Productivity grid (Tauchen method)
            par.Alen = 7;
            par.m = 3;
            [Agrid, pmat] = model.tauchen(par.mu, par.rho, par.sigma_eps, par.Alen, par.m);
            par.Agrid = exp(Agrid);
            par.pmat = pmat;
            
            %% Debt grid (explicitly added determinant)
            par.Blen = 50;
            par.Bgrid = linspace(-par.Bmax, par.Bmax, par.Blen)';
        end

        %% Tauchen method
        function [y, pi] = tauchen(mu, rho, sigma, N, m)
            ar_mean = mu/(1-rho);
            ar_sd = sigma / sqrt(1-rho^2);
            y1 = ar_mean - m * ar_sd;
            yn = ar_mean + m * ar_sd;
            y = linspace(y1, yn, N);
            d = y(2) - y(1);

            ymatk = repmat(y, N, 1);
            ymatj = mu + rho * ymatk';

            pi = normcdf(ymatk, ymatj - d/2, sigma) - normcdf(ymatk, ymatj + d/2, sigma);
            pi(:,1) = normcdf(y(1), mu + rho*y - d/2, sigma);
            pi(:,N) = 1 - normcdf(y(N), mu + rho*y + d/2, sigma);
        end

        %% Revenue (production) function
        function output = production(A, k, par)
            output = A .* k .^ par.alpha;
        end

        %% Cost function (adjustment costs + investment)
        function cost = total_cost(k_next, k, par)
            invest = k_next - (1 - par.delta) .* k;
            adj_cost = (par.gamma / 2) .* ((invest ./ k) .^ 2) .* k;
            cost = adj_cost + invest;
        end

        %% Firm optimization problem (including borrowing constraint)
        function [profit, borrowing_next] = firm_profit(A, k, B, k_next, B_next, par)
            rev = model.production(A, k, par);
            cost = model.total_cost(k_next, k, par);
            financial_cost = (1 + par.r) .* B - B_next;

            profit = rev - cost - financial_cost;

            % Enforce borrowing constraint
            borrowing_next = min(max(B_next, -par.Bmax), par.Bmax);
        end

    end
end
