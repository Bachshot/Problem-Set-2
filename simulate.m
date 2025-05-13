%% File Info.

%{

    simulate.m
    ----------
    This code simulates the model.

%}

%% Simulate class.

classdef simulate
    methods(Static)
        %% Simulate the model. 
        
        function sim = firm_dynamics(par, sol)
            %% Set initial conditions
            rng(par.seed);
            
            T = par.T;
            Asim = zeros(T,1);
            ksim = zeros(T,1);
            Bsim = zeros(T,1);
            isim = zeros(T,1);

            % Initialize from stationary distribution
            Aind = randsample(par.Alen, 1, true, par.pmat(1,:));
            kind = round(par.klen/2);
            Bind = round(par.Blen/2);

            Asim(1) = par.Agrid(Aind);
            ksim(1) = par.kgrid(kind);
            Bsim(1) = par.Bgrid(Bind);

            %% Simulate model
            for t = 1:T-1

                % Find optimal next-period capital and debt
                k_next = sol.k(kind, Aind, Bind);
                B_next = sol.B(kind, Aind, Bind);

                % Record investment
                isim(t) = k_next - (1-par.delta)*ksim(t);

                % Draw next-period productivity
                cum_prob = cumsum(par.pmat(Aind,:));
                draw = rand;
                Aind_next = find(draw <= cum_prob, 1);

                % Update indices for next period
                kind = find(abs(par.kgrid - k_next) == min(abs(par.kgrid - k_next)), 1);
                Bind = find(abs(par.Bgrid - B_next) == min(abs(par.Bgrid - B_next)), 1);
                Aind = Aind_next;

                % Store states
                Asim(t+1) = par.Agrid(Aind);
                ksim(t+1) = k_next;
                Bsim(t+1) = B_next;

            end

            %% Store simulation results
            sim = struct();
            sim.Asim = Asim;
            sim.ksim = ksim;
            sim.Bsim = Bsim;
            sim.isim = isim;

        end

    end
end
