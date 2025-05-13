%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions and the time path of the variables.

%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function plot_results(par, sol, sim, firm_type)

            %% Plot Capital Policy Function
            figure;
            surf(par.Agrid, par.kgrid, squeeze(sol.k(:,:,round(par.Blen/2))));
            xlabel('Productivity (A)'); ylabel('Current Capital (k_t)'); zlabel('Next-period Capital (k_{t+1})');
            title(['Capital Policy Function - ', firm_type, ' firms']);

            %% Plot Debt Policy Function
            figure;
            surf(par.Agrid, par.kgrid, squeeze(sol.B(:,:,round(par.Blen/2))));
            xlabel('Productivity (A)'); ylabel('Current Capital (k_t)'); zlabel('Next-period Debt (B_{t+1})');
            title(['Debt Policy Function - ', firm_type, ' firms']);

            %% Plot Simulated Capital
            figure;
            plot(sim.ksim);
            xlabel('Time Period'); ylabel('Capital');
            title(['Simulated Capital Dynamics - ', firm_type, ' firms']);

            %% Plot Simulated Debt
            figure;
            plot(sim.Bsim);
            xlabel('Time Period'); ylabel('Debt');
            title(['Simulated Debt Dynamics - ', firm_type, ' firms']);

            %% Plot Simulated Investment
            figure;
            plot(sim.isim);
            xlabel('Time Period'); ylabel('Investment');
            title(['Simulated Investment Dynamics - ', firm_type, ' firms']);

            %% Plot Productivity Shocks
            figure;
            plot(sim.Asim);
            xlabel('Time Period'); ylabel('Productivity Shock');
            title(['Simulated Productivity Shocks - ', firm_type, ' firms']);
        end

        %% Heatmap of Average Capital
        function heatmap_capital(gamma_values, delta_values, avg_capital)
            figure;
            imagesc(gamma_values, delta_values, avg_capital);
            colorbar;
            xlabel('\gamma (Adjustment Cost)'); ylabel('\delta (Depreciation Rate)');
            title('Heatmap of Average Simulated Capital');
        end

    end
end
