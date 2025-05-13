%% File Info.

%{

    solve.m
    -------
    This code solves the model.

%}

%% Solve class.

classdef solve
    methods(Static)
        %% Solve the model using VFI. 
        
        function sol = firm_problem(par)
            %% Initialization
            v0 = zeros(par.klen, par.Alen, par.Blen); % Initial guess

            v1 = nan(par.klen, par.Alen, par.Blen);
            k1 = nan(par.klen, par.Alen, par.Blen);
            B1 = nan(par.klen, par.Alen, par.Blen);
            
            crit = 1e-6;
            maxiter = 1000;
            diff = 1;
            iter = 0;

            fprintf('--- Starting Value Function Iteration with Financial Constraints ---\n\n')
            
            %% Value Function Iteration
            while diff > crit && iter < maxiter
                for ki = 1:par.klen
                    for Ai = 1:par.Alen
                        for Bi = 1:par.Blen

                            k_current = par.kgrid(ki);
                            A_current = par.Agrid(Ai);
                            B_current = par.Bgrid(Bi);

                            v_temp = -inf;
                            k_opt = 0;
                            B_opt = 0;

                            for kj = 1:par.klen
                                for Bj = 1:par.Blen

                                    k_next = par.kgrid(kj);
                                    B_next = par.Bgrid(Bj);

                                    % Check borrowing constraints
                                    if abs(B_next) <= par.Bmax

                                        % Profit calculation
                                        [profit, borrowing_next] = model.firm_profit(A_current, k_current, B_current, k_next, B_next, par);

                                        % Expected continuation value
                                        EV = 0;
                                        for Aj = 1:par.Alen
                                            EV = EV + par.pmat(Ai, Aj) * v0(kj, Aj, Bj);
                                        end

                                        value = profit + par.beta * EV;

                                        % Update if higher value is found
                                        if value > v_temp
                                            v_temp = value;
                                            k_opt = k_next;
                                            B_opt = borrowing_next;
                                        end
                                    end
                                end
                            end

                            v1(ki, Ai, Bi) = v_temp;
                            k1(ki, Ai, Bi) = k_opt;
                            B1(ki, Ai, Bi) = B_opt;

                        end
                    end
                end
                
                diff = norm(v1(:)-v0(:));
                v0 = v1;
                iter = iter + 1;

                if mod(iter,10)==0
                    fprintf('Iteration: %d, Diff: %g\n', iter, diff)
                end
            end

            fprintf('\nConverged after %d iterations.\n', iter)

            %% Store results
            sol = struct();
            sol.v = v1;
            sol.k = k1;
            sol.B = B1;
        end
        
    end
end
