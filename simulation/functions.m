
classdef functions

    methods(Static)

        function sensor_list = gen_sensor_list(type, options)

            arguments
                type;

                options.beta_max = c.beta_max;
                options.bandwidth = c.B_total;
                options.p_k = c.p_k;
            end

            K = 4;
            sensor_list = struct([]);
            v = [0.25, 0.75, 1.25, 1.75];
            for k = 1:K
                sensor_list(k).beta_max = options.beta_max;
                sensor_list(k).B_total = options.bandwidth;
                sensor_list(k).d_k = c.d_k;
                sensor_list(k).theta_k = randn(1,1) * 2 * pi;
                CN = c.sigma_k * randn(1, 1, 'like', 1j);
                NLoS = sqrt(1 / (c.kappa + 1)) * CN;
                LoS = sqrt(c.kappa / (c.kappa + 1)) * c.sigma_k * exp(1j * sensor_list(k).theta_k);
                sensor_list(k).h_k = LoS + NLoS;
                sensor_list(k).H_k = sqrt(c.A0) * sensor_list(k).d_k^(-0.5*c.alpha) * sensor_list(k).h_k;

                % Initial variables awaiting for optimization
                sensor_list(k).f_dt_k = c.C_DT / K;
                sensor_list(k).b_k = options.bandwidth / K;
                sensor_list(k).lam = [0.1, 0.1];

                % sensor constants
                sensor_list(k).D_k = 4e3;
                sensor_list(k).f_s_k = 200e3;
                sensor_list(k).p_k = options.p_k;

                switch type
                    case scenario.uniform
                        
                    case scenario.p_k
                        sensor_list(k).p_k = v(k) * options.p_k;

                    case scenario.D_k
                        sensor_list(k).D_k = v(k) * 4e3;

                    case scenario.f_s_k
                        sensor_list(k).f_s_k = v(k) * 200e3;
                    otherwise
                        
                end
                
            end
        end

        function r_k = r_k(sensor)
            r_k = sensor.b_k * log(1 + functions.c_r_k(sensor) / sensor.b_k) / log(2);
        end

        function beta_star = best_beta(sensor)
            beta_star = min(...
                (1/c.epsilon) * log((sensor.f_s_k / (c.epsilon * sensor.D_k)) * (sensor.lam(1) + sensor.lam(2))), ...
                sensor.beta_max...
                );
            beta_star = max(1.0, beta_star);
        end

        function mu_star = best_mu(sensor)
            if sensor.lam(2) <= 0
                mu_star = 1.0;
            else
                mu_star = sqrt(sensor.lam(2) * functions.r_k(sensor) / sensor.D_k);
            end
        end

        function eta = comp_complexity(beta)
            eta = exp(beta * c.epsilon) - exp(c.epsilon);
        end

        function T_comp = T_comp(sensor)
            beta_star = functions.best_beta(sensor);
            T_comp = sensor.D_k * functions.comp_complexity(beta_star) / sensor.f_s_k;
        end

        function T_tr = T_tr(sensor)
            r_k = functions.r_k(sensor);
            beta_star = functions.best_beta(sensor);
            T_tr = sensor.D_k / (r_k * beta_star);

        end

        function T_DT = T_DT(sensor)
            T_DT = c.c_k * (sensor.D_k / sensor.f_dt_k);
        end

        function T_bs = T_bs_sensor(sensor)
            T_bs = functions.T_comp(sensor) + functions.T_tr(sensor) + functions.T_DT(sensor);
        end

        function T_bs = T_bs(sensor_list)
            T_bs = functions.T_bs_sensor(sensor_list(1));
            for k = 1:length(sensor_list)
                if functions.T_bs_sensor(sensor_list(k)) > T_bs
                    T_bs = functions.T_bs_sensor(sensor_list(k));
                end
            end
        end

        function re = c_r_k(sensor)
            % constant part of r_k function: |H_k|^2 * p_k / N_0
            re = (abs(sensor.H_k) ^ 2) * sensor.p_k / c.N0;
        end

        function sensor_list = leader_optimization(sensor_list, options)

            arguments
                sensor_list 
                options.fixed_b_k = false
                options.fixed_f_DT_k = false
            end

            K = length(sensor_list);
            cvx_begin
            cvx_solver mosek;
            cvx_quiet TRUE;

            variable T_DT(K) nonnegative;
            variable T_tr(K) nonnegative;
            variable T nonnegative;
            if ~ options.fixed_b_k
                variable b(K) nonnegative;
            else
                b = repelem(sensor_list(1).B_total/K, K);
            end
            
            if ~ options.fixed_f_DT_k
                variable f(K) nonnegative;
            else
                f = repelem(c.C_DT/K, K);
            end

            minimize(T);

            subject to
            sum(b) <= sensor_list(1).B_total;
            sum(f) <= c.C_DT;
            for k = 1:K
                f(k) >= c.c_k * sensor_list(k).D_k * inv_pos(T_DT(k));
                -rel_entr(b(k), b(k) + functions.c_r_k(sensor_list(k))) * inv_pos(log(2)) >= sensor_list(k).D_k * inv_pos(functions.best_beta(sensor_list(k)) * T_tr(k));
                T >= T_DT(k) + T_tr(k) + functions.T_comp(sensor_list(k));
            end
            cvx_end
            for k = 1:K
                sensor_list(k).b_k = b(k);
                sensor_list(k).f_dt_k = f(k);
            end
        end

        function g_sgn = g_sgn(sensor)
            g = sensor.g;
            g_sgn = [0,0];
            for k = 1:2
                if g(k) >= 0
                    g_sgn(k) = 1;
                else
                    g_sgn(k) = -1;
                end
            end
        end

        function g = g(sensor)
            beta_star = functions.best_beta(sensor);
            mu_star = functions.best_mu(sensor);
            g = [1 - beta_star, ...
                1/mu_star - beta_star];
        end

        function sensor = update_lambda(sensor, h, g, original_lam)

            % update lambda
            for k = 1:2
                sensor.lam(k) = max(0, original_lam(k) + h * g(k));
            end

        end

        function dual = dual(sensor)
            beta = functions.best_beta(sensor);
            mu = functions.best_mu(sensor);
            g = functions.g(sensor);
            dual = (sensor.D_k / sensor.f_s_k) * exp(beta * c.epsilon) + (sensor.D_k / functions.r_k(sensor)) * mu + dot(sensor.lam, g);
        end

        function sensor = linesearch(sensor)
            h = 0.1;

            g = functions.g(sensor);
            original_lam = sensor.lam;

            dual = functions.dual(sensor);

            sensor = functions.update_lambda(sensor, h, g, original_lam);
            new_dual = functions.dual(sensor);

            while new_dual > dual 
                dual = functions.dual(sensor);
            
                h = 2 * h;
                sensor = functions.update_lambda(sensor, h, g, original_lam);
                new_dual = functions.dual(sensor);
            end

            sensor = functions.update_lambda(sensor, h/2, g, original_lam);
        end

        function beta_list = beta_list(sensor_list)
            K = length(sensor_list);
            beta_list = [];
            for k = 1:K
                beta_list(end+1) = functions.best_beta(sensor_list(k));
            end
        end

        function T_comp_list = T_comp_list(sensor_list)
            K = length(sensor_list);
            T_comp_list = [];
            for k = 1:K
                T_comp_list(end+1) = functions.T_comp(sensor_list(k));
            end
        end

        function T_tr_list = T_tr_list(sensor_list)
            K = length(sensor_list);
            T_tr_list = [];
            for k = 1:K
                T_tr_list(end+1) = functions.T_tr(sensor_list(k));
            end
        end

        function T_DT_list = T_DT_list(sensor_list)
            T_DT_list = [];
            K = length(sensor_list);
            for k = 1:K
                T_DT_list(end+1) = functions.T_DT(sensor_list(k));
            end
        end

        function T_bs_list = T_bs_list(sensor_list)
            T_bs_list = [];
            K = length(sensor_list);
            for k = 1:K
                T_bs_list(end+1) = functions.T_bs_sensor(sensor_list(k));
            end
        end

        function [b, f, beta, tStruct] = run_one_time(sensor_list, options)
            arguments
                sensor_list 
                options.fixed_b_k = false
                options.fixed_f_DT_k = false
            end
            K = length(sensor_list);

            t = [0];

            delta_t = 10000;
            l = 1;

            while delta_t > c.xi && l < c.L_max
                sensor_list = functions.leader_optimization(sensor_list, "fixed_b_k", options.fixed_b_k, "fixed_f_DT_k", options.fixed_f_DT_k);

                t(end + 1) = functions.T_bs(sensor_list);

                for k = 1:K
                    sensor_list(k) = functions.linesearch(sensor_list(k));
                end

                t(end + 1) = functions.T_bs(sensor_list);

                delta_t = abs(t(end) - t(end-1));
                l = l+1;
            end

            t(1) = [];
            tStruct.total = t(end);
            tStruct.comp_list = functions.T_comp_list(sensor_list);
            tStruct.tr_list = functions.T_tr_list(sensor_list);
            tStruct.DT_list = functions.T_DT_list(sensor_list);
            b = [sensor_list.b_k];
            f = [sensor_list.f_dt_k];
            beta = functions.beta_list(sensor_list);
        end

        function [b_avg, f_avg, beta_avg, t_avg] = run_many_times(n, options)

            arguments
                n

                options.type = scenario.uniform;
                options.beta_max = c.beta_max;
                options.bandwidth = c.B_total;
                options.p_k = c.p_k;
                options.fixed_b_k = false
                options.fixed_f_DT_k = false
            end

            b_list = zeros(4, n);
            f_list = zeros(4, n);
            beta_list = zeros(4, n);
            t_total = zeros(1, n);
            t_comp_list = zeros(4, n);
            t_DT_list = zeros(4, n);
            t_tr_list = zeros(4, n);
            parfor k = 1:n

                % print the current iteration
                fprintf('iteration %d; beta_max: %.2f; bandwidth: %e; p_k: %e', k, options.beta_max, options.bandwidth, options.p_k);
                if options.fixed_b_k
                    fprintf('; fixed b_k') 
                end

                if options.fixed_f_DT_k
                    fprintf('; fixed f_DT_k') 
                end

                fprintf('\n')
               
                sensor_list = functions.gen_sensor_list(options.type, "beta_max", options.beta_max, "bandwidth", options.bandwidth, "p_k", options.p_k);
                
                [b, f, beta, tStruct] = functions.run_one_time(sensor_list, "fixed_b_k", options.fixed_b_k, "fixed_f_DT_k", options.fixed_f_DT_k);
                t_total(1, k) = tStruct.total;
                t_comp_list(:, k) = tStruct.comp_list;
                t_tr_list(:, k) = tStruct.tr_list;
                t_DT_list(:, k) = tStruct.DT_list;
                b_list(:, k) = b;
                f_list(:, k) = f;
                beta_list(:, k) = beta;
            end
            b_avg = mean(b_list, 2);
            f_avg = mean(f_list, 2);
            beta_avg = mean(beta_list, 2);
            t_avg.total = mean(t_total, 2);
            t_avg.comp_list = mean(t_comp_list, 2)';
            t_avg.tr_list = mean(t_tr_list, 2)';
            t_avg.DT_list = mean(t_DT_list, 2)';
        end

        function [t_max, beta_avg] = varied_beta_max(beta_max_list, options)
            arguments
                beta_max_list;

                options.runs = 100;
                options.fixed_b_k = false
                options.fixed_f_DT_k = false
            end

            K = length(beta_max_list);

            if K == 0
                t_max = null;
                beta_avg = null;
                return
            end

            % 1. Define all the field names in a cell array
            field_names = {
                'default', ...
                'fixed_b_k', ...
                'fixed_f_DT_k', ...
                'fixed_both', ...
                'beta_max_1'
            };

            % 2. Define the initial value once
            init_array = zeros(1, K);

            % 3. Initialize the structs (good practice)
            beta_avg = struct();
            t_max = struct();

            % 4. Loop through the field names and assign the value
            for i = 1:length(field_names)
                field = field_names{i};
                beta_avg.(field) = init_array;
                t_max.(field) = init_array;
            end

            for k = 1:K
                [b_avg, f_avg, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", beta_max_list(k));
                beta_avg.default(k) = mean(beta_avg_k);
                t_max.default(k) = t_avg.total;
                [b_avg, f_avg, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", beta_max_list(k), "fixed_b_k", TRUE);
                beta_avg.fixed_b_k(k) = mean(beta_avg_k);
                t_max.fixed_b_k(k) = t_avg.total;
                [b_avg, f_avg, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", beta_max_list(k), "fixed_f_DT_k", TRUE);
                beta_avg.fixed_f_DT_k(k) = mean(beta_avg_k);
                t_max.fixed_f_DT_k(k) = t_avg.total;
                [b_avg, f_avg, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", beta_max_list(k), "fixed_b_k", TRUE, "fixed_f_DT_k", TRUE);
                beta_avg.fixed_both(k) = mean(beta_avg_k);
                t_max.fixed_both(k) = t_avg.total;
                [b_avg, f_avg, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", 1);
                beta_avg.beta_max_1(k) = mean(beta_avg_k);
                t_max.beta_max_1(k) = t_avg.total;
            end
        end

        function [t_max, beta_avg] = run_one_vary(vary_list, vary_type, options)
            arguments
                vary_list;
                vary_type;
                options.runs = 100;
            end

            K = length(vary_list);
            if K == 0
                t_max = null;
                beta_avg = null;
                return
            end

            beta_avg = zeros(1, K);
            t_max = zeros(1, K);
            
            beta_max_list = repelem(c.beta_max, K);
            p_k_list = repelem(c.p_k, K);
            bandwidth_list = repelem(c.B_total/4, K);

            switch vary_type
                case VaryType.beta_max
                    beta_max_list = vary_list;
                case VaryType.B
                    bandwidth_list = vary_list;
                case VaryType.p_k
                    p_k_list = vary_list;
            end

            K = length(vary_list);

            % 1. Define all the field names in a cell array
            field_names = {
                'default', ...
                'fixed_b_k', ...
                'fixed_f_DT_k', ...
                'fixed_both', ...
                'beta_max_1'
            };

            % 2. Define the initial value once
            init_array = zeros(1, K);

            % 3. Initialize the structs (good practice)
            beta_avg = struct();
            t_max = struct();

            % 4. Loop through the field names and assign the value
            for i = 1:length(field_names)
                field = field_names{i};
                beta_avg.(field) = init_array;
                t_max.(field) = init_array;
            end

            for k = 1:K
                [~, ~, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", beta_max_list(k), "bandwidth", bandwidth_list(k), "p_k", p_k_list(k));
                beta_avg.default(k) = mean(beta_avg_k);
                t_max.default(k) = t_avg.total;
                [~, ~, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", beta_max_list(k), "bandwidth", bandwidth_list(k), "p_k", p_k_list(k), "fixed_b_k", true);
                beta_avg.fixed_b_k(k) = mean(beta_avg_k);
                t_max.fixed_b_k(k) = t_avg.total;
                [~, ~, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", beta_max_list(k), "bandwidth", bandwidth_list(k), "p_k", p_k_list(k), "fixed_f_DT_k", true);
                beta_avg.fixed_f_DT_k(k) = mean(beta_avg_k);
                t_max.fixed_f_DT_k(k) = t_avg.total;
                [~, ~, beta_avg_k, t_avg] = functions.run_many_times(options.runs, "beta_max", beta_max_list(k), "bandwidth", bandwidth_list(k), "p_k", p_k_list(k), "fixed_b_k", true, "fixed_f_DT_k", true);
                beta_avg.fixed_both(k) = mean(beta_avg_k);
                t_max.fixed_both(k) = t_avg.total;
                if vary_type == VaryType.beta_max
                    beta_avg.beta_max_1(k) = 1.0;
                    t_max.beta_max_1(k) = t_max.default(1);
                    continue
                end
                [~, ~, ~, t_avg] = functions.run_many_times(options.runs, "beta_max", 1, "bandwidth", bandwidth_list(k), "p_k", p_k_list(k));
                beta_avg.beta_max_1(k) = 1.0;
                t_max.beta_max_1(k) = t_avg.total;
            end
            
        end

        function [t_max, beta_avg] = run_multiple_varies(n, options)
            arguments (Input)
                n 
                options.beta_max_list = c.beta_max_list;
                options.p_k_list = c.p_k_list;
                options.bandwidth_list = c.bandwidth_list;
            end

            [t_max.p_k, beta_avg.p_k] = functions.run_one_vary(options.p_k_list, VaryType.p_k, "runs", n);

            [t_max.B, beta_avg.B] = functions.run_one_vary(options.bandwidth_list, VaryType.B, "runs", n);

            [t_max.beta, beta_avg.beta] = functions.run_one_vary(options.beta_max_list, VaryType.beta_max, "runs", n);

        end


    end
end