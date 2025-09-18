
classdef functions

    methods(Static)

        function r_k = r_k(sensor)
            r_k = sensor.b_k * log(1 + functions.c_r_k(sensor) / sensor.b_k) / log(2);
        end

        function beta_star = best_beta(sensor)
            beta_star = min(...
                (1/c.epsilon) * log((c.f_S / (c.epsilon * sensor.D_k)) * (sensor.lam1 + sensor.lam2)), ...
                1.0...
                );
        end

        function mu_star = best_mu(sensor)
            if sensor.lam2 <= 0
                mu_star = 1.0;
            else
                mu_star = sqrt(sensor.lam2 * functions.rate(sensor) / sensor.D_k);
            end
        end

        function eta = comp_complexity(beta)
            eta = exp(beta * c.epsilon) - exp(c.epsilon);
        end

        function T_comp = T_comp(sensor)
            beta_star = functions.best_beta(sensor);
            T_comp = sensor.D_k * functions.comp_complexity(beta_star) / c.f_S;
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
            re = (abs(sensor.H_k) ^ 2) * c.p_k / c.N0;
        end

        function sensor_list = leader_optimization(sensor_list)
            K = length(sensor_list);
            cvx_begin
                variable T_DT(K) nonnegative;
                variable T_tr(K) nonnegative;
                variable T nonnegative;
                variable b(K) nonnegative;
                variable f(K) nonnegative;

                minimize(T);

                subject to
                    sum(b) <= c.B_total;
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

        function sensor_list = T_DT_optimization(sensor_list)
            K = length(sensor_list);
            cvx_begin
                variable T nonnegative
                variable f(K) nonnegative

                minimize(T);

                subject to
                    sum(f) <= c.C_DT;
                    for k = 1:K
                        f(k) >= c.c_k * sensor_list(k).D_k * inv_pos(T); %#ok<*VUNUS>
                    end
            cvx_end
            for k = 1:K
                sensor_list(k).f_dt_k = f(k);
            end
        end

        function sensor_list = T_tr_optimization(sensor_list)
            K = length(sensor_list);
            cvx_begin
                variable T nonnegative;
                variable b(K) nonnegative;

                minimize(T);

                subject to
                   for k = 1:K
                        -rel_entr(b(k), b(k) + functions.c_r_k(sensor_list(k))) >= sensor_list(k).D_k * inv_pos(functions.best_beta(sensor_list(k)) * T);
                   end
                   sum(b) <= c.B_total;
                
                
            cvx_end
            for k = 1:K
                sensor_list(k).b_k = b(k);
            end
        end

        function g = g(sensor)
            beta_star = functions.best_beta(sensor);
            mu_star = functions.best_mu(sensor);
            g = [1 - beta_star, ...
                1/mu_star - beta_star];
        end

    end
end

% prop_inv(x, y)
% -rel_entr(x, y)