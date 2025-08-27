
classdef functions

    methods(Static)

        function rate = rate(sensor)
            SNR = (abs(sensor.H_k).^2) .* c.p_k ./ (c.N0 .* inv_pos(sensor.b_k));
            rate = sensor.b_k .* log(1 + SNR) / log(2);
        end

        function beta_star = best_beta(sensor)
            beta_star = min(...
                (1/c.epsilon) * log((c.f_S / (c.epsilon + sensor.D_k)) * (sensor.lam1 + sensor.lam2)), ...
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

        function T_comp = T_comp(sensor)
            beta_star = functions.best_beta(sensor);
            T_comp = sensor.D_k * (exp(beta_star * c.epsilon) - exp(c.epsilon)) / c.f_S;
        end

        function T_tr = T_tr(sensor)
            r = functions.rate(sensor);
            mu_star = functions.best_mu(sensor);
            T_tr = sensor.D_k * inv_pos(r * mu_star);

        end

        function T_DT = T_DT(sensor)
            T_DT = c.c_k * (sensor.D_k * inv_pos(sensor.f_dt_k));
        end

        function T_total_bs = T_total_bs(sensor)
            T_total_bs = functions.T_comp(sensor) + functions.T_tr(sensor) + functions.T_DT(sensor);
        end

        function T_total_bs = T_total_bs_from_b_f(b, f, sensor)
            sensor.b_k = b;
            sensor.f_dt_k = f;
            T_total_bs = functions.T_comp(sensor) + functions.T_tr(sensor) + functions.T_DT(sensor);
        end

        function g = g(sensor)
            beta_star = functions.best_beta(sensor);
            mu_star = functions.best_mu(sensor);
            g = [1 - beta_star, ...
                1/mu_star - beta_star];
        end

        function return_sensor_list = leader_optimization(sensor_list)
            K = length(sensor_list);
            cvx_begin
                variable T nonnegative;
                variable b(K) nonnegative;
                variable f(K) nonnegative;

                for k = 1:K
                    sensor_list(k).f_dt_k = f(k);
                    sensor_list(k).b_k = b(k);
                end

                minimize(T);
                subject to
                sum(b) <= c.B_total;
                sum(f) <= c.C_DT;
                for k = 1:K
                    T >= functions.T_total_bs(sensor_list(k));
                    T >= 0;
                end
            cvx_end
            for k = 1:K
                sensor_list(k).b_k = b(k);
                sensor_list(k).f_dt_k = f(k);
            end
            return_sensor_list = sensor_list;
        end

        function return_sensor_list = T_DT_optimization(sensor_list)
            K = length(sensor_list);
            cvx_begin
                variable T nonnegative;
                variable f(K) nonnegative;
                
                minimize(T);
                subject to
                    sum(f) <= c.C_DT;
                    for k = 1:K
                        T >= functions.T_DT(sensor_list(k));
                        T >= 0;
                    end
            cvx_end
            for k = 1:K
                    sensor_list(k).f_dt_k = f(k);
            end
            return_sensor_list = sensor_list;
        end

        function return_sensor_list = T_tr_optimization(sensor_list)
            K = length(sensor_list);
            cvx_begin
                cvx_precision low
                cvx_solver sedumi
                variable T nonnegative;
                variable b(K) nonnegative;

                minimize(T);

                subject to
                   for k = 1:K
                    -rel_entr(c.N0 * b(k), c.N0 * b(k) + (abs(sensor_list(k).H_k) ^ 2) * c.p_k) * inv_pos(log(2)) >= sensor_list(k).D_k * inv_pos(functions.best_beta(sensor_list(k)) * T * c.N0);
                   T >= 0;
                   end
                   sum(b) <= c.B_total;
                
                
            cvx_end
            for k = 1:K
                    sensor_list(k).b_k = b(k);
            end
            return_sensor_list = sensor_list;
        end

    end
end

% prop_inv(x, y)
% -rel_entr(x, y)