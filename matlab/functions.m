
classdef functions

    methods(Static)

        function rate = rate(sensor)
            SNR = (abs(sensor.H_k).^2) .* c.p_k ./ (c.N0 .* sensor.b_k);
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
            T_tr = sensor.D_k / (r * mu_star);

        end

        function T_DT = T_DT(sensor)
            T_DT = c.c_k * (sensor.D_k / sensor.f_dt_k);
        end

        function T_total_bs = T_total_bs(sensor)
            T_total_bs = functions.T_comp(sensor) + functions.T_tr(sensor) + functions.T_DT(sensor);
        end

        function g = g(sensor)
            beta_star = functions.best_beta(sensor);
            mu_star = functions.best_mu(sensor);
            g = [1 - beta_star, ...
                1/mu_star - beta_star];
        end

        function return_sensor_list = T_DT_optimization(sensor_list)
            K = length(sensor_list);
            D = zeros(K,1);
            for k=1:K
                D(k) = sensor_list(k).D_k;
            end

            soft = @(x) (exp(x - max(x)) ./ sum(exp(x - max(x)))) * c.C_DT;

            function v = objective(x)
                f = soft(x);
                T = D .* c.c_k ./ f;
                v = max(T);
            end

            x0 = zeros(K,1); opts = optimset('Display','off');
            best_val = inf; best_x = x0;

            for r = 1:6
                xr = (r==1) * x0 + (r>1) * (0.3 * randn(K,1));
                [xt, vt] = fminsearch(@objective, xr, opts);
                if vt < best_val, best_val = vt; best_x = xt; end
            end

            f_opt = soft(best_x);

            for k = 1:K
                    sensor_list(k).f_dt_k = f_opt(k);
            end
            return_sensor_list = sensor_list;
        end

        function return_sensor_list = T_tr_optimization(sensor_list)
            K = length(sensor_list);
            D = zeros(K,1); beta = zeros(K,1); r_k_const = zeros(K,1);
            for k=1:K
                D(k) = sensor_list(k).D_k;
                beta(k) = functions.best_beta(sensor_list(k));
                r_k_const(k) = abs(sensor_list(k).H_k)^2 .* c.p_k ./ c.N0;
            end
            rate_fun = @(b,ck) b .* log(1 + ck ./ max(b,1e-12));
            
            soft = @(x) (exp(x - max(x)) ./ sum(exp(x - max(x)))) * c.B_total;
            function v = objective(x)
                b = soft(x);
                T = D ./ (beta .* rate_fun(b, r_k_const(k)));
                v = max(T);
            end

            x0 = zeros(K,1); opts = optimset('Display','off');
            best_val = inf; best_x = x0;

            for r = 1:6
                xr = (r==1) * x0 + (r>1) * (0.3 * randn(K,1));
                [xt, vt] = fminsearch(@objective, xr, opts);
                if vt < best_val, best_val = vt; best_x = xt; end
            end

            b_opt = soft(best_x);
            for k = 1:K
                sensor_list(k).b_k = b_opt(k);
            end
            return_sensor_list = sensor_list;
        end

    end
end