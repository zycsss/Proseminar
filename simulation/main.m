
rng("default");

K = 1;

sensor_list = functions.gen_sensor_list(K);

t = [0];

delta_t = 10000;
l = 1;

while delta_t > c.xi && l < c.L_max
    % sensor_list = functions.leader_optimization(sensor_list);

    for k = 1:K
        sensor_list(k) = functions.update_lambda(sensor_list(k));
    end

    t(end + 1) = functions.T_bs(sensor_list);

    delta_t = t(end) - t(end-1);
    l = l+1;
    disp({sensor_list.h})
end

t(1) = [];

disp(t)