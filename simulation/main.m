
rng("default");

sensor_list = functions.gen_sensor_list_not_uniform();
K = length(sensor_list);

t = [0];

delta_t = 10000;
l = 1;

while delta_t > c.xi && l < c.L_max
    sensor_list = functions.leader_optimization(sensor_list);

    t(end + 1) = functions.T_bs(sensor_list);

    for k = 1:K
        sensor_list(k) = functions.linesearch(sensor_list(k));
    end

    t(end + 1) = functions.T_bs(sensor_list);

    delta_t = abs(t(end) - t(end-1));
    l = l+1;
end

t(1) = [];

sensor_name = categorical({'s_1', 's_2', 's_3', 's_4'});
stack_name = {'comp', 'tr', 'DT'};

time = [functions.T_comp_list(sensor_list); functions.T_tr_list(sensor_list); functions.T_DT_list(sensor_list)];

figure;
bar(sensor_name, time, 'stacked')
legend(stack_name)

xlabel('sensor')
ylabel('time (s)');