
% rng("default");

sensor_list = functions.gen_sensor_list_not_uniform();

[b, f, beta, tStruct] = functions.run_one_time(sensor_list);
t_avg = functions.run_many_times(functions.gen_sensor_list_not_uniform(), 50);
plot_bar(b, f, beta, tStruct)