
% rng("shuffle");

% sensor_list = functions.gen_sensor_list_not_uniform();

% [b, f, beta, tStruct] = functions.run_one_time(sensor_list);

% t_avg = tStruct;
[b_avg, f_avg, beta_avg, t_avg] = functions.run_many_times(100);
plot_bar(b_avg, f_avg, beta_avg, t_avg)