
% bar_charts
% n = 500;
% seed_list = randi([1, inf], 1, n);
% [b_avg, f_avg, beta_avg, t_avg] = functions.run_many_times(n, seed_list, "type", scenario.D_k);
% plot_bar(b_avg, f_avg, beta_avg, t_avg)

% 3 varied
[t_max, beta_avg] = functions.run_multiple_varies(50);
plot_all_vary(t_max, beta_avg)

