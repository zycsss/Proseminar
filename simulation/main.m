
% bar_charts
% [b_avg, f_avg, beta_avg, t_avg] = functions.run_many_times(100, "type", scenario.uniform, "beta_max", 2, "fixed_b_k", true);
% plot_bar(b_avg, f_avg, beta_avg, t_avg)

% 3 varied
[t_max, beta_avg] = functions.run_multiple_varies(50);
plot_all_vary(t_max, beta_avg)

% 
% [t_max, beta_avg] = functions.run_one_vary(c.beta_max_list, VaryType.beta_max, 'runs', 1000);
% 
% tiledlayout(1, 1);
% 
% 
% 
%     name = {'proposed scheme', 'fixed b_k', 'fixed f^{DT}_k', 'fixed b_k and f^{DT}_k', '\beta = 1'};
% 
%     field_names = {
%                 'default', ...
%                 'fixed_b_k', ...
%                 'fixed_f_DT_k', ...
%                 'fixed_both', ...
%                 'beta_max_1'
%             };
% 
% 
%     nexttile;
%     for i = 1: length(field_names)
%         field = field_names{i};
%         plot(c.beta_max_list, t_max.(field))
%         hold on
%     end
%     xlabel('max compresion ratio');
%     ylabel('t_{max} [s]');
%     legend(name);
% 