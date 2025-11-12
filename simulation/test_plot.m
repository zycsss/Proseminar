function test_plot(beta_max_list, t_max, b_avg)
legend_names = {'proposed scheme', 'fixed b_k', 'fixed f^{DT}_k', 'fixed b_k and f^{DT}_k', '\beta = 1'};

    field_names = {
                'default', ...
                'fixed_b_k', ...
                'fixed_f_DT_k', ...
                'fixed_both', ...
                'beta_max_1'
            };


    for i = 1: length(field_names)
        field = field_names{i};
        plot(beta_max_list, t_max.(field))
        hold on
    end
    xlabel('max compresion ratio');
    ylabel('beta_{avg} [s]');
    legend(legend_names);
end