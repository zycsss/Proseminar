function plot_all_vary(t_max, beta_avg, p_k_list, bandwidth_list, beta_max_list)
    arguments (Input)
        t_max
        beta_avg
        p_k_list = c.p_k_list;
        bandwidth_list = c.bandwidth_list;
        beta_max_list = c.beta_max_list;
    end

    tiledlayout(2, 3);
    
    

    name = {'proposed scheme', 'fixed b_k', 'fixed f^{DT}_k', 'fixed b_k and f^{DT}_k', '\beta = 1'};

    field_names = {
                'default', ...
                'fixed_b_k', ...
                'fixed_f_DT_k', ...
                'fixed_both', ...
                'beta_max_1'
            };


    nexttile;
    for i = 1: length(field_names)
        field = field_names{i};
        plot(beta_max_list, t_max.beta.(field))
        hold on
    end
    xlabel('max compresion ratio');
    ylabel('t_{max} [s]');
    legend(name);
    
            
    nexttile;
    for i = 1: length(field_names)
        field = field_names{i};
        semilogy(p_k_list, t_max.p_k.(field));
        hold on
    end
    xlabel('p_k [W]');
    ylabel('t_{max} [s]');
    legend(name);

    nexttile;
    for i = 1: length(field_names)
        field = field_names{i};
        semilogx(bandwidth_list, t_max.B.(field));
        hold on
    end
    xlabel('bandwidth [Hz]');
    ylabel('t_{max} [s]');
    legend(name);

    nexttile;
    for i = 1: length(field_names)
        field = field_names{i};
        plot(beta_max_list, beta_avg.beta.(field))
        hold on
    end
    xlabel('max compresion ratio');
    ylabel('beta_{avg}');
    legend(name);
    
    nexttile;
    for i = 1: length(field_names)
        field = field_names{i};
        plot(p_k_list, beta_avg.p_k.(field));
        hold on
    end
    xlabel('p_k [W]');
    ylabel('beta_{avg}');
    legend(name);
    
    nexttile;
    for i = 1: length(field_names)
        field = field_names{i};
        semilogx(bandwidth_list, beta_avg.B.(field));
        hold on
    end
    xlabel('bandwidth [Hz]');
    ylabel('beta_{avg}');
    legend(name);
end