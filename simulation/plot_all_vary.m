function plot_all_vary(t_max, beta_avg, p_k_list, bandwidth_list, beta_max_list)
    arguments (Input)
        t_max
        beta_avg
        p_k_list = c.p_k_list;
        bandwidth_list = c.bandwidth_list;
        beta_max_list = c.beta_max_list;
    end

    tiledlayout(2, 3);

    name = {'proposed scheme', '\beta = 1'};

    nexttile;
    plot(beta_max_list, t_max.beta, 'o--');
    hold on;
    plot(beta_max_list, repelem(t_max.beta(1), length(t_max.beta))', '.--');
    xlabel('max compresion ratio');
    ylabel('t_{max} [s]');
    legend(name);
            
    nexttile;
    semilogy(p_k_list, t_max.p_k, 'o--');
    xlabel('p_k [W]');
    ylabel('t_{max} [s]');
    legend(name);

    nexttile;
    semilogx(bandwidth_list, t_max.B, 'o--');
    xlabel('bandwidth [Hz]');
    ylabel('t_{max} [s]');
    legend(name);

    nexttile;
    plot(beta_max_list, beta_avg.beta, 'o--');
    hold on;
    plot(beta_max_list, repelem(beta_avg.beta(1), length(beta_avg.beta))', '.--');
    xlabel('max compresion ratio');
    ylabel('beta_{avg} [s]');
    legend(name);
    
    nexttile;
    plot(p_k_list, beta_avg.p_k, 'o--');
    xlabel('p_k [W]');
    ylabel('beta_{avg}');
    legend(name);
    
    nexttile;
    semilogx(bandwidth_list, beta_avg.B, 'o--');
    xlabel('bandwidth [Hz]');
    ylabel('beta_{avg}');
    legend(name);
end