function [t, beta_avg] = run_one_vary_test(vary_list, vary_type, options)
    arguments
        vary_list;
        vary_type;
        options.runs = 100;
    end

    K = length(vary_list);
    if K == 0
        t = null;
        beta_avg = null;
        return
    end
    
    beta_max_list = repelem(c.beta_max, K);
    p_k_list = repelem(c.p_k, K);
    bandwidth_list = repelem(c.B_total/4, K);

    switch vary_type
        case VaryType.beta_max
            beta_max_list = vary_list;
        case VaryType.B
            bandwidth_list = vary_list;
        case VaryType.p_k
            p_k_list = vary_list;
    end

    K = length(vary_list);

    % 1. Define all the field names in a cell array
    field_names = {
        'default', ...
        'fixed_b_k', ...
        'fixed_f_DT_k', ...
        'fixed_both', ...
        'beta_max_1'
    };

    % 2. Define the initial value once
    init_array = zeros(1, K);

    % 3. Initialize the structs (good practice)
    beta_avg = struct();
    t = struct();

    % 4. Loop through the field names and assign the value
    for i = 1:length(field_names)
        field = field_names{i};
        beta_avg.(field) = init_array;
        t.(field) = init_array;
    end
    seed_list = rand(1,options.runs);

    for k = 1:K
        [~, ~, beta_avg_k, t_avg] = functions.run_many_times(options.runs, seed_list, "beta_max", beta_max_list(k), "bandwidth", bandwidth_list(k), "p_k", p_k_list(k), "fixed_b_k", true);
        beta_avg.fixed_b_k(k) = mean(beta_avg_k);
        t.fixed_b_k(k) = t_avg.total;
    end
    
end


[t_max, beta_avg] = run_one_vary_test(c.beta_max_list, VaryType.beta_max, 'runs', 50);

plot(c.beta_max_list, t_max.fixed_b_k)