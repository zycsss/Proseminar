classdef c
    properties (Constant = true)
        N0 = 1e-12;
        epsilon = 3.5;
        kappa = 10.0;
        alpha = 2.0;
        sigma_k = 0.8;
        A0 = 10^(-30/10);
        beta_max = 3.0;
        f_S = 200e3; 
        p_k = 31.6e-3;
        B_total = 100e3;   
        C_DT = 50e6;
        c_k = 1e4;
        L_max = 300;
        xi = 1e-6;
        d_k = 378;
        D_k = 4e3;
        h0 = 1e-2;
        h_factor = 2;
        c_linesearch = 1e-4;
        beta_linesearch = 0.5;
        beta_max_list = 1:0.25:3;
        p_k_list = linspace(1.6e-3, 100e-3, 10);
        bandwidth_list = logspace(2, 5, 10);
    end
end