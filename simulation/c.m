classdef c
    properties (Constant = true)
        N0 = 1e-12;
        epsilon = 3.5;
        kappa = 10.0;
        alpha = 2.0;
        sigma_k = 1.0;
        A0 = 10^(-30/10);
        beta_max = 3.0;
        f_S = 200e3; 
        p_k = 10^((15-30)/10);
        B_total = 100e3;   
        C_DT = 50e6;
        c_k = 1e4;
        L_max = 300;
        xi = 1e-6;
        d_k = 378;
        D_k = 4e3;
    end
end