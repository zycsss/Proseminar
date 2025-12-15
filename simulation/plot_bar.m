function plot_bar(b, f, beta, tStruct)
    K = length(b);

    tiledlayout(2, 2);

    nexttile
    bar(1:K, b)
    title('Bandwidth')
    xlabel('sensor');
    ylabel('Bandwidth (Hz)');

    nexttile
    bar(1:K, f)
    title('f_{dt}')
    xlabel('sensor')
    ylabel('op-rate (Hz)');

    nexttile
    bar(1:K, beta)
    title('beta')
    xlabel('sensor')
    ylabel('comp-ratio');
    ylim([1 inf])

    nexttile
    stack_name = {'comp', 'tr', 'DT'};

    time = [tStruct.comp_list; tStruct.tr_list; tStruct.DT_list];
    bar(1:K, time, 'stacked')
    title('Time used for each process on difference sensors')
    legend(stack_name)

    xlabel('sensor')
    ylabel('time (s)');
end