function plot_bar(b, f, beta, tStruct)
    K = length(b);
    figure(1)
    bar(1:K, b)
    title('Bandwidth')
    xlabel('sensor');
    ylabel('Bandwidth (Hz)');

    figure(2)
    bar(1:K, f)
    title('f_{dt}')
    xlabel('sensor')
    ylabel('op-rate (Hz)');

    figure(3)
    bar(1:K, beta)
    title('beta')
    xlabel('sensor')
    ylabel('comp-ratio');

    figure(4)
    stack_name = {'comp', 'tr', 'DT'};

    time = [tStruct.comp_list; tStruct.tr_list; tStruct.DT_list];
    bar(1:K, time, 'stacked')
    legend(stack_name)

    xlabel('sensor')
    ylabel('time (s)');
end