close all; clear; clc;

snr_array = 1:1:10;
err_max = 50;
N = 100;
ber_array_a = zeros(1, length(snr_array));
ber_array_b = zeros(1, length(snr_array));

for sim_i = 1:length(snr_array)
    snr = snr_array(sim_i);
    npower = 1 / (2 * (10 ^ (snr / 10)));
    err_num = 0;
    total_num = 0;
    while err_num < err_max
        bits = randi([0, 1],[1, N]);
        syms = 1 - 2 * bits;
        noise = sqrt(npower) * randn(1, N);
        receiver = syms + noise;
        demodulated = double(receiver < 0);
        err = sum(abs(bits - demodulated));
        err_num = err_num + err;
        total_num = total_num + N;
    end
    ber_array_a(sim_i) = err_num / total_num;
    
    err_num = 0;
    total_num = 0;
    while err_num < err_max
        bits = double(myrnd(N) > 0.5);
        syms = 1 - 2 * bits;
        noise = sqrt(npower) * mygaussrnd(1, N);
        receiver = syms + noise;
        demodulated = double(receiver < 0);
        err = sum(abs(bits - demodulated));
        err_num = err_num + err;
        total_num = total_num + N;
    end
    ber_array_b(sim_i) = err_num / total_num;
end

[snr_array', ber_array_a']
[snr_array', ber_array_b']

lw = 2;
ms = 16;
figure;
semilogy(snr_array, ber_array_a, 'r.-', 'linewidth', lw, 'markersize', ms);
hold on;
grid on;
semilogy(snr_array, ber_array_b, 'b.-', 'linewidth', lw, 'markersize', ms);
ber_theoretical = berawgn(snr_array, 'psk', 2, 'nodiff');
semilogy(snr_array, ber_theoretical, 'g.-', 'linewidth', lw, 'markersize', ms);
axis([min(snr_array), max(snr_array), 0.0001, 1]);
xlabel('SNR(dB)');
ylabel('BER');
title('BPSK');
legend('Simulation - Matlab Random', 'Simulation - My Random', 'Theoretical');
