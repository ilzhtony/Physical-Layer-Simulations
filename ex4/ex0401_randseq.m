close all; clear; clc;

snr_array = 1:1:10;
err_max = 50;
N = 100;
ber_array = zeros(1, length(snr_array));

SP = 4;
seq = randi([0, 1], [1, SP]);
seq = 1 - seq * 2;
seq = seq / sqrt(SP);
spread = repmat(seq, 1, N);

for sim_i = 1:length(snr_array)
    snr = snr_array(sim_i);
    err_num = 0;
    total_num = 0;
    npower = 1 / (2 * (10 ^ (snr / 10)));
    while err_num < err_max
        bits = randi([0, 1],[1, N]);
        syms_t = 1 - 2 * bits;
        syms_t = repmat(syms_t, SP, 1);
        syms_t = reshape(syms_t, 1, []);
        syms = syms_t .* spread;
        noise = sqrt(npower) * randn(1, SP * N);
        receiver = syms + noise;
        receiver_t = receiver .* spread;
        receiver_t = reshape(receiver_t, SP, []);
        receiver_t = sum(receiver_t);
        demodulated = double(receiver_t < 0);
        err = sum(abs(bits - demodulated));
        err_num = err_num + err;
        total_num = total_num + N;
    end
    ber_array(sim_i) = err_num / total_num;
end

[snr_array', ber_array']

lw = 2;
ms = 16;
figure;
semilogy(snr_array, ber_array, 'r.-', 'linewidth', lw, 'markersize', ms);
hold on;
grid on;
ber_theoretical = berawgn(snr_array, 'psk', 2, 'nodiff');
semilogy(snr_array, ber_theoretical, 'b.-', 'linewidth', lw, 'markersize', ms);
% axis([min(snr_array), max(snr_array), 0.0001, 1]);
xlabel('SNR(dB)');
ylabel('BER');
legend('Simulation', 'Theoretical');
