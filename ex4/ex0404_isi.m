close all; clear; clc;

snr_array = -1:8;
err_max = 100;
N = 10000;
ber_array_n = zeros(1, length(snr_array));
ber_array_o = zeros(1, length(snr_array));
ber_array_m = zeros(1, length(snr_array));
ber_array_g = zeros(1, length(snr_array));

SP = 4;
M = 4;
MP = 2^M - 1;

mseq = idinput(MP, 'prbs')';
mseq = mseq / sqrt(SP);
mseqlen = length(mseq);
spread_m = [repmat(mseq, 1, floor(SP * N / mseqlen)), mseq(1:mod(SP * N, mseqlen))];

gseq = idinput([MP, 2], 'prbs');
gseq = gseq(:, 1) .* gseq(:, 2);
gseq = gseq';
gseq = gseq / sqrt(SP);
gseqlen = length(gseq);
spread_g = [repmat(gseq, 1, floor(SP * N / gseqlen)), gseq(1:mod(SP * N, gseqlen))];

for sim_i = 1:length(snr_array)
    snr = snr_array(sim_i);
    npower = 1 / (2 * (10 ^ (snr / 10)));
    
    if snr <= 8
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
        ber_array_n(sim_i) = err_num / total_num;
    end
    
    err_num = 0;
    total_num = 0;
    while err_num < err_max
        bits = randi([0, 1],[1, N]);
        syms = 1 - 2 * bits;
        noise = sqrt(npower) * randn(1, N);
        fad = abs((randn(1, N) + randn(1, N) * 1j) / sqrt(2));
        receiver = syms .* fad + noise;
        demodulated = double(receiver < 0);
        err = sum(abs(bits - demodulated));
        err_num = err_num + err;
        total_num = total_num + N;
    end
    ber_array_o(sim_i) = err_num / total_num;
    
    err_num = 0;
    total_num = 0;
    spread = spread_m;
    while err_num < err_max
        bits = randi([0, 1],[1, N]);
        syms_t = 1 - 2 * bits;
        syms_t = repmat(syms_t, SP, 1);
        syms_t = reshape(syms_t, 1, []);
        syms = syms_t .* spread;
        noise = sqrt(npower) * randn(1, SP * N);
        fad = abs((randn(1, SP * N) + randn(1, SP * N) * 1j) / sqrt(2));
        receiver = syms .* fad + noise;
        receiver_t = receiver .* spread;
        receiver_t = reshape(receiver_t, SP, []);
        receiver_t = sum(receiver_t);
        demodulated = double(receiver_t < 0);
        err = sum(abs(bits - demodulated));
        err_num = err_num + err;
        total_num = total_num + N;
    end
    ber_array_m(sim_i) = err_num / total_num;
    
    err_num = 0;
    total_num = 0;
    spread = spread_g;
    while err_num < err_max
        bits = randi([0, 1],[1, N]);
        syms_t = 1 - 2 * bits;
        syms_t = repmat(syms_t, SP, 1);
        syms_t = reshape(syms_t, 1, []);
        syms = syms_t .* spread;
        noise = sqrt(npower) * randn(1, SP * N);
        fad = abs((randn(1, SP * N) + randn(1, SP * N) * 1j) / sqrt(2));
        receiver = syms .* fad + noise;
        receiver_t = receiver .* spread;
        receiver_t = reshape(receiver_t, SP, []);
        receiver_t = sum(receiver_t);
        demodulated = double(receiver_t < 0);
        err = sum(abs(bits - demodulated));
        err_num = err_num + err;
        total_num = total_num + N;
    end
    ber_array_g(sim_i) = err_num / total_num;
end

ber_array_n
ber_array_o
ber_array_m
ber_array_g

% figure;
lw = 2;
ms = 16;
semilogy(snr_array, ber_array_n, 'r.-', 'linewidth', lw, 'markersize', ms);
hold on;
grid on;
semilogy(snr_array, ber_array_o, 'b.-', 'linewidth', lw, 'markersize', ms);
semilogy(snr_array, ber_array_m, 'g.-', 'linewidth', lw, 'markersize', ms);
semilogy(snr_array, ber_array_g, 'm.-', 'linewidth', lw, 'markersize', ms);
% axis([min(snr_array), max(snr_array), 0.0001, 1]);
xlabel('SNR(dB)');
ylabel('BER');
title('DSSS ISI');
legend('Normal', 'ISI', 'M Series ISI', 'Gold Series ISI');
