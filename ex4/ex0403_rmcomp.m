close all; clear; clc;

snr_array = -1:8;
err_max = 100;
N = 100;
ber_array_m = zeros(1, length(snr_array));
ber_array_r = zeros(1, length(snr_array));
ber_array_g = zeros(1, length(snr_array));

SP = 4;
M = 10;
MP = 2^M - 1;

mseq = idinput(MP, 'prbs')';
mseq = mseq / sqrt(SP);
mseqlen = length(mseq);
spread_m = [repmat(mseq, 1, floor(SP * N / mseqlen)), mseq(1:mod(SP * N, mseqlen))];

rseq = randi([0, 1], [1, SP]);
rseq = 1 - rseq * 2;
rseq = rseq / sqrt(SP);
rseqlen = length(rseq);
spread_r = [repmat(rseq, 1, floor(SP * N / rseqlen)), rseq(1:mod(SP * N, rseqlen))];

gseq = idinput([MP, 2], 'prbs') > 0;
gseq = xor(gseq(:, 1), gseq(:, 2));
gseq = gseq';
gseq = 1 - 2 * gseq;
gseq = gseq / sqrt(SP);
gseqlen = length(gseq);
spread_g = [repmat(gseq, 1, floor(SP * N / gseqlen)), gseq(1:mod(SP * N, gseqlen))];


for sim_i = 1:length(snr_array)
    snr = snr_array(sim_i);
    npower = 1 / (2 * (10 ^ (snr / 10)));
    
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
        receiver = syms + noise;
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
        receiver = syms + noise;
        receiver_t = receiver .* spread;
        receiver_t = reshape(receiver_t, SP, []);
        receiver_t = sum(receiver_t);
        demodulated = double(receiver_t < 0);
        err = sum(abs(bits - demodulated));
        err_num = err_num + err;
        total_num = total_num + N;
    end
    ber_array_g(sim_i) = err_num / total_num;
    
    err_num = 0;
    total_num = 0;
    spread = spread_r;
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
    ber_array_r(sim_i) = err_num / total_num;
end

ber_array_r
ber_array_m
ber_array_g

lw = 2;
ms = 16;
figure;
semilogy(snr_array, ber_array_m, 'r.-', 'linewidth', lw, 'markersize', ms);
hold on;
grid on;
semilogy(snr_array, ber_array_r, 'b.-', 'linewidth', lw, 'markersize', ms);
semilogy(snr_array, ber_array_g, 'g.-', 'linewidth', lw, 'markersize', ms);
ber_theoretical = berawgn(snr_array, 'psk', 2, 'nodiff');
semilogy(snr_array, ber_theoretical, 'c.-', 'linewidth', lw, 'markersize', ms);
% axis([min(snr_array), max(snr_array), 0.0001, 1]);
xlabel('SNR(dB)');
ylabel('BER');
title('DSSS');
legend('M Series', 'Random Series', 'Gold Series', 'Theoretical BPSK');
