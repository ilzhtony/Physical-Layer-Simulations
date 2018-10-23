close all; clear; clc;

snr_array = 1:1:10;
err_max = 100;
N = 100;
ber_array_awgn = zeros(1, length(snr_array));
ber_array_fad = zeros(1, length(snr_array));
demodulated = zeros(1, 4 * N);

c1 = 1;
c2 = 3;

enc_coord = [
    c1, -c2;
    c1, -c1;
    c1, c2;
    c1, c1;
    c2, -c2;
    c2, -c1;
    c2, c2;
    c2, c1;
    -c1, -c2;
    -c1, -c1;
    -c1, c2;
    -c1, c1;
    -c2, -c2;
    -c2, -c1;
    -c2, c2;
    -c2, c1;
] / sqrt(10);

for sim_i = 1:length(snr_array)
    snr = snr_array(sim_i);
    npower = 1 / (2 * 4 * (10 ^ (snr / 10)));
    err_num_awgn = 0;
    total_num_awgn = 0;
    while err_num_awgn < err_max
        bits = randi([0, 1],[1, N * 4]);
        syms_i = zeros(1, N);
        syms_q = zeros(1, N);
        for i = 1:N
            s = char('0' + bits(i*4-3:i*4));
            s = bin2dec(s);
            syms_i(i) = enc_coord(s+1, 1);
            syms_q(i) = enc_coord(s+1, 2);
        end
        syms = syms_i + syms_q * 1j;
        noise_i = (sqrt(npower) * randn(1, N));
        noise_q = (sqrt(npower) * randn(1, N));
        noise = noise_i + noise_q * 1j;
        receiver = syms + noise;
        receiver_i = real(receiver);
        receiver_q = imag(receiver);
        for i = 1:N
            g = distance([receiver_i(i), receiver_q(i)], enc_coord);
            s = find(g == min(g));
            da = (dec2bin(s-1) - '0');
            demodulated((i*4-3):(i*4)) = [zeros(1, 4-length(da)), da];
        end
        err = sum(abs(bits - demodulated));
        err_num_awgn = err_num_awgn + err;
        total_num_awgn = total_num_awgn + 4 * N;
    end
    ber_array_awgn(sim_i) = err_num_awgn / total_num_awgn;
    
    err_num_fad = 0;
    total_num_fad = 0;
    while err_num_fad < err_max
        bits = randi([0, 1],[1, N * 4]);
        syms_i = zeros(1, N);
        syms_q = zeros(1, N);
        for i = 1:N
            s = char('0' + bits(i*4-3:i*4));
            s = bin2dec(s);
            syms_i(i) = enc_coord(s+1, 1);
            syms_q(i) = enc_coord(s+1, 2);
        end
        syms = syms_i + syms_q * 1j;
        fad = (randn(1, N) + randn(1, N) * 1j) / sqrt(2);
        noise_i = (sqrt(npower) * randn(1, N));
        noise_q = (sqrt(npower) * randn(1, N));
        noise = noise_i + noise_q * 1j;
        receiver = syms .* abs(fad) + noise;
        receiver = receiver ./ abs(fad);
        receiver_i = real(receiver);
        receiver_q = imag(receiver);
        for i = 1:N
            g = distance([receiver_i(i), receiver_q(i)], enc_coord);
            s = find(g == min(g));
            da = (dec2bin(s-1) - '0');
            demodulated((i*4-3):(i*4)) = [zeros(1, 4-length(da)), da];
        end
        err = sum(abs(bits - demodulated));
        err_num_fad = err_num_fad + err;
        total_num_fad = total_num_fad + 4 * N;
    end
    ber_array_fad(sim_i) = err_num_fad / total_num_fad;
end

[snr_array', ber_array_awgn']
[snr_array', ber_array_fad']

lw = 2;
ms = 16;
figure;
semilogy(snr_array, ber_array_awgn, 'r.-', 'linewidth', lw, 'markersize', ms);
hold on;
grid on;
semilogy(snr_array, ber_array_fad, 'b.-', 'linewidth', lw, 'markersize', ms);
ber_theoretical_awgn = berawgn(snr_array, 'qam', 16);
semilogy(snr_array, ber_theoretical_awgn, 'g.-', 'linewidth', lw, 'markersize', ms);
ber_theoretical_fad = berfading(snr_array, 'qam', 16, 1);
semilogy(snr_array, ber_theoretical_fad, 'c.-', 'linewidth', lw, 'markersize', ms);
[snr_array', ber_theoretical_fad']
% axis([min(snr_array), max(snr_array), 0.0001, 1]);
xlabel('SNR(dB)');
ylabel('BER');
title('16QAM');
legend('Simulation AWGN', 'Simulation Fading', 'Theoretical AWGN', 'Theoretical Fading');

