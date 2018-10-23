close all; clear; clc;

snr_array = 1:1:10;
err_max = 100;
N = 100;
ber_array_awgn = zeros(1, length(snr_array));
ber_array_fad = zeros(1, length(snr_array));
demodulated = zeros(1, 3 * N);
enc_p = [0, pi/4, 3*pi/4, pi/2, 7*pi/4, 3*pi/2, pi, 5*pi/4];

for sim_i = 1:length(snr_array)
    snr = snr_array(sim_i);
    npower = 1 / (6 * (10 ^ (snr / 10)));
    err_num_awgn = 0;
    total_num_awgn = 0;
    while err_num_awgn < err_max
        bits = randi([0, 1],[1, N * 3]);
        syms_p = zeros(1, N);
        for i = 1:N
            s = char('0' + bits(i*3-2:i*3));
            s = bin2dec(s);
            syms_p(i) = enc_p(s+1);
        end
        syms = exp(syms_p * 1j);
        noise_i = (sqrt(npower) * randn(1, N));
        noise_q = (sqrt(npower) * randn(1, N));
        noise = noise_i + noise_q * 1j;
        receiver = syms + noise;
        receiver_p = angle(receiver);
        receiver_p = receiver_p + (receiver_p < 0) * 2 * pi;
        for i = 1:N
            s = receiver_p(i);
            r = -1;
            for di = 1:length(enc_p)
                d = enc_p(di);
                upper_bound = d + pi / 8;
                if (s < upper_bound) && (s >= d)
                    r = di;
                else
                    lower_bound = d - pi / 8;
                    if(d == 0)
                        lower_bound = lower_bound + 2 * pi;
                        d = 2 * pi;
                    end
                    if (s >= lower_bound) && (s < d)
                        r = di;
                    end
                end
                if r >= 0
                    da = (dec2bin(r-1) - '0');
                    demodulated((i*3-2):(i*3)) = [zeros(1, 3-length(da)), da];
                    break;
                end
            end
        end
        err = sum(abs(bits - demodulated));
        err_num_awgn = err_num_awgn + err;
        total_num_awgn = total_num_awgn + 3 * N;
    end
    ber_array_awgn(sim_i) = err_num_awgn / total_num_awgn;
    
    err_num_fad = 0;
    total_num_fad = 0;
    while err_num_fad < err_max
        bits = randi([0, 1],[1, N * 3]);
        syms_p = zeros(1, N);
        for i = 1:N
            s = char('0' + bits(i*3-2:i*3));
            s = bin2dec(s);
            syms_p(i) = enc_p(s+1);
        end
        syms = exp(syms_p * 1j);
        fad = (randn(1, N) + randn(1, N) * 1j) / sqrt(2);
        noise_i = (sqrt(npower) * randn(1, N));
        noise_q = (sqrt(npower) * randn(1, N));
        noise = noise_i + noise_q * 1j;
        receiver = syms .* abs(fad) + noise;
        receiver = receiver ./ abs(fad);
        receiver_p = angle(receiver);
        receiver_p = receiver_p + (receiver_p < 0) * 2 * pi;
        for i = 1:N
            s = receiver_p(i);
            r = -1;
            for di = 1:length(enc_p)
                d = enc_p(di);
                upper_bound = d + pi / 8;
                if (s < upper_bound) && (s >= d)
                    r = di;
                else
                    lower_bound = d - pi / 8;
                    if(d == 0)
                        lower_bound = lower_bound + 2 * pi;
                        d = 2 * pi;
                    end
                    if (s >= lower_bound) && (s < d)
                        r = di;
                    end
                end
                if r >= 0
                    da = (dec2bin(r-1) - '0');
                    demodulated((i*3-2):(i*3)) = [zeros(1, 3-length(da)), da];
                    break;
                end
            end
        end
        err = sum(abs(bits - demodulated));
        err_num_fad = err_num_fad + err;
        total_num_fad = total_num_fad + 3 * N;
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
ber_theoretical_awgn = berawgn(snr_array, 'psk', 8, 'nondiff');
semilogy(snr_array, ber_theoretical_awgn, 'g.-', 'linewidth', lw, 'markersize', ms);
ber_theoretical_fad = berfading(snr_array, 'psk', 8, 1);
semilogy(snr_array, ber_theoretical_fad, 'c.-', 'linewidth', lw, 'markersize', ms);
% axis([min(snr_array), max(snr_array), 0.0001, 1]);
xlabel('SNR(dB)');
ylabel('BER');
title('8PSK');
legend('Simulation AWGN', 'Simulation Fading', 'Theoretical AWGN', 'Theoretical Fading');

