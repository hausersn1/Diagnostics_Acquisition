% Plots from a pre-loaded TEOAE data structure
freq = data.resp.freq;
Resp = data.resp.output_Pa_per_20uPa;
NoiseFloor = data.resp.noise_Pa_per_20uPa; 

figure(1);
hold on;
plot(freq*1e-3, db(abs(Resp)), 'linew', 2);
ylabel('Response (dB SPL)', 'FontSize', 16);

uplim = max(db(abs(Resp)));
hold on;
semilogx(freq*1e-3, db(abs(NoiseFloor)), 'linew', 2);
xlabel('Frequency (kHz)', 'FontSize', 16);
legend('TEOAE', 'NoiseFloor');
xlim([0.4, 16]);
ticks = [0.5, 1, 2, 4, 8, 16];
set(gca, 'XTick', ticks, 'FontSize', 14, 'xscale', 'log');
ylim([-60, uplim + 5]);