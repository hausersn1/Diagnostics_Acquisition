% SFOAEtrials = stim.ProbeBuffs + stim.SuppBuffs - stim.BothBuffs;
% 
% SFOAE = median(SFOAEtrials, 1) * stim.VoltageToPascal ...
%     * stim.PascalToLinearSPL;
% 
% Noisetrials = SFOAEtrials;
% Noisetrials(2:2:end) = Noisetrials(2:2:end) * -1;
% 
% NOISE = median(Noisetrials, 1)  * stim.VoltageToPascal ...
%     * stim.PascalToLinearSPL;

SFOAE = cos(2*pi*stim.phiProbe_inst); 
 NOISE = cos(2*pi*stim.phiProbe_inst); 
 
 

figure; 
[Pxx, ~] = pmtm(SFOAE, 1.5, [], stim.Fs);
[Pxx_n, f] = pmtm(NOISE, 1.5, [], stim.Fs);
plot(f, pow2db(Pxx), 'linew', 2);
hold on;
plot(f, pow2db(Pxx_n), 'linew', 2);
xlim([stim.fmin, stim.fmax]);
