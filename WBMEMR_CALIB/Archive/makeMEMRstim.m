function stim = makeMEMRstim(rawstim)

% rawstim (structure) should contain fields clickwin, noisewin,
% Fs, VtoPforH, rampdur

if ~exist('rawstim', 'var')
    stim.clickwin = 41.92; % ms
    stim.noisewin = 120; % ms
    stim.Fs = 48828.125;
    stim.rampdur = 5; %ms
    % stim.fc = 6500; % 3 - 10 Kz
    stim.fc = 4500; % 1 - 8 kHz
    % stim.bw  = 7000; % along with fc, this decides the band
    stim.bw = 8400;
    stim.noiseramp = 25;
    stim.nreps = 7; % How many reps per trial
else
    stim = rawstim;
end

% Additional stimulation parameters
stim.clickatt = 30;
stim.noiseatt = 40:-10:0; %Note makeNBNoiseFFT returns RMS of -20 dB re: 1
stim.ThrowAway = 0;
stim.Averages = 32;
stim.pad = 256; % Number of samples extra to read in after stim ends
stim.nLevels = numel(stim.noiseatt);

clicksamps = ceil(stim.clickwin * 1e-3  * stim.Fs);


template = makeNBNoiseFFT(stim.bw, stim.fc,...
    stim.noisewin * 1e-3, stim.Fs, stim.noiseramp * 1e-3, 0);

tokenlength = numel(template);
totalsamps = clicksamps + tokenlength + stim.pad;
stim.noise = zeros(stim.nLevels, stim.ThrowAway + stim.Averages, totalsamps);

for L = 1:stim.nLevels
    for m =  1: (stim.ThrowAway + stim.Averages)
        stim.noise(L, m, clicksamps + (1:tokenlength)) = ...
            makeNBNoiseFFT(stim.bw, stim.fc,...
            stim.noisewin * 1e-3, stim.Fs,...
            stim.noiseramp * 1e-3, 0);
        
    end
end



nsampsclick = 5;
initbuff =  2;
stim.click = zeros(1, totalsamps);
stim.click(initbuff + (1:nsampsclick)) = 0.95;
stim.t = (0:(totalsamps - 1)) / stim.Fs;





