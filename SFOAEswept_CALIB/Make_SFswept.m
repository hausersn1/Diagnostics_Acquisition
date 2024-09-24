function stim = Make_SFswept()

% Stimulus Parameters
stim.fmin = 500; %stim.cf/sqrt(2); % 1/2 octave below
stim.fmax = 16000; %stim.cf*sqrt(2); % 1/2 octave above 
stim.speed = -1; % oct/sec downsweep
stim.scale = 'log'; % octave sweep, alt is 'linear'
stim.diff = 50; % Hz (Fprobe - 50 = Fsupp; Probe is higher)
stim.buffdur = 0.25; %seconds; for either side of sweep 
stim.Fs = 48828.125;

stim.drop_Probe = 54; % 75; % for 40dB probe, was 60
stim.drop_Supp = 34; % 55 for 60dB suppressor, was 40

% Trial Parameters
stim.ThrowAway = 0;
stim.SNRcriterion = 6; 
stim.minTrials = 12; 
stim.maxTrials = 50; 

% Add useful info to structure
gain = 30; 
stim.mic_sens = 50e-3; % mV/Pa
stim.mic_gain = db2mag(gain + 6); % +6 for balanced cable
stim.P_ref = 20e-6;
stim.DR_onesided = 1;
stim.VoltageToPascal = 1 / (stim.DR_onesided * stim.mic_gain * stim.mic_sens);
stim.PascalToLinearSPL = 1 /  stim.P_ref;

%% For live analysis
stim.windowdur = 0.06;
stim.npoints = 128; 

%% Create the stimulus
buffdur = stim.buffdur; 
Fs = stim.Fs;

if stim.speed < 0 %downsweep
    f_start = stim.fmax; 
    f_end = stim.fmin; 
    stim.nearfreqs = [1.10, 1.12, 1.14, 1.16];
else 
    f_start = stim.fmin; 
    f_end = stim.fmax; 
    stim.nearfreqs = [.90, .88, .86, .84];
end 

if strcmp(stim.scale, 'log')
    dur = log2(stim.fmax/stim.fmin) / abs(stim.speed) + (2*buffdur);
else
    dur = abs(f_start - f_end) / abs(stim.speed) + (2*buffdur);
end

t = 0: (1/Fs): (dur - 1/Fs);
stim.t = t;

buffinst1 = find(t < buffdur, 1, 'last');
buffinst2 = find(t > (dur-buffdur) , 1, 'first');

if strcmp(stim.scale, 'log')
    % Create probe
    start_probe = f_start*t(1:buffinst1);
    buffdur_exact = t(buffinst1);
    phiProbe_inst = f_start * (2.^( (t-buffdur_exact) * stim.speed) - 1) / (stim.speed * log(2)) + start_probe(end);
    end_probe = f_end*t(1:(length(t)-buffinst2+1)) + phiProbe_inst(buffinst2);
    phiProbe_inst(1:buffinst1) = start_probe;
    phiProbe_inst(buffinst2:end) = end_probe;
else % linear sweep
    start_probe = f_start*t(1:buffinst1);
    buffdur_exact = t(buffinst1);
    phiProbe_inst = f_start*(t-buffdur_exact) + stim.speed*((t-buffdur_exact).^2)/2 + start_probe(end); % Cycles
    end_probe = f_end*t(1:(length(t)-buffinst2+1)) + phiProbe_inst(buffinst2);
    phiProbe_inst(1:buffinst1) = start_probe;
    phiProbe_inst(buffinst2:end) = end_probe;
end

% Create suppressor
phiSupp_inst = phiProbe_inst - stim.diff*t;

% Save stim
stim.yProbe = scaleSound(rampsound(cos(2 * pi * phiProbe_inst), stim.Fs, 0.005));
stim.ySupp = scaleSound(rampsound(cos(2 * pi * phiSupp_inst), stim.Fs, 0.005));
stim.phiProbe_inst = phiProbe_inst;
stim.phiSupp_inst = phiSupp_inst;
