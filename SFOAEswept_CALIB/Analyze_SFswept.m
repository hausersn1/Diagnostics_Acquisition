% DPOAE swept Analysis
% Author: Samantha Hauser
% Created: May 2023
% Last Updated: Sept 16, 2023
% Purpose:
% Helpful info: Ensure a datafile is loaded

% For quick visualization after data collection.

%%%%%%%%% Set these parameters %%%%%%%%%%%%%%%%%%
windowdur = 0.040; % 40ms in paper
offsetwin = 0.0; % 20ms in paper
npoints = 512; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpfilter = 1; 

% to be compatible with old data
if exist('data', 'var')
    stim = data.stim;
    subj = data.info.subj.ID;
    ear = data.info.subj.ear;
    SFOAEtrials = data.resp.ProbeBuffs + data.resp.SuppBuffs - data.resp.BothBuffs; 
else
    SFOAEtrials = stim.ProbeBuffs + stim.SuppBuffs - stim.BothBuffs; 
end

%% Set variables needed from stim.
phiProbe_inst = 2*pi*stim.phiProbe_inst;
t = stim.t;

% downward vs upward sweeps
if stim.speed < 0
    f1 = stim.fmax;
    f2 = stim.fmin;
else
    f1 = stim.fmin;
    f2 = stim.fmax;
end

% set freq we're testing and the timepoints when they happen.
if abs(stim.speed) < 20 %linear sweep
    testfreq = 2 .^ linspace(log2(f1), log2(f2), npoints);
    t_freq = log2(testfreq/f1)/stim.speed + stim.buffdur;
else % log sweep
    testfreq = linspace(f1, f2, npoints);
    t_freq = (testfreq-f1)/stim.speed + stim.buffdur;
end

%duration changes w/ frequency
durs = .038*(2.^(-0.3*t_freq)-1)/ (-0.3*log(2)) + 0.038;

%% Artifact rejection
trials = size(SFOAEtrials,1);

% high pass filter the response (can also be done on ER10X hardware)
filtcutoff = 300;
b = fir1(1000, filtcutoff*2/stim.Fs, 'high');
SFOAEtrials= filtfilt(b, 1, SFOAEtrials')';

% Set empty matricies for next steps
coeffs = zeros(npoints, 6);
a_temp = zeros(trials, npoints);
b_temp = zeros(trials, npoints);

% Least Squares fit of SF Only for AR
for x = 1:trials
    SFOAE = SFOAEtrials(x, :);
    fprintf(1, 'Checking trial %d / %d for artifact\n', x, trials);
    
    for k = 1:npoints
        windowdur = durs(k);
        win = find( (t > (t_freq(k) - windowdur/2)) & ...
            (t < (t_freq(k) + windowdur/2)));
        taper = hanning(numel(win))';
        model_sf = [cos(phiProbe_inst(win)) .* taper;
            -sin(phiProbe_inst(win)) .* taper];
        resp = SFOAE(win) .* taper;
        coeffs(k, 1:2) = model_sf' \ resp';
    end
    a_temp(x,:) = coeffs(:, 1);
    b_temp(x,:) = coeffs(:, 2);
end

oae = abs(complex(a_temp, b_temp));

median_oae = median(oae);
std_oae = std(oae);
resp_AR = SFOAEtrials;
for j = 1:trials
    for k = 1:npoints
        if oae(j,k) > median_oae(1,k) + 4*std_oae(1,k)
            win = find( (t > (t_freq(k) - durs(k).*.1)) & ...
                (t < (t_freq(k) + durs(k).*.1)));
            resp_AR(j,win) = NaN;
        end
    end
end

SFOAE = mean(resp_AR, 'omitNaN'); % mean SFOAE after artifact rejection

%% LSF Analysis

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2);
coeffs_n = zeros(npoints, 2);
tau = zeros(npoints, 1);
coeffs_noise = zeros(npoints,8);

durs = .038*(2.^(-0.3*t_freq)-1)/ (-0.3*log(2)) + 0.038;

% Generate model of chirp and test against response
for k = 1:npoints
    fprintf(1, 'Running window %d / %d\n', k, (npoints));
    %windowdur = durs(k);
    win = find( (t > (t_freq(k)-windowdur/2)) & ...
        (t < (t_freq(k)+windowdur/2)));
    taper = hanning(numel(win))';
    
    % SF probe frequency
    model = [cos(phiProbe_inst(win)) .* taper;
        -sin(phiProbe_inst(win)) .* taper];
    
    % nearby frequencies for nf calculation
    model_noise = [cos(1.1*phiProbe_inst(win)) .* taper;
        -sin(1.1*phiProbe_inst(win)) .* taper;
        cos(1.12*phiProbe_inst(win)) .* taper;
        -sin(1.12*phiProbe_inst(win)) .* taper;
        cos(1.14*phiProbe_inst(win)) .* taper;
        -sin(1.14*phiProbe_inst(win)) .* taper;
        cos(1.16*phiProbe_inst(win)) .* taper;
        -sin(1.16*phiProbe_inst(win)) .* taper];
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 2);
    resid = zeros(maxoffset, 1);
    coeff_noise = zeros(maxoffset, 8);
    
    for offset = 0:maxoffset
        resp = SFOAE(win+offset) .* taper;
        
        coeff(offset + 1, :) = model' \ resp';
        coeff_noise(offset + 1, :) = model_noise' \ resp';
        
        resid(offset +1) = sum( (resp - coeff(offset+1, :) * model).^2);
    end
    
    [~, ind] = min(resid);
    
    coeffs(k, :) = coeff(ind, :);
    coeffs_noise(k,:) = coeff_noise(ind,:);
    
    tau(k) = (ind - 1) * (1/stim.Fs); % delay in sec
    
end

%% Amplitude and delay calculations
a = coeffs(:, 1);
b = coeffs(:, 2);

phi = tau.*testfreq'; % cycles (from delay/offset)
phasor = exp(-1j * phi* 2 * pi);

% for noise
noise = zeros(npoints,4);
for i = 1:2:8
    noise(:,ceil(i/2)) = complex(coeffs_noise(:,i), coeffs_noise(:,i+1));
end

oae_complex = complex(a, b).*phasor;
noise_complex = mean(noise,2);
VtoSPL = stim.VoltageToPascal.* stim.PascalToLinearSPL;

%% Plot resulting figure
figure(3); hold on; 
plot(testfreq/1000, db(abs(oae_complex).*VtoSPL), 'linew', 2,'Color', [0.4940 0.1840 0.5560]);
hold on;
plot(testfreq/1000, db(abs(noise_complex).*VtoSPL), '--', 'linew', 2, 'Color', [0.6350 0.0780 0.1840]);
title('SFOAE', 'FontSize', 14 )
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-50, 40])
xticks([.5, 1, 2, 4, 8, 16])
xlabel('Frequency (kHz)','FontWeight','bold')
ylabel('Amplitude (dB SPL)','FontWeight','bold')
legend('SFOAE', 'NF')
drawnow; 
