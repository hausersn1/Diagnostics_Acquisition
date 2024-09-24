%% SFOAEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Samantha Hauser, AuD
% Modified from: Hari Bharadwaj, PhD (SNAP Lab)
% Created: November 2021
% Last revision: 6-Sep-2023 (added artifact checking)
%               24-Sep-2024 (correctly applied FPL calib)
%
% References:
%   SNR based endpoint: Abdala, C., Luo, P. & Shera, C.A. Characterizing
%       the Relationship Between Reflection and Distortion Otoacoustic
%       Emissions in Normal-Hearing Adults. JARO 23, 647–664 (2022).
%       https://doi.org/10.1007/s10162-022-00857-z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data Storage
% Measure-General info
info.measure = 'SFOAEswept';
info.version = 'Auto_v02';

% v01 - original version
% v02 - correctly applies FPL ear calibration

% Visit info
if exist('C:\Experiments\Sam\current_visit.mat','file')
    load('C:\Experiments\Sam\current_visit.mat', 'visit')
    ask = questdlg(sprintf('Is this subject %s?', visit.subj.ID), ...
        'Check Subject', 'Yes', 'No', 'No');
else
    ask = 'No';
end

if strcmp(ask, 'No')
    cd ..
    startVisit
    cd('SFOAEswept_CALIB')
end

subj = visit.subj;
info.room = visit.room;
info.univ = visit.univ;
info.researcher = visit.researcher;

% Get ear info
ear = questdlg('Which ear?', 'Ear', 'L', 'R', 'R');
earCode = ear; 
subj.ear = ear; 

% Get date/time
datetag = datestr(clock);
info.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '-';

% Make directory to save results
paraDir = 'C:\Experiments\Sam\SFOAEswept_CALIB\Results\';
respDir = strcat(paraDir,filesep,visit.subj.ID,filesep);
addpath(genpath(paraDir));
if(~exist(respDir,'dir'))
    mkdir(respDir);
end

fname = strcat(respDir, info.measure, '_', ...
    subj.ID, '_', subj.ear, '_', datetag, '.mat');

% Check to make sure there are EarCal files for both drivers
calib1_file = ['C:\Experiments\FPLclick\EARCAL\', subj.ID, '\', 'Calib_Ph1ER-10X_', subj.ID, earCode,'Ear_', info.date(1:11), '*.mat']; 
calib2_file = ['C:\Experiments\FPLclick\EARCAL\', subj.ID, '\', 'Calib_Ph2ER-10X_', subj.ID, earCode,'Ear_', info.date(1:11), '*.mat']; 

if (isempty(dir(calib1_file)) | isempty(dir(calib2_file)))
    fprintf('------No Ear Calibration Found -- Go run calibEar!!------\n')
    return
else
    % load the calibration files in
    get_calib_1 = uigetfile(calib1_file, 'Get Driver 1 calib file');
    get_calib_2 = uigetfile(calib2_file, 'Get Driver 2 calib file');
    calib_1 = load(get_calib_1);
    calib_2 = load(get_calib_2);
end

%% Run Test
tic;

try
    
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    initializeER10X;
    % initializeER10X_300Hz_Highpass;
    
    % Initializing TDT and specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    % Get stimulus structure; Change to _linear for linear sweep
    stim = Make_SFswept;
    
    % Set live analysis parameters
    windowdur = stim.windowdur;
    SNRcriterion = stim.SNRcriterion;
    maxTrials = stim.maxTrials;
    minTrials = stim.minTrials;
    
    phiProbe_inst = stim.phiProbe_inst*2*pi;
    t = stim.t;
    npoints = stim.npoints; 
    nearfreqs = stim.nearfreqs;
    VtoSPL = stim.VoltageToPascal .* stim.PascalToLinearSPL;
    
    edges = 2 .^ linspace(log2(stim.fmin), log2(stim.fmax), 21);
    bandEdges = edges(2:2:end-1);
    centerFreqs = edges(3:2:end-2);
    
    if stim.speed < 0
        f_start = stim.fmax;
        f_end = stim.fmin;
    else
        f_start = stim.fmin;
        f_end = stim.fmax;
    end
    
    testfreq = 2 .^ linspace(log2(f_start), log2(f_end), npoints); 
     
    if strcmp(stim.scale, 'log')
        t_freq = log2(testfreq/f_start)/stim.speed + stim.buffdur;
    else
        t_freq = (testfreq-f_start)/stim.speed + stim.buffdur;
    end
    
    k = 0;
    doneWithTrials = 0;
    figure;
    
    % Make arrays to store measured mic outputs
    ProbeBuffs = zeros(maxTrials, numel(stim.yProbe));
    SuppBuffs = zeros(maxTrials, numel(stim.yProbe));
    BothBuffs = zeros(maxTrials, numel(stim.yProbe));
    flip = -1;

    % generate each filter 
    h1 = flatFPLfilter(calib_1.calib); 
    h2 = flatFPLfilter(calib_2.calib); 
    
    % get the fplh value (FPL1k) at 1000 hZ for each calib file
    FPL1k_1 = interp1(calib_1.calib.freq, db(abs(calib_1.calib.Pfor)), 1000); 
    FPL1k_2 = interp1(calib_2.calib.freq, db(abs(calib_2.calib.Pfor)), 1000); 
    
    % Filter the stimulus and drop by 30 dB to prevent clipping (will
    % adjust in attenuation (i.e. drop_f1 and drop_f2)
    filt_yProbe = filter(h1, 1, stim.yProbe) * db2mag(FPL1k_1) * db2mag(-30); 
    filt_ySupp = filter(h2, 1, stim.ySupp)* db2mag(FPL1k_2) * db2mag(-30); 
    
    while doneWithTrials == 0
        k = k+1;
        k_kept = k - stim.ThrowAway;
        
        % alternate phase of the suppressor
        flip = flip .* -1;
        delayComp = 1 + 128; % Always 1, plus 128 for filter 
        
        % Do probe only
        dropSupp = 120;
        dropProbe = stim.drop_Probe - 30;
        buffdata = zeros(2, numel(filt_yProbe));
        buffdata(1, :) = filt_yProbe;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        if k > stim.ThrowAway
            ProbeBuffs(k_kept,  :) = vins;
        end
        
        WaitSecs(0.1);
        
        % Do suppressor only
        dropProbe = 120;
        dropSupp = stim.drop_Supp - 30;
        buffdata = zeros(2, numel(filt_ySupp));
        buffdata(2, :) = flip.*filt_ySupp;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        if k > stim.ThrowAway
            SuppBuffs(k_kept,  :) = vins;
        end
        
        WaitSecs(0.1);
        
        % Do both
        dropProbe = stim.drop_Probe - 30;
        dropSupp = stim.drop_Supp - 30;
        buffdata = zeros(2, numel(filt_yProbe));
        buffdata(1, :) = filt_yProbe;
        buffdata(2, :) = flip.*filt_ySupp;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        if k > stim.ThrowAway
            BothBuffs(k_kept,  :) = vins;
        end
        
        WaitSecs(0.1);
        
        if k > stim.ThrowAway
            % Set empty matricies for next steps
            coeffs_resp = zeros(npoints, 2);
            coeffs_noise = zeros(npoints, 8);
            
            for p = 1:npoints
                win = find( (t > (t_freq(p) - windowdur/2)) & ...
                    (t < (t_freq(p) + windowdur/2)));
                taper = hanning(numel(win))';
                
                a_plus_b_minus_ab = ProbeBuffs(k_kept,:) ...
                    + SuppBuffs(k_kept,:) - BothBuffs(k_kept,:);
                resp_trial = a_plus_b_minus_ab(win).* taper; % just curr trial
                
                model_probe = [cos(phiProbe_inst(win)) .* taper;
                    -sin(phiProbe_inst(win)) .* taper];
                model_noise = ...
                    [cos(nearfreqs(1)*phiProbe_inst(win)) .* taper;
                    -sin(nearfreqs(1)*phiProbe_inst(win)) .* taper;
                    cos(nearfreqs(2)*phiProbe_inst(win)) .* taper;
                    -sin(nearfreqs(2)*phiProbe_inst(win)) .* taper;
                    cos(nearfreqs(3)*phiProbe_inst(win)) .* taper;
                    -sin(nearfreqs(3)*phiProbe_inst(win)) .* taper;
                    cos(nearfreqs(4)*phiProbe_inst(win)) .* taper;
                    -sin(nearfreqs(4)*phiProbe_inst(win)) .* taper];
                
                coeffs_resp(p,:) = model_probe' \ resp_trial';
                coeffs_noise(p,:) = model_noise' \ resp_trial';
            end

            % calculate amplitudes
            oae_trials(k_kept,:) = abs(complex(coeffs_resp(:, 1),  coeffs_resp(:, 2)));
            median_oae = median(oae_trials,1);
            sfoae_full = db(median_oae.*VtoSPL); 
            
            noise_trial = zeros(npoints,4);
            for i = 1:2:8
                noise_trial(:,ceil(i/2)) = complex(coeffs_noise(:,i), coeffs_noise(:,i+1));
            end
            noise_trials(k_kept,:) = abs(mean(noise_trial, 2));
            median_noise = median(noise_trials,1); 
            nf_full = db(median_noise.*VtoSPL); 
            
            % Get summary points (weighted by SNR)
            sfoae = zeros(length(centerFreqs),1);
            nf = zeros(length(centerFreqs),1);
            sfoae_w = zeros(length(centerFreqs),1);
            nf_w = zeros(length(centerFreqs),1);
            
                        % weighted average around 9 center frequencies
            for z = 1:length(centerFreqs)
                band = find( testfreq >= bandEdges(z) & testfreq < bandEdges(z+1));
          
                % TO DO: NF from which SNR was calculated included median of 7 points
                % nearest the target frequency.
                SNR = sfoae_full(band) - nf_full(band);
                weight = (10.^(SNR./10)).^2;
                
                sfoae(z, 1) = mean(sfoae_full(band));
                nf(z,1) = mean(nf_full(band));
                
                sfoae_w(z,1) = sum(weight.*sfoae_full(band))/sum(weight);
                nf_w(z,1) = sum(weight.*nf_full(band))/sum(weight);
                
            end
                        
            % median SNR
            SNR_temp = sfoae_w - nf_w;
            
            noisy_trials = 0;
            % artifact check
            if k_kept >= stim.minTrials
                std_oae = std(oae_trials,1);
                for r = 1:k_kept
                    for q = 1:npoints
                        if oae_trials(r,q) > median_oae(1,q) + 3*std_oae(1,q)
                            noisy_trials = noisy_trials+1;
                            break;
                        end
                    end
                end
            end
            
            % if SNR is good enough and we've hit the minimum number of
            % trials, then stop.
            if SNR_temp(1:8) >= SNRcriterion
                if k_kept >= minTrials + noisy_trials
                    doneWithTrials = 1;
                end
            elseif k == maxTrials
                doneWithTrials = 1;
            end
            
            pass = (SNR_temp>=SNRcriterion);
            oae_pass = sfoae_w;
            oae_fail = sfoae_w;
            oae_pass(~pass) = NaN;
            oae_fail(pass) = NaN;
            
            % Plot amplitudes from live analysis
            hold off;
            plot(centerFreqs./1000,oae_pass, 'o', 'linew', 2, 'color', [0 0.4470 0.7410]);
            hold on;
            plot(centerFreqs/1000,oae_fail, 'o', 'linew', 2, 'color', 'k'),
            plot(centerFreqs./1000,nf_w, 'x', 'linew', 2, 'color', [0.6350 0.0780 0.1840]);
            hold off;
            legend('SFOAE', '', 'NOISE', 'location', 'northeast');
            title('SFOAE')
            xlabel('Frequency (Hz)')
            ylabel('Median Amplitude (dB SPL)')
            set(gca, 'XScale', 'log', 'FontSize', 14)
            xlim([0.5, 16]);
            xticks([.5, 1, 2, 4, 8, 16]);
            ylim([-45, 45]);
            yticks((-45:15:45))
            grid on;
            drawnow;
            
            fprintf(1, 'Trials run: %d / Noisy Trials: %d \n', k_kept, noisy_trials);
            
        end
    end
    %% Save full data struct
    data.info = info;
    data.stim = stim;
    data.info.subj = subj;
    data.resp.trialsCollected = k_kept;
    data.resp.ProbeBuffs = ProbeBuffs(1:k_kept,:);
    data.resp.SuppBuffs = SuppBuffs(1:k_kept,:);
    data.resp.BothBuffs = BothBuffs(1:k_kept,:);
    data.resp.testDur_s = toc;
    data.calib.Ph1 = calib_1.calib;
    data.calib.Ph2 = calib_2.calib; 

    save(fname,'data');
    
    %% Close TDT, ER-10X connections etc. and cleanup
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    
catch me
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    rethrow(me);
    
end

