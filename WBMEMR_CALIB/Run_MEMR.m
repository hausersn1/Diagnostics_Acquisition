%% WB-MEMR for Humans (ER-10X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Hari Bharadwaj, PhD (SNAP Lab)
% Modified by: Samantha Hauser, AuD
% Created:
% Last revision: 16-Sep-2023 (added artifact checking)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

% Measure-General info
info.version = 'v02';
info.measure = 'WBMEMR';

% v01 - original version
% v02 - with proper calibration, started 9/25/24

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
    cd('WBMEMR_CALIB')
end

subj = visit.subj;
info.room = visit.room;
info.univ = visit.univ;
info.researcher = visit.researcher;

% Get ear info
subj.ear = questdlg('Which ear?', 'Ear', 'L', 'R', 'R');
earCode = subj.ear;

% Get date/time
datetag = datestr(clock);
info.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '-';

% Make directory to save results
paraDir = 'C:\Experiments\Sam\WBMEMR_CALIB\Results\';
respDir = strcat(paraDir,filesep,visit.subj.ID,filesep);
addpath(genpath(paraDir));
if(~exist(respDir,'dir'))
    mkdir(respDir);
end

% Check to make sure there are EarCal files for both drivers
calib1_file = ['C:\Experiments\FPLclick\EARCAL\', subj.ID, '\', 'Calib_Ph1ER-10X_', subj.ID, earCode,'Ear_', info.date(1:11), '*.mat'];
calib2_file = ['C:\Experiments\FPLclick\EARCAL\', subj.ID, '\', 'Calib_Ph2ER-10X_', subj.ID, earCode,'Ear_', info.date(1:11), '*.mat'];

if (isempty(dir(calib1_file)) || isempty(dir(calib2_file)))
    fprintf('------No Ear Calibration Found -- Go run calibEar!!------\n')
    return
else
    % load the calibration files in
    get_calib_1 = uigetfile(calib1_file, 'Get Driver 1 calib file');
    get_calib_2 = uigetfile(calib2_file, 'Get Driver 2 calib file');
    calib_1 = load(get_calib_1);
    calib_2 = load(get_calib_2);
    
    % generate each filter
    h1 = flatFPLfilter(calib_1.calib);
    h2 = flatFPLfilter(calib_2.calib);

    % get the fplh value (FPL1k) at 1000 hZ for each calib file
    FPL1k_1 = interp1(calib_1.calib.freq, db(abs(calib_1.calib.Pfor)), 1000);
    FPL1k_2 = interp1(calib_2.calib.freq, db(abs(calib_2.calib.Pfor)), 1000);

end

%% Plays clicks and noise, records response
% Initializing TDT
fig_num=99;
GB_ch=1;
FS_tag = 3;
Fs = 48828.125;
[f1RZ,RZ,~]=load_play_circuit(FS_tag,fig_num,GB_ch);

% Initialize ER-10X  (Also needed for ER-10C for calibrator)
initializeER10X;


%% DO FULL BAND FIRST
info.measure = 'WBMEMR_Full';
fname = strcat(respDir, info.measure, '_', ...
    subj.ID, '_', subj.ear, '_', datetag, '.mat');

stim = makeMEMRstim_500to8500Hz_shortversion;

pause(3);

tic;
%Set the delay of the sound
invoke(RZ, 'SetTagVal', 'onsetdel',0); % onset delay is in ms
playrecTrigger = 1;
filtdelay = 128; 

%% Set attenuation and play
resplength = numel(stim.t);
resp = zeros(stim.nLevels, stim.Averages, stim.nreps, resplength);
for L = 1:stim.nLevels
    invoke(RZ, 'SetTagVal', 'attA', stim.clickatt - 30);
    invoke(RZ, 'SetTagVal', 'attB', stim.noiseatt(L)-30);
    invoke(RZ, 'SetTagVal', 'nsamps', resplength);

    for n = 1: (stim.Averages + stim.ThrowAway)

        % buffdata = [stim.click; squeeze(stim.noise(L, n, :))'];

        buffdataL_raw = stim.click;
        buffdataR_raw = squeeze(stim.noise(L, n, :))';

        % Filter the stimulus and drop by 30 dB to prevent clipping (will
        % adjust in attenuation (i.e. drop_f1 and drop_f2)
        buffdataL = filter(h1, 1, buffdataL_raw) * db2mag(FPL1k_1) * db2mag(-30);
        buffdataR = filter(h2, 1, buffdataR_raw)* db2mag(FPL1k_2) * db2mag(-30);

        % Check for clipping and load to buffer
        if(any(abs(buffdataL(:)) > 1) || any(abs(buffdataR(:)) > 1))
            error('What did you do!? Sound is clipping!! Cannot Continue!!\n');
        end
        %Load the 2ch variable data into the RZ6:
        %invoke(RZ, 'WriteTagVEX', 'datain', 0, 'I16', (buffdata*2^15));
        invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', buffdataL);
        invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', buffdataR);
        pause(1.5);
        for k = 1:stim.nreps
            %Start playing from the buffer:
            invoke(RZ, 'SoftTrg', playrecTrigger);
            currindex = invoke(RZ, 'GetTagVal', 'indexin');
            while(currindex < resplength)
                currindex=invoke(RZ, 'GetTagVal', 'indexin');
            end

            vin = invoke(RZ, 'ReadTagVex', 'dataout', 0, resplength,...
                'F32','F64',1);
            %Accumluate the time waveform - no artifact rejection
            if (n > stim.ThrowAway)
%                 response = zeros(size(vin));
%                 response(1:end-filtdelay) = vin(filtdelay+1:end);
%                 resp(L, n-stim.ThrowAway, k, :) = response;
                resp(L, n-stim.ThrowAway, k, :) = vin; 
            end

            % Get ready for next trial
            invoke(RZ, 'SoftTrg', 8); % Reset OAE buffer

            fprintf(1, 'Done with Level #%d, Trial # %d \n', L, n);
        end

    end
    pause(2);
end


%% Info for conversion.. no averaging or conversion done online

mic_sens = 0.05; % V / Pa-RMS
mic_gain = db2mag(36);
P_ref = 20e-6; % Pa-RMS

DR_onesided = 1;

stim.mat2Pa = 1 / (DR_onesided * mic_gain * mic_sens * P_ref);


%% Save results

data.info = info;
data.stim = stim;
data.info.subj = subj;
data.resp.AllBuffs = resp;
data.resp.testDur_s = toc;
data.calib.Ph1 = calib_1.calib;
data.calib.Ph2 = calib_2.calib;

save(fname,'data');

%% DO HIGH BAND NEXT
clear stim; clear data;

info.measure = 'WBMEMR_HP';
fname = strcat(respDir, info.measure, '_', ...
    subj.ID, '_', subj.ear, '_', datetag, '.mat');

stim = makeMEMRstim_3kto11k_shortversion;
pause(3);

%Set the delay of the sound
tic;
invoke(RZ, 'SetTagVal', 'onsetdel',0); % onset delay is in ms
playrecTrigger = 1;


%% Set attenuation and play
resplength = numel(stim.t);
resp = zeros(stim.nLevels, stim.Averages, stim.nreps, resplength);
for L = 1:stim.nLevels
    invoke(RZ, 'SetTagVal', 'attA', stim.clickatt-30);
    invoke(RZ, 'SetTagVal', 'attB', stim.noiseatt(L)-30);
    invoke(RZ, 'SetTagVal', 'nsamps', resplength);

    for n = 1: (stim.Averages + stim.ThrowAway)

        % buffdata = [stim.click; squeeze(stim.noise(L, n, :))'];

        buffdataL_raw = stim.click;
        buffdataR_raw = squeeze(stim.noise(L, n, :))';

        % Filter the stimulus and drop by 30 dB to prevent clipping (will
        % adjust in attenuation (i.e. drop_f1 and drop_f2)
        buffdataL = filter(h1, 1, buffdataL_raw) * db2mag(FPL1k_1) * db2mag(-30);
        buffdataR = filter(h2, 1, buffdataR_raw)* db2mag(FPL1k_2) * db2mag(-30);

        % Check for clipping and load to buffer
        if(any(abs(buffdataL(:)) > 1) || any(abs(buffdataR(:)) > 1))
            error('What did you do!? Sound is clipping!! Cannot Continue!!\n');
        end
        %Load the 2ch variable data into the RZ6:
        %invoke(RZ, 'WriteTagVEX', 'datain', 0, 'I16', (buffdata*2^15));
        invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', buffdataL);
        invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', buffdataR);
        pause(1.5);
        for k = 1:stim.nreps
            %Start playing from the buffer:
            invoke(RZ, 'SoftTrg', playrecTrigger);
            currindex = invoke(RZ, 'GetTagVal', 'indexin');
            while(currindex < resplength)
                currindex=invoke(RZ, 'GetTagVal', 'indexin');
            end

            vin = invoke(RZ, 'ReadTagVex', 'dataout', 0, resplength,...
                'F32','F64',1);
            %Accumluate the time waveform - no artifact rejection
            if (n > stim.ThrowAway)
%                 response = zeros(size(vin));
%                 response(1:end-filtdelay) = vin(filtdelay+1:end);
%                 resp(L, n-stim.ThrowAway, k, :) = response;
                resp(L, n-stim.ThrowAway, k, :) = vin;
            end

            % Get ready for next trial
            invoke(RZ, 'SoftTrg', 8); % Reset OAE buffer

            fprintf(1, 'Done with Level #%d, Trial # %d \n', L, n);
        end

    end
    pause(2);
end


%% Info for conversion.. no averaging or conversion done online

mic_sens = 0.05; % V / Pa-RMS
mic_gain = db2mag(36);
P_ref = 20e-6; % Pa-RMS

DR_onesided = 1;

stim.mat2Pa = 1 / (DR_onesided * mic_gain * mic_sens * P_ref);


%% Save results

data.info = info;
data.stim = stim;
data.info.subj = subj;
data.resp.AllBuffs = resp;
data.resp.testDur_s = toc;
data.calib.Ph1 = calib_1.calib;
data.calib.Ph2 = calib_2.calib;

save(fname,'data');
%% Close and clean up
close_play_circuit(f1RZ, RZ);
closeER10X;
