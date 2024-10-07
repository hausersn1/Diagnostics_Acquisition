%% CLICK OAE using traditional windowing method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Samantha Hauser, AuD
% Modified from: Hari Bharadwaj, PhD (SNAP Lab)
% Created: November 2021
% Last revision: 16-Sep-2023 (added metadata saving)
%               24-Sep-2024 (correctly applied FPL calib)
%
% References:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up data storage and subject info
clear;
clc;

% Measure-General info
info.measure = 'TEOAE';
info.version = 'v02';

% v01 - original version 
% v02 - started 9/25/24 - w/ calibration 

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
    cd('TEOAE_CALIB')
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
paraDir = 'C:\Experiments\Sam\TEOAE_CALIB\Results\';
respDir = strcat(paraDir,filesep,visit.subj.ID,filesep);
addpath(genpath(paraDir));
if(~exist(respDir,'dir'))
    mkdir(respDir);
end

fname = strcat(respDir, info.measure, '_', ...
    subj.ID, '_', subj.ear, '_', datetag, '.mat');

% load stimulus parameters
click = clickSetDefaults();
driver = 1; % 1 or 2
dr = string(driver);
dr = char(dr);

% Check to make sure there are EarCal files for both drivers
calib1_file = ['C:\Experiments\FPLclick\EARCAL\', subj.ID, '\', 'Calib_Ph', dr , 'ER-10X_', subj.ID, earCode,'Ear_', info.date(1:11), '*.mat']; 

if isempty(dir(calib1_file)) 
    fprintf('------No Ear Calibration Found -- Go run calibEar!!------\n')
    return
else
    % load the calibration files in
    addpath(['C:\Experiments\FPLclick\EARCAL\' subj.ID])
    get_calib_1 = uigetfile(calib1_file, 'Get Calib file');
    calib_1 = load(get_calib_1);
    rmpath(['C:\Experiments\FPLclick\EARCAL\' subj.ID])
end

scaledB = 10; 
   

tic; 
try
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    initializeER10X;
    
    % Initializing TDT and specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    drivername = strcat('Ph',num2str(driver));
    click.drivername = drivername;
    
    % Make click
    vo = clickStimulus(click.BufferSize + click.StimWin);
    
    % generate each filter 
    h1 = flatFPLfilter(calib_1.calib); 
    
    % get the fplh value (FPL1k) at 1000 hZ for each calib file
    FPL1k_1 = interp1(calib_1.calib.freq, db(abs(calib_1.calib.Pfor)), 1000); 
    
    % Filter the stimulus and drop by 30 dB to prevent clipping (will
    % adjust in attenuation (i.e. drop_f1 and drop_f2)
    filt_vo = filter(h1, 1, vo) * db2mag(FPL1k_1) * db2mag(94-FPL1k_1) * db2mag(-scaledB);     
    
    buffdata = zeros(2, numel(filt_vo));
    buffdata(driver, :) = filt_vo; % The other source plays nothing
    click.vo = vo;
    click.filt_vo = filt_vo; 

    delayComp = 1; % 1 always
    filtdelay = 128; % delay for 256 order filter
    odd = 1:2:click.Averages;
    even = 2:2:click.Averages;

    drop = click.Attenuation - scaledB;
    dropOther = 120;
    fprintf('Starting Stimulation...\n')
    
    if driver == 1
        vins = playCapture2(buffdata, card, click.Averages, ...
            click.ThrowAway, drop, dropOther, delayComp);
    else
        vins = playCapture2(buffdata, card, click.Averages, ...
            click.ThrowAway, dropOther, drop, delayComp);
    end
    
    %compute the average
    % window for the stimulus
    vins = vins(:, (click.StimWin+1):(click.StimWin + click.RespDur)); % Remove stimulus by windowing
    % window to handle filter delay        
    raw_nodelay = vins; 
    response = zeros(size(raw_nodelay)); 
    response(:, 1:end-filtdelay) = vins(:, filtdelay+1:end); 
    vins = response; 
    
        
    if click.doFilt
        % High pass at 200 Hz using IIR filter
        [b, a] = butter(4, 200 * 2 * 1e-3/click.SamplingRate, 'high');
        vins = filtfilt(b, a, vins')';
    end
    
    vavg_odd = trimmean(vins(odd, :), 20, 1);
    vavg_even = trimmean(vins(even, :), 20, 1);
    rampdur = 0.2e-3; %seconds
    Fs = click.SamplingRate/2 * 1e3;
    resp.vavg = rampsound((vavg_odd + vavg_even)/2, Fs, rampdur);
    resp.noisefloor = rampsound((vavg_odd - vavg_even)/2, Fs, rampdur);
    
    Vavg = rfft(resp.vavg);
    Vavg_nf = rfft(resp.noisefloor);
    
    % Apply calibrations to convert voltage to pressure
    % For ER-10X, this is approximate
    mic_sens = 50e-3; % mV/Pa. TO DO: change after calibration
    mic_gain = db2mag(gain + 6); % +6 for balanced cable
    P_ref = 20e-6;
    DR_onesided = 1;
    factors = DR_onesided * mic_gain * mic_sens * P_ref;
    resp.output_Pa_per_20uPa = Vavg / factors; % unit: 20 uPa / Vpeak
    resp.noise_Pa_per_20uPa = Vavg_nf / factors;
    
    resp.freq = 1000*linspace(0,click.SamplingRate/2,length(Vavg))';
    
    
    %% Plot data
    %PlotResults_simple;
    
    %% Save Ear Measurements
    data.info = info;
    data.stim = click;
    data.resp = resp; 
    data.info.subj = subj;
    data.resp.allTrials = vins;
    data.resp.testDur_s = toc;
    data.calib = calib_1.calib; 
    %data.resp.raw_nodelay = raw_nodelay; 
    
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