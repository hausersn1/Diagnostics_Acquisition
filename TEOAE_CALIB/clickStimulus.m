function y = clickStimulus(N, noiseon, fs, ipsicontrabi)
% function y = clickStimulus (N, noiseon, fs, ipsicontrabi)
% Generates a 1 sample click with a total buffer duration of N samples

if ~exist('noiseon', 'var')
    noiseon = 0;
end

if ~exist('ipsicontrabi', 'var')
    ipsicontrabi = 1;
end

if ~exist('fs', 'var')
    fs = 48828.125;
end
if noiseon %Add 50 ms noise
    bw = 7500;
    fc = 4250;
    tmax = 0.05;
    ramp = 0.005;
    noise = makeEqExNoiseFFT(bw, fc, tmax, fs, ramp, 0);
    tail = zeros(1024, 1);
    noiseatt = 0;
    noise = [noise(:)*db2mag(-noiseatt); tail];
    
    noise2 = makeEqExNoiseFFT(bw, fc, tmax, fs, ramp, 0);
    noise2 = [noise2(:)*db2mag(-noiseatt); tail];
end
nsampsclick = 4; % 82 microsecond
initbuff =  2;
y = zeros(1, N);
y(initbuff + (1:nsampsclick)) = 0.95;

if noiseon
    switch ipsicontrabi
        case 1 % Ipsi
            y = [noise; y(:)];
            y(:, 2) = zeros(size(y));
        case 2 % Contra
            y = [zeros(size(noise)); y(:)];
            y(:, 2) = [noise; zeros(N, 1)];
        case 3
            y = [noise; y(:)];
            y(:, 2) = [noise2; zeros(N, 1)];
        otherwise
            error('You messed up! Ipsi should be 1, contra 2 and bi 3!\n');
    end
else
    y = y(:); % Just in case
end