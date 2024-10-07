function click = clickSetDefaults(att)

%set default parameters

if ~exist('att', 'var')
    click.Attenuation = 15 +6; % changed from 60 on 10/4/24; % assume max is like 115 is bw=7500, spectrum level = 40
else
    click.Attenuation = att;
end
click.Vref  = 1;
click.BufferSize = 2048;
click.RespDur = 1024;
click.SamplingRate = 48.828125; %kHz
click.Averages = 2048;
click.ThrowAway = 8;
click.doFilt = 1;
click.StimWin = 128;
click.device = 'ER10X';
