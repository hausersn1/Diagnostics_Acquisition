load('C:\Experiments\TEOAE\calibdata\dat.mat');
y = y(:);
x = x(:);
X = rfft(x);
Y = rfft(y);
f = linspace(0,fs/2,length(X))';
H = Y ./ X;
Htemp = H((f > 8e3) & (f < 10e3));
[~, maxind] = max(abs(Htemp));
ftemp = f((f > 8e3) & (f < 10e3));
fmax = ftemp(maxind);
Hfmax = Htemp(maxind);

H(f > fmax) = Hfmax;

Hinv = 1./H;
b = fir2(1024, f/(fs/2), Hinv);
z = filter(b, 1, y);
plot(x);
hold on;
plot(z);



