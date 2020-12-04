%% Pre-defined 
syms f f0 f1 bups t n Xbase
fc = 0; % carrier frequency (Hz)
fs = 20; % sampling rate (Hz)
M = 1; % group size

%% Request user input
prompt = 'Please write down the inputs as an array [Transimission rate(bps), "b[k]", "f0", "f1", "f2", "f3"]: ';
% [2, "00011011","sinc(2*t)", "-sinc(2*t)", "-sinc(4*t)", "sinc(4*t)"]
in = input(prompt);
bps = str2double(in(1)); % data rate (bps)
bits = num2str(in(2))-'0'; % k groups of M bits (textscan(num2str(in(4)), '%1s'))'
f0 = str2sym(in(3)); % f0(t)
f0 = rcosdesign(0.25, 20, 10);
f1 = str2sym(in(4)); % f1(t)
%f2 = str2sym(in(5));
%f3 = str2sym(in(6));
% str2num(num2str(reshape(bits,2,[])', '%d'))   

%% Analyze the inputs
k = length(bits)/M;% number of groups
bits = reshape(bits,M,[])'; % reshape to k*M matrix
T = M/bps; % time shifting parameter (sec)
% k*T = overall transmission time (sec)
N = floor(fs * T); % period of xbase[n]
% k*N = overall entries of Xbase[n]

%% Sampling input signals
n = -4*fs : 1 : 4*fs; % take one main lobe and four side-lobes of sinc function (4 to 10 side-lobes is good)
n(4*fs + 1) = 0.00000001; % "subs" function changes sinc(t) to sin(2pit)/(2pit), and 0 cannot be divided.
f = double(subs(f0, t, n/N)); % f0[n/N]
f = f0;
%% Upsampling
bups = zeros(1, k*N);
for i = 1:N:k*N
    bups(i) = symbol_mapper(bits((floor((i+N-1)/N)),:),M);
end

%% Convolution sum
Xbase = conv(bups,f); %convolution sum of b_upsampling[n] and f0[n]
Xbase2 = conv(Xbase,f); 
%% Plot
%t = 0:1/fs:(k*N-1)/fs;
temp = length(Xbase2);
t = 0: 1/fs: (temp-1)/fs;
subplot(2,1,1)
plot(t, real(Xbase2))
xlabel('t');ylabel('I'); title('Real part of Xbase(t)');
grid on;  
xticks(0:N/fs:k*N/fs);
yticks(-1:1:1);

subplot(2,1,2)
plot(t, imag(Xbase2))
xlabel('t');ylabel('Q'); title('Imaginary part of Xbase(t)');
grid on;  
xticks(0:N/fs:k*N/fs);
yticks(-1:1:1);

%% Symbol mapper function
function [f_sym] = symbol_mapper(x,M)
f_sym = (-1)^x(1)+1i*(-1)^x(M);
%f_sym = 0;
    for idx = 1:1:(M/2-1)%M=4; idx= 1,2,3
        f_sym = f_sym + 2^(idx)*(-1)^(x(M/2-idx))+ 1i * 2^(idx)*(-1)^(x(M-idx));
    end
end
