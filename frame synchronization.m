%transmitter
preamble = [1 0 1 0 1 1 0 0 1 1 0 1 1 1 0 1 1 0 1 0 0 1 0 0 1 1 1 0 0 0 ...
    1 0 1 1 1 1 0 0 1 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 1 1 1 1 1 1 0 0];
lenP = length(preamble);
modPreamble = modulator(preamble);
message = randi([0 1], 1, 100); %generate random 100 bits
crc = CRC32generator(message);
info = [preamble message crc];
len = length(info);
signal = modulator(info);

%channel
signal = AWGN(signal,50);

%receiver
recover = correlator(modPreamble, signal, len);
recover = demodulator(recover, len-lenP);
disp(CRC32checker(recover));


function [crc] = CRC32generator(message)
    % CRC-32-IEEE = x^{32}+x^{26}+x^{23}+x^{22}+x^{16}+x^{12}+x^{11}
    %               +x^{10}+x^{8}+x^{7}+x^{5}+x^{4}+x^{2}+x+1
    divider = [1 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 1];
    padMessage = [message zeros(1, length(divider)-1)];
    [~, remainder] = deconv(padMessage, divider); % [quotient, remainder] = padMessage/divider
    remainder = mod(abs(remainder), 2); % convert to binary
    crc = remainder((length(message)+1) : end);
end

function [check] = CRC32checker(info) % check errors in info after transmission
    %info = (message + CRC)
    divider = [1 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 1];
    [~, remainder] = deconv(info, divider);
    remainder = mod(abs(remainder), 2);
    if remainder == 0, check = "error free";
    else, check = "errors"; end
end

function [signal] = correlator(preamble, signal, len)
    lenP = length(preamble);
    r = xcorr(preamble, signal);
    i = max(abs(r))-lenP;
    t = i+1;
    k = i+((len-64)/2);
    signal = signal(t:k);
end

function r = AWGN(s,SNRdB)
    gamma = 10^(SNRdB/10); %SNR to linear scal
    P = sum(abs(s).^2)/length(s);%Actual power in the vector
    N0 = P/gamma; %Find the noise spectral density
    n = sqrt(N0/2) * (randn(size(s)) + 1i * randn(size(s)));%computed noise
    r = s + n; %received signal
end

function [signal] = modulator(info)
    k = 2;
    len = length(info);
    signal = zeros(1, floor(len/2));
    for i = 1:1:len/2
        signal(i) = symbol_mapper((info(floor(((i - 1) * k/ 1) + 1))),...
        (info(floor(((i - 1) * k/1) + 2)))); 
    end
end

function [mapped_value] = symbol_mapper(bit1, bit2)
    mapped_value = (-1)^(bit1) + 1i * (-1)^(bit2);
end

function [received_bits] = demodulator(signal, len)
    k = 2;
    received_bits = zeros(1, len/2 * k);
    % Bit extraction Array
    j = 1;
    for i = 1:1:len/2
        Re = real(signal(i));
        Im = imag(signal(i));
        [received_bits(j), received_bits(j+1)] = demap_vals(Re, Im);
        j = j + k;
    end
end

function [bit0, bit1] = demap_vals(Re, Im)
    if(abs(Re - 1) < 0.01), bit0 = 0;
    else, bit0 = 1; end
    if(abs(Im - 1) < 0.01), bit1 = 0;
    else, bit1 = 1; end
end

% function [signal] = modulation(info) %4QAM
%     N = length(info);
%     T = 1;
%     fc = 2;
%     Fs = 100;
%     info = 2*info-1; % 1 0 -> 1 -1
%     I = []; Q = [];
%     for i = 1:N
%         if mod(i,2)~=0
%             I = [I, info(i)];
%         else
%             Q = [Q, info(i)];
%         end
%     end
%     bitdata = [];
%     for i = 1:N
%         bitdata = [bitdata, info(i)*ones(1,T*Fs)];
%     end
%     Idata = []; Qdata = [];
%     for i = 1:N/2
%         Idata = [Idata, I(i)*ones(1, T*Fs*2)];
%         Qdata = [Qdata, Q(i)*ones(1, T*Fs*2)];
%     end
% %    t = 0:1/Fs:N*T-1/Fs;
% %     figure();
% %     subplot(3,1,1)
% %     plot(t,bitdata);legend('Bitstream')
% %     subplot(3,1,2)
% %     plot(t,Idata);legend('Bitstream')
% %     subplot(3,1,3)
% %     plot(t,Qdata);legend('Bitstream')
%     
%     %carrier
%     bit_t=0:1/Fs:2*T-1/Fs;
%     Icarrier = []; Qcarrier = [];
%     for i = 1:N/2
%         Icarrier=[Icarrier,I(i)*cos(2*pi*fc*bit_t)];
%         Qcarrier=[Qcarrier,Q(i)*cos(2*pi*fc*bit_t+pi/2)];
%     end
%     signal = Icarrier+Qcarrier;
% %     figure();
% %     subplot(3,1,1)
% %     plot(t,Icarrier);legend('I signal')
% %     subplot(3,1,2)
% %     plot(t,Qcarrier);legend('Q signal')
% %     subplot(3,1,3)
% %     plot(t,signal);legend('QPSK signal')
%     %signal = Idata+Qdata;
% end

% function [a] = demodulation(receive, len)    
%     N = len-12;
%     T = 1;
%     fc = 2;
%     Fs = 100;
%     bit_t=0:1/Fs:2*T-1/Fs;
%     for i=1:N/2
%         I_output=receive(1,(i-1)*length(bit_t)+1:i*length(bit_t)).*cos(2*pi*fc*bit_t);
%         if sum(I_output)>0 
%             I_recover(i)=1;
%         else
%             I_recover(i)=-1;
%         end
%          Q_output=receive(1,(i-1)*length(bit_t)+1:i*length(bit_t)).*cos(2*pi*fc*bit_t+ pi/2);
%         if sum(Q_output)>0
%             Q_recover(i)=1;
%         else
%             Q_recover(i)=-1;
%         end
%     end
%     bit_recover=[];
%     for i=1:N
%         if mod(i,2)~=0
%             bit_recover=[bit_recover,I_recover((i-1)/2+1)];
%         else
%             bit_recover=[bit_recover,Q_recover(i/2)];
%         end
%     end
%     
%     a = [];
%     for i=1:N
%         if bit_recover(i) == -1
%            a(i) = 0;
%         end
%     end
%     
%     recover_data=[];
%     for i=1:N
%         recover_data=[recover_data,bit_recover(i)*ones(1,T*Fs)];
%     end
% 
%     I_recover_data=[]; Q_recover_data=[];
%     for i=1:N/2
%         I_recover_data=[I_recover_data,I_recover(i)*ones(1,T*Fs*2)];
%         Q_recover_data=[Q_recover_data,Q_recover(i)*ones(1,T*Fs*2)];
%     end
%     %plot
% %     figure();
% %     t=0:1/Fs:N*T-1/Fs;
% %     subplot(3,1,1)
% %     plot(t,recover_data);legend('Bitstream')
% %     subplot(3,1,2)
% %     plot(t,I_recover_data);legend('I Bitstream')
% %     subplot(3,1,3)
% %     plot(t,Q_recover_data);legend('Q Bitstream')
%     
% end

function [conv_sum] = convert(info)
    Fs = 20;
    T = 0.5;
    len = size(info);
    num_bits = len(2);
    k = 2;
    total_time = num_bits / k * Fs * T;
    b_upsample = zeros(1, total_time);
    trans_rate = floor(Fs * T);
    for i = 1:trans_rate:total_time
        b_upsample(i) = symbol_mapper((info(floor(((i - 1) * k/ trans_rate) + 1))),...
        (info(floor(((i - 1) * k / trans_rate) + 2)))) ; 
    end
    n = -4 * Fs :1: 4 * Fs;
    filter_t = sinc(n / 10);
    conv_sum = conv(filter_t, b_upsample);
    temp = size(conv_sum);
    t = 1:1:temp(2);
    t = t * 1/Fs;
    %conv_sum = conv_sum +0.4*randn(size(t));
    a = fix(conv_sum);
    for i=1:length(a)
        if a(i) <= -1
           a(i) = 1;
        elseif a(i) >=1
            a(i) = 0;
        end
    end
    plot(t, real(conv_sum), 'r')
    hold on
    plot(t, a, 'g')
    hold off
end
%https://blog.csdn.net/qq_44830040/article/details/106433066
%https://www.gaussianwaves.com/2015/06/how-to-generate-awgn-noise-in-matlaboctave-without-using-in-built-awgn-function/
%https://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/
