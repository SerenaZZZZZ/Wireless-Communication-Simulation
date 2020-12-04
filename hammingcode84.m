%(8,4)hamming code
p = 4; % # of parity bits per block
n = 8; % # of total bits per block = 2^(p-1)
k = 4; % # of message bits per block = n-p

%% demonstartion part
M = [0 1 1 1]; % message array
C = generator(M); % generate 7 bit wordcode
C= [0     0     0      0     1     1     1     1]; % wordcode with one error
P = check(C); % check the error position
correct(C,P); % correct the error and print the right wordcode and message array

%% scheme:     
%              column1    column2    column3    column4
%       row1   bit0(p0)   bit1(p1)   bit2(p2)   bit3(a)
%       row2   bit4(p3)   bit5(b)    bit6(c)    bit7(d)

%% generator: generate 8 bits hamming code before transmission
function [ C ] = generator(M)
    %M = [a,b,c,d];
    % generator matrix G
    G = [1     1     1     0     0     0     0; %   
         1     0     0     1     1     0     0; %   
         0     1     0     1     0     1     0; %  
         1     1     0     1     0     0     1];%   
    C = mod(M*G, 2); % [a+b+d, a+c+d, a, b+c+d, b, c, d]
    bit0 = mod(C*[1 1 1 1 1 1 1]', 2);
    C = [bit0, C];% [parity0, a+b+d, a+c+d, a, b+c+d, b, c, d]
end

%% check wrong bit position after transmission
function [ P ] = check(C)
    % check matrix H
    H = [0     0     0     0     1     1     1     1   % parity3: check bit 4567 (row 2)
         0     0     1     1     0     0     1     1   % parity2: check bit 2367 (column 3 4)
         0     1     0     1     0     1     0     1   % parity1: check bit 1357 (column 2 4)
         1     1     1     1     1     1     1     1]; % parity0: check all
    S = mod(H*C', 2); 
    p0 = S(4);
    p123 = S(1:3);
    if p0==0 % parity0 = 0
        if p123==0
            P = 8; % no error
        else
            P = -1;% two error
        end
    else % parity0 = 1
        if p123==0
            P = 0;% the wrong position is bit 0
        else
            P = bin2dec(num2str(p123'));% the wrong position
        end
    end
end

%% ccorrected the wrong bit
function [ C ] = correct(C,P)
    M = [C(3+1),C(5+1),C(6+1),C(7+1)];
    if P == -1
        fprintf('Two errors with unknowing positions\nReceived codeword:\n'); disp(C)
    elseif P == 8
        fprintf('Nothing wrong with this block:\n'); disp(C)
        fprintf('Code should be: '); disp(M)
    else
        fprintf('Bit %i is incorrect.\nReceived codeword:\n',P); disp(C)
        fprintf('Corrected codeword:\n');
        C(P+1) = ~C(P+1); disp(C)
        fprintf('Code should be: '); disp(M)
    end
end

