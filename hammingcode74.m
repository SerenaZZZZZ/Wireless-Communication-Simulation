%(7,4)hamming code
p = 3; % # of parity bits per block
n = 7; % # of total bits per block n = 2^p-1
k = 4; % # of message bits per block k = n-p

%% demonstarte part
M = [0 1 1 1]; % message array
C = generator(M); % generate 7 bit wordcode
%C = [0     0     0     1     1     1     1]; % wordcode with one error
P = check(C); % check the error position
correct(C,P); % correct the error and print the right wordcode and message array

%% scheme:     
%              column1    column2    column3    column4
%       row1              bit1(p1)   bit2(p2)   bit3(a)
%       row2   bit4(p3)   bit5(b)    bit6(c)    bit7(d)

%%
function [ C ] = generator(M)
    %M = [a,b,c,d];

    % generator matrix G
    G = [1     1     1     0     0     0     0; % a
         1     0     0     1     1     0     0; % b
         0     1     0     1     0     1     0; % c
         1     1     0     1     0     0     1];% d
    C = mod(M*G, 2); % [a+b+d, a+c+d, a, b+c+d, b, c, d]
end

%%
function [ P ] = check(C)
    % check matrix H
    H = [0     0     0     1     1     1     1   % parity3: check bit 4567 (row 2)
         0     1     1     0     0     1     1   % parity2: check bit 2367 (column 3 4)
         1     0     1     0     1     0     1]; % parity1: check bit 1357 (column 2 4)
    % so that mod(H*G',2)=0
    S = mod(H*C', 2); 
    P = bin2dec(num2str(S'));% the wrong position
end

%%
function [ C ] = correct(C,P)
    if (P ~= 0)
        fprintf('Bit %i is incorrect.\n Received codeword:\n',P)
        disp(C)
        C(P) = ~C(P);
        fprintf('Corrected codeword:\n')
        disp(C)
    else
        fprintf('Nothing wrong with this block:\n')
        disp(C)  
    end
    fprintf('Code should be [%i %i %i %i].\n',C(3),C(5),C(6),C(7))
end

