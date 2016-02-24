% BER simulation of VHT SU 802.11 ac 
% Derived from a skeleton BER script for a wireless link simulation encased in human
% function [ber] = WirelessLinkSimulation(MCS)

clear all; close all; clc; format compact;

% warnings that occur when running remotely
warning('off','MATLAB:xlswrite:AddSheet');
warning('off','MATLAB:xlswrite:NoCOMServer');
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');


%%%%%%%%%%%%%%%%%%%%%%%%%
%        INPUTS         %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify SNR range in dB
SNR_Vec = 0:16; %0:16

% Specify modulation and coding scheme (ranges from 0-9)
MCS = 9;
%{
0: BPSK Rate 1/2
1: QPSK Rate 1/2
2: QPSK Rate 3/4
3: 16-QAM Rate 1/2
4: 16-QAM Rate 3/4
5: 64-QAM Rate 2/3
6: 64-QAM Rate 3/4
7: 64-QAM Rate 5/6
8: 256-QAM Rate 3/4
9: 256-QAM Rate 5/4
%}

% Specify encoding method (BCC or LDPC)
encType = 'BCC';

% Specify debug mode (0 if running without encoding)
debug = 0;

% Specify number of iterations for simulation
numIter = 1; %1e6

%fprintf('MCS:\tType:\tIter:\t\n%d\t%s\t%d\n', MCS, encType, numIter);
%%%%%%%%%%%%%%%%%%%%%%%%%
%       ENCODING        %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Set modulation and coding scheme
[modType, M, k, R, k_TCB, puncpat, hMod, htDemod] = SetMCS(MCS);

% Set encoding method (BCC or LDPC, debug or no debug)
[htEnc, htDec, htErrorCalc, ...
    N_Pre_Pad, N_Data_Bits, N_Post_Pad] = ...
    SetEncoder(encType, debug, k, R, k_TCB, puncpat);

% Run simulation and retrieve BERs
[ber, berHypo] = Simulation(numIter, SNR_Vec, encType, debug, ...
    modType, k, R, M, hMod, htDemod, ...
    htEnc, htDec, htErrorCalc, ...
    N_Pre_Pad, N_Data_Bits, N_Post_Pad)

check = ((berHypo - ber(1,:)) > 0);

%%%%%%%%%%%%%%%%%%%%%%%%%
%        GRAPHS         %
%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%d\t%s\t%d\t%d\t%d\n', MCS, encType, numIter, ber(3,1),  toc);
% Graph theoretical and actual BERs
figure
    if (debug ~= 0) % for encoding
        berHypo = berawgn(SNR_Vec - 10*log10(k*R), modType, M, 'nondiff');
    else 
        berHypo = berawgn(SNR_Vec - 10*log10(k), modType, M, 'nondiff');
    end
semilogy(SNR_Vec,berHypo,'r')
hold on
semilogy(SNR_Vec,ber(1,:));
hold off
xlabel('SNR (dB)')
legend('Theoretical BER', 'BER');
title(strcat(num2str(M), ' ', modType));
