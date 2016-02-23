% BER simulation of VHT SU 802.11 ac 
% Derived from a skeleton BER script for a wireless link simulation encased in human
% function [ber] = WirelessLinkSimulation(MCS)

clear all; close all; clc; format compact;

% warnings that occur when running remotely
warning('off','MATLAB:xlswrite:AddSheet');
warning('off','MATLAB:xlswrite:NoCOMServer');


%%%%%%%%%%%%%%%%%%%%%%%%%
%        INPUTS         %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify SNR range in dB
SNR_Vec = 0:5:30;

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

% Specify debug mode(0 if running without encoding, else -1)
debug = -1;

% Specify number of iterations for simulation
numIter = 1e1; %1e6


%%%%%%%%%%%%%%%%%%%%%%%%%
%       ENCODING        %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Set modulation and coding scheme
[modType, M, k, R, puncpat hMod, htDemod] = SetMCS(MCS);

% Set encoding method (BCC or LDPC, debug or no debug)
[htConvEnc, htVitDec, htErrorCalc, ...
    N_Pre_Pad, N_Data_Bits, N_Post_Pad] = SetEncoder(encType, debug, k, R, puncpat);

% Run simulation and retrieve BERs
[ber, berHypo] = Simulation(numIter, SNR_Vec, encType, debug, ...
    modType, k, M, hMod, htDemod, ...
    htConvEnc, htVitDec, htErrorCalc, ...
    N_Pre_Pad, N_Data_Bits, N_Post_Pad)


%%%%%%%%%%%%%%%%%%%%%%%%%
%        GRAPHS         %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Graph theoretical and actual BERs
figure
semilogy(SNR_Vec, berHypo, 'r')
hold on
semilogy(SNR_Vec, ber);
hold off
xlabel('SNR (dB)')
legend('Theoretical BER', 'BER');
title(strcat(num2str(M), ' ', modType));
