% BER simulation of VHT SU 802.11 ac 
% Derived from a skeleton BER script for a wireless link simulation encased in human
% function [ber] = WirelessLinkSimulation(MCS)

% This is the main function that runs simulations with all the modulation
% and coding schemes for a given encoding type (or no encoding).

clear all; close all; clc; format compact;

% warnings that occur when running remotely
warning('off','MATLAB:xlswrite:AddSheet');
warning('off','MATLAB:xlswrite:NoCOMServer');
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');


%%%%%%%%%%%%%%%%%%%%%%%%%
%        INPUTS         %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify SNR range in dB
SNR_Vec = 0:20; %0:16

% Specify encoding method (BCC or LDPC)
encType = 'BCC';

% Specify debug mode (0 if running without encoding)
debug = 0;

% Specify number of iterations for simulation
numIter = 10; %1e6


%%%%%%%%%%%%%%%%%%%%%%%%%
%      SIMULATION       %
%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIATING BATTLE MODE.
% waitbars are invalid in parallel pools.
h = waitbar(0, 'Initializing data cannon...');

% BATTLE START!
for MCS = 0:9
    WLinkSim(MCS, SNR_Vec, encType, debug, numIter, h);
end

% Label graph
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
if(debug == 0)
    title('802.11ac: AWGN BER w/o FEC');
    axis([0 20 1e-7 1e0]);
else
    title(strcat('802.11ac:',{' '}, encType));
end

% INITIATING RECOVERY MODE.
% close waitbar
close(h);
