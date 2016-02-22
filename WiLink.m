% BER simulation of VHT SU 802.11 ac 
% Derived from a skeleton BER script for a wireless link simulation encased in human
%function [ber] = WirelessLinkSimulation(MCS)
clear all;close all;clc;format compact; MCS = 3;
% warnings that occur when running remotely
warning('off','MATLAB:xlswrite:AddSheet');
warning('off','MATLAB:xlswrite:NoCOMServer');

SNR_Vec = 0:1:20; % in dB
% The M-ary number, 2 corresponds to binary modulation
msgM = [2 4 16 64 256]; 
k = log2(msgM);   % # of information bits per symbol
puncpat = [1; 1; 1; 0;]; % Rate 2/3 Figure 18-9
puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
puncpat = -1; 
switch MCS
    case 0 
        disp('BPSK Rate 1/2')
        msgM = 2;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.BPSKModulator;
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(1/2);
    case 1 
        disp('QPSK Rate 1/2')
        msgM = 4;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.QPSKModulator;
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(1/2);
    case 2
        disp('QPSK Rate 3/4')
        msgM = 4;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.QPSKModulator;
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(3/4);
    case 3 
        disp('16-QAM Rate 1/2')
        msgM = 16;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(1/2);
    case 4 
        disp('16-QAM Rate 3/4')
        msgM = 16;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(3/4);
    case 5
        disp('64-QAM Rate 2/3')
        msgM = 64;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        puncpat = [1; 1; 1; 0;]; % Rate 2/3 Figure 18-9
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(2/3);
    case 6
        disp('64-QAM Rate 3/4')
        msgM = 64;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(3/4);
    case 7
        disp('64-QAM Rate 5/6')
        msgM = 64;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(5/6);
    case 8
        disp('256-QAM Rate 3/4')
        msgM = 256;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(3/4);
    case 9
        disp('256-QAM Rate 5/6')
        msgM = 256;
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
        EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(5/6);
    otherwise 
        warning('Unexpected MCS.')
end

length_param = 4095;  % 4095 is the max LENGTH parameter. See 18.2.2.2
N_DBPS = k;
N_SYM = ceil((16+8*length_param+6)/N_DBPS); % 18.3.5.4
N_DATA = N_SYM * N_DBPS;

constlen=7;
codegen = [171 133]; 
trellis = poly2trellis(constlen, codegen); % Industry standard 18.3.5.6
htConvEnc = comm.ConvolutionalEncoder(trellis); 
htVitDec = comm.ViterbiDecoder(trellis, 'InputFormat', 'Unquantized');
htVitDec.TracebackDepth = 96; % TODO: find reciever TBD
htErrorCalc = comm.ErrorRate('ReceiveDelay', htVitDec.TracebackDepth);

if (puncpat ~= -1) 
    htConvEnc.PuncturePatternSource = 'Property';
    htVitDec.PuncturePatternSource = 'Property';
    htConvEnc.PuncturePattern = puncpat;
    htVitDec.PuncturePattern = puncpat;
end

%{
hConvEnc = comm.ConvolutionalEncoder(trellis); 
hVitDec = comm.ViterbiDecoder(trellis, 'InputFormat', 'Unquantized');
hVitDec.TracebackDepth = 96; % TODO: find reciever TBD
hErrorCalc = comm.ErrorRate('ReceiveDelay', hVitDec.TracebackDepth);

if (puncpat ~= -1) 
    hConvEnc.PuncturePatternSource = 'Property';
    hVitDec.PuncturePatternSource = 'Property';
    hConvEnc.PuncturePattern = puncpat;
    hVitDec.PuncturePattern = puncpat;
end
%}

% Create a vector to store the BER computed during each iteration
BERVec = zeros(3,length(EbNoEncoderOutput)); % Allocate memory to store results
env_c = length(EbNoEncoderOutput);
frameLength = N_DATA;         % this value must be an integer multiple of 3
numIter = 100 %1e6;

% ScatterPlot! 
%sPlotFig = scatterplot(qammod(0:msgM-1,msgM,0,'gray'),1,0,'k*');

tic;
%waitbars are invalid in parallel pools.
%h = waitbar(0, 'Initializing data cannon...');

% Run the simulation numIter amount of times
% Note that using a parallel pool will not output graphs. 
% Graphs will be generated and can be saved using print
%parfor 
parfor n=1:env_c
  %reset(hErrorCalc)
  %reset(hConvEnc)
  %reset(hVitDec)
  hErrorCalc = htErrorCalc.clone;
  hConvEnc = htConvEnc.clone;
  hVitDec = htVitDec.clone;
  hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)',...
  'SignalPower', 1, 'SamplesPerSymbol', 1);
  hChan.EbNo = EbNoEncoderOutput(n); % Set the channel EbNo value for simulation
  for i = 1:numIter
    % Generate binary frames of size specified by the frameLength variable
    bits = randi([0 1], frameLength, 1);
    % Interleave the bits 
    txdata = bits; %randintrlv(bits,sum(double('Keenesus')));
    % Convolutionally encode the data
    encData = step(hConvEnc, txdata);
    % Modulate the encoded data
    modData = step(hMod, encData);
    % Pass the modulated signal through an AWGN channel
    channelOutput = step(hChan, modData);
    % Pass the real part of the channel complex outputs as the unquantized
    % input to the Viterbi decoder.
    decData = step(hVitDec, real(channelOutput));
    % Deinterleave the bits
    data = decData; %randdeintrlv(decData,sum(double('Keenesus')));
    % Compute and accumulate errors
    BERVec(:,n) = step(hErrorCalc, bits, data);
  end
  
  %waitbar(n/env_c,h,sprintf('Evaluating %d',EbNoEncoderOutput(n)))
end % End iteration
toc;
%close(h);

ber = BERVec(1,:)


%{
% Compute the theoretical BER for this scenario
SNR_Vec = EbNoEncoderInput;
berHypo = berawgn(SNR_Vec, 'qam', msgM, 'nondiff');
figure
semilogy(SNR_Vec,berHypo,'r')
hold on
semilogy(SNR_Vec,BERVec(1,:));
hold off
legend('Theoretical BER', 'BER');
title(strcat(num2str(msgM), ' QAM'));


filename = sprintf('snr%d_ber%dQAM%dfw%dbk%dk%.3fff%.2fis.csv',snr,QAM,nfwdweights,nfbkweights,Kappa,ff,is);
xlswrite(filename,berVec,'berVec');
filename = sprintf('snr%d_br%dQAM%dfw%dbk%dk%.3fff%.2fis.csv',snr,QAM,nfwdweights,nfbkweights,Kappa,ff,is);
xlswrite(filename,brVec,'brVec');
filename = sprintf('snr%d_ew%dQAM%dfw%dbk%dk%.3fff%.2fis.csv',snr,QAM,nfwdweights,nfbkweights,Kappa,ff,is);
xlswrite(filename,ewVec,'ewVec');
filename = sprintf('snr%d_eqchan%dQAM%dfw%dbk%dk%.3fff%.2fis.csv',snr,QAM,nfwdweights,nfbkweights,Kappa,ff,is);
xlswrite(filename,eqchan,'eqchan');
filename = sprintf('snr%d_eqweights%dQAM%dfw%dbk%dk%.3fff%.2fis.csv',snr,QAM,nfwdweights,nfbkweights,Kappa,ff,is);
xlswrite(filename,eqweights,'eqweights');
%}
%endtime = datetime('now');
