% BER simulation of VHT SU 802.11 ac 
% Derived from a skeleton BER script for a wireless link simulation encased in human
%function [ber] = WirelessLinkSimulation(MCS)
clear all;close all;clc;format compact; 
% warnings that occur when running remotely
warning('off','MATLAB:xlswrite:AddSheet');
warning('off','MATLAB:xlswrite:NoCOMServer');

% Inputs
MCS = 3;
type = 'BCC';
numIter = 1e3 %1e6;
SNR_Vec = 0:1:10; % in dB

% Reference Materials
% The M-ary number, 2 corresponds to binary modulation
msgM = [2 4 16 64 256]; 
k = log2(msgM);   % # of information bits per symbol
puncpat = [1; 1; 1; 0;]; % Rate 2/3 Figure 18-9
puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
puncpat = -1; % Rate 1/2 Default Rate; No puncture 

% Choosing which Modulation and Coding Scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Functionalize this shit. Add LDPC switch cases
% if (strcomp(type,'LDPC'))
% else if (strcomp(type'BCC')) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch MCS
    case 0 
        disp('BPSK Rate 1/2')
        R = 1/2;
        msgM = 2;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'psk', msgM, 'nondiff');
        hMod = comm.BPSKModulator;
        hDeMod = comm.BPSKDemodulator;
    case 1 
        disp('QPSK Rate 1/2')
        R = 1/2;
        msgM = 4;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'psk', msgM, 'nondiff'); 
        hMod = comm.QPSKModulator;
        hDeMod = comm.QPSKDemodulator;
    case 2
        disp('QPSK Rate 3/4')
        R = 3/4;
        msgM = 4;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'psk', msgM, 'nondiff');
        hMod = comm.QPSKModulator;
        hDeMod = comm.QPSKDemodulator;
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
    case 3 
        disp('16-QAM Rate 1/2')
        R = 1/2;
        msgM = 16;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec + 0*log10(k), 'qam', msgM, 'nondiff');
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        hDeMod = comm.RectangularQAMDemodulator(msgM);
    case 4 
        disp('16-QAM Rate 3/4')
        R = 3/4;
        msgM = 16;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'qam', msgM, 'nondiff');
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        hDeMod = comm.RectangularQAMDemodulator(msgM);
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
    case 5
        disp('64-QAM Rate 2/3')
        R = 2/3;
        msgM = 64;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'qam', msgM, 'nondiff');
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        hDeMod = comm.RectangularQAMDemodulator(msgM);
        puncpat = [1; 1; 1; 0;]; % Rate 2/3 Figure 18-9
    case 6
        disp('64-QAM Rate 3/4')
        R = 3/4;
        msgM = 64;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'qam', msgM, 'nondiff');
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        hDeMod = comm.RectangularQAMDemodulator(msgM);
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
    case 7
        disp('64-QAM Rate 5/6')
        R = 5/6;
        msgM = 64;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'qam', msgM, 'nondiff');
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        hDeMod = comm.RectangularQAMDemodulator(msgM);
        puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
    case 8
        disp('256-QAM Rate 3/4')
        R = 3/4;
        msgM = 256;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'qam', msgM, 'nondiff');
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        hDeMod = comm.RectangularQAMDemodulator(msgM);
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
    case 9
        disp('256-QAM Rate 5/6')
        R = 5/6;
        msgM = 256;
        k = log2(msgM);   % # of information bits per symbol
        berHypo = berawgn(SNR_Vec - 10*log10(k), 'qam', msgM, 'nondiff');
        hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
        hDeMod = comm.RectangularQAMDemodulator(msgM);
        puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
    otherwise 
        warning('Unexpected MCS.')
end

% Math Setup for # of bits
length_param = 4095;  % 4095 is the max LENGTH parameter. See 18.2.2.2
N_DBPS = k;
N_Scrambler_Init_Bits = 7;
N_Reserved_Service_Bits = 9;
N_Tail_bits = 6;
N_SYM = ceil((N_Scrambler_Init_Bits+N_Reserved_Service_Bits+...
    8*length_param+N_Tail_bits)/N_DBPS); % 18.3.5.4
N_DATA = N_SYM * N_DBPS;
N_PAD = N_DATA - (N_Scrambler_Init_Bits+N_Reserved_Service_Bits+...
    8 * length_param + N_Tail_bits);
N_Punc_Pad = mod(length(puncpat) - N_PAD - mod(N_DATA,length(puncpat)), ...
    length(puncpat)); % Padding for computation
N_Punc_Pad = N_Punc_Pad + mod(N_DATA + N_Punc_Pad + N_PAD + N_Tail_bits, lcm(k,length(puncpat)));

% Convolutional Encoding Setup
constlen=7;
codegen = [171 133]; 
trellis = poly2trellis(constlen, codegen); % Industry standard 18.3.5.6
htConvEnc = comm.ConvolutionalEncoder(trellis); 
htVitDec = comm.ViterbiDecoder(trellis, 'InputFormat', 'Hard'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
htVitDec.TracebackDepth = 96; % TODO: find reciever TBD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
htErrorCalc = comm.ErrorRate('ReceiveDelay', htVitDec.TracebackDepth);
htErrorCalc = comm.ErrorRate; % for debug no encoding

if (puncpat ~= -1) 
    htConvEnc.PuncturePatternSource = 'Property';
    htVitDec.PuncturePatternSource = 'Property';
    htConvEnc.PuncturePattern = puncpat;
    htVitDec.PuncturePattern = puncpat;
end

% Create a vector to store the BER computed during each iteration
BERVec = zeros(3,length(SNR_Vec)); % Allocate memory to store results
env_c = length(SNR_Vec);
frameLength = N_DATA;

tic;
%waitbars are invalid in parallel pools.
%h = waitbar(0, 'Initializing data cannon...');

% Run the simulation numIter amount of times
% Note that using a parallel pool will not output graphs. 
% Graphs will be generated and can be saved using print
for n=1:env_c
  %reset(hErrorCalc)
  %reset(hConvEnc)
  %reset(hVitDec)
  hErrorCalc = htErrorCalc.clone;
  hConvEnc = htConvEnc.clone;
  hVitDec = htVitDec.clone;
  hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)',...
  'SignalPower', 1, 'SamplesPerSymbol', 1, 'BitsPerSymbol', k);
  hChan.EbNo = SNR_Vec(n); % Set the channel EbNo value for simulation
  %{
  hChan = comm.AWGNChannel(...
            'NoiseMethod','Signal to noise ratio (SNR)','SNR', ...
            SNR_Vec(n) - 10*log10(k/msgM) - 10*log10(k));% + 10*log10(R));
  % http://www.mathworks.com/examples/matlab-communications/mw/comm-ex70334664-punctured-convolutional-coding
  %}
  for i = 1:numIter
    % Generate binary frames of size specified by the frameLength variable
    bits = [zeros(N_Scrambler_Init_Bits+N_Reserved_Service_Bits,1); ...
        randi([0 1], frameLength-N_Scrambler_Init_Bits-N_Reserved_Service_Bits, 1); ...
        zeros(N_Tail_bits+N_PAD+N_Punc_Pad,1)];
    % Interleave the bits   % Not interleaving because parity bit math mess
    txdata = bits; %randintrlv(bits,sum(double('Keenesus'))); 
    % Convolutionally encode the data
    encData = txdata;%step(hConvEnc, txdata);
    % Convert bits to symbols
    txsyms = bi2de(reshape(encData,k,numel(encData)/k).','left-msb');
    % Modulate the encoded data
    modData = qammod(txsyms,msgM,0,'gray');%step(hMod, txsyms);
    % Pass the modulated signal through an AWGN channel
    channelOutput = step(hChan, modData);
    % Demodulate the signal 
    rxsyms = qamdemod(channelOutput,msgM,0,'gray');%step(hDeMod, channelOutput);
    % Convert symbols to bits
    rxbits = de2bi(rxsyms,k,'left-msb').'; % Map Symbols to Bits
    % Pass the demodulated channel outputs as input to the Viterbi decoder
    decData = rxbits(:);%step(hVitDec, rxbits(:));
    % Deinterleave the bits % Not interleaving because parity bit math mess
    data = decData; %randdeintrlv(decData,sum(double('Keenesus')));
    % Compute and accumulate errors
    BERVec(:,n) = step(hErrorCalc, bits, data);
  end
  
  %waitbar(n/env_c,h,sprintf('Evaluating %d',EbNoEncoderOutput(n)))
end % End iteration
toc;
%close(h);

ber = BERVec(1,:)



% Compute the theoretical BER for this scenario
%berHypo = berawgn(SNR_Vec, 'qam', msgM, 'nondiff');
%berHypo = berawgn(SNR_Vec, 'psk', msgM, 'nondiff');
figure
semilogy(SNR_Vec,berHypo,'r')
hold on
semilogy(SNR_Vec,BERVec(1,:));
hold off
legend('Theoretical BER', 'BER');
title(strcat(num2str(msgM), ' QAM'));

%{
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
%filename = sprintf('snr%d:%d_MCS%d_Type%s.csv',snr(1),snr(end),MCS,type);
%xlswrite(filename,ber,'berVec');
%endtime = datetime('now');
