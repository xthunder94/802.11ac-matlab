% BER simulation of VHT SU 802.11 ac 
% Derived from a skeleton BER script for a wireless link simulation encased in human
%function [ber] = WirelessLinkSimulation(MCS)
clear all;close all;clc;format compact; 
% warnings that occur when running remotely
warning('off','MATLAB:xlswrite:AddSheet');
warning('off','MATLAB:xlswrite:NoCOMServer');
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

%% Inputs
MCS = 6 % 0:9;
type = 'BCC'; % ['BCC' 'LDPC'];
numIter = 1e2 %1e6; %TODO: Include suggested values or std bits simulated
SNR_Vec = 0:5:30; % in dB
debug = 1; % If 0, running without encoding

%% Choosing which Modulation and Coding Scheme 
switch MCS
    case 0 
        disp('BPSK Rate 1/2')
        modType = 'PSK';
        R = 1/2; % The Encoding Rate
        msgM = 2; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.BPSKModulator;
        htDeMod = comm.BPSKDemodulator;
        puncpat = -1; % Rate 1/2 Default Rate; No puncture 
    case 1 
        disp('QPSK Rate 1/2')
        modType = 'PSK';
        R = 1/2; % The Encoding Rate
        msgM = 4; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.QPSKModulator('BitInput', true);
        htDeMod = comm.QPSKDemodulator('BitOutput', true);
        puncpat = -1; % Rate 1/2 Default Rate; No puncture 
    case 2
        disp('QPSK Rate 3/4')
        modType = 'PSK';
        R = 3/4; % The Encoding Rate
        msgM = 4; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.QPSKModulator('BitInput', true);
        htDeMod = comm.QPSKDemodulator('BitOutput', true);
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
    case 3 
        disp('16-QAM Rate 1/2')
        modType = 'QAM';
        R = 1/2; % The Encoding Rate
        msgM = 16; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator('ModulationOrder', msgM, 'BitInput', true); % See 22.3.10.9
        htDeMod = comm.RectangularQAMDemodulator('ModulationOrder', msgM, 'BitOutput', true);
        puncpat = -1; % Rate 1/2 Default Rate; No puncture 
    case 4 
        disp('16-QAM Rate 3/4')
        modType = 'QAM';
        R = 3/4; % The Encoding Rate
        msgM = 16; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator('ModulationOrder', msgM, 'BitInput', true); % See 22.3.10.9
        htDeMod = comm.RectangularQAMDemodulator('ModulationOrder', msgM, 'BitOutput', true);
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
    case 5
        disp('64-QAM Rate 2/3')
        modType = 'QAM';
        R = 2/3; % The Encoding Rate
        msgM = 64; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator('ModulationOrder', msgM, 'BitInput', true); % See 22.3.10.9
        htDeMod = comm.RectangularQAMDemodulator('ModulationOrder', msgM, 'BitOutput', true);
        puncpat = [1; 1; 1; 0;]; % Rate 2/3 Figure 18-9
    case 6
        disp('64-QAM Rate 3/4')
        modType = 'QAM';
        R = 3/4; % The Encoding Rate
        msgM = 64; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator('ModulationOrder', msgM, 'BitInput', true); % See 22.3.10.9
        htDeMod = comm.RectangularQAMDemodulator('ModulationOrder', msgM, 'BitOutput', true);
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
    case 7
        disp('64-QAM Rate 5/6')
        modType = 'QAM';
        R = 5/6; % The Encoding Rate
        msgM = 64; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator('ModulationOrder', msgM, 'BitInput', true); % See 22.3.10.9
        htDeMod = comm.RectangularQAMDemodulator('ModulationOrder', msgM, 'BitOutput', true);
        puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
    case 8
        disp('256-QAM Rate 3/4')
        modType = 'QAM';
        R = 3/4; % The Encoding Rate
        msgM = 256; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator('ModulationOrder', msgM, 'BitInput', true); % See 22.3.10.9
        htDeMod = comm.RectangularQAMDemodulator('ModulationOrder', msgM, 'BitOutput', true);
        puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
    case 9
        disp('256-QAM Rate 5/6')
        modType = 'QAM';
        R = 5/6; % The Encoding Rate
        msgM = 256; % The M-ary number, 2 corresponds to binary modulation
        k = log2(msgM);   % # of information bits per symbol
        hMod = comm.RectangularQAMModulator('ModulationOrder', msgM, 'BitInput', true); % See 22.3.10.9
        htDeMod = comm.RectangularQAMDemodulator('ModulationOrder', msgM, 'BitOutput', true);
        puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
    otherwise 
        warning('Unexpected MCS.')
end

if (debug == 0)
    N_Pre_Pad = 0;
    N_Data_Bits = 1e4;
    N_Post_Pad = 0;
    N_Data_Bits = N_Data_Bits + k - mod(N_Data_Bits,k);
    htErrorCalc = comm.ErrorRate; % for debug no encoding
elseif (strcmp(type,'BCC'))
    % Convolutional Encoding Setup
    constlen=7;
    codegen = [171 133]; 
    trellis = poly2trellis(constlen, codegen); % Industry standard 18.3.5.6
    htConvEnc = comm.ConvolutionalEncoder(trellis); 
    htVitDec = comm.ViterbiDecoder(trellis, 'InputFormat', 'Hard'); 
    htVitDec.TracebackDepth = 96;
    htErrorCalc = comm.ErrorRate('ReceiveDelay', htVitDec.TracebackDepth);

    if (puncpat ~= -1) 
        htConvEnc.PuncturePatternSource = 'Property';
        htVitDec.PuncturePatternSource = 'Property';
        htConvEnc.PuncturePattern = puncpat;
        htVitDec.PuncturePattern = puncpat;
    end
    
    % Math Setup for # of bits
    length_param = 4095;  % 4095 is the max LENGTH parameter. See 18.2.2.2
    [numerator,denominator] = rat(R);
    N_DBPS = k;
    N_Scrambler_Init_Bits = 7;
    N_Reserved_Service_Bits = 9;
    N_Tail_bits = 6;
    N_SYM = ceil((N_Scrambler_Init_Bits+N_Reserved_Service_Bits+...
        8*length_param + N_Tail_bits)/lcm(N_DBPS,length(puncpat)))*numerator; % 18.3.5.4
    N_DATA = N_SYM * lcm(N_DBPS,length(puncpat));
    N_PAD = N_DATA - (N_Scrambler_Init_Bits+N_Reserved_Service_Bits+...
        8 * length_param + N_Tail_bits);
    N_Pre_Pad = N_Scrambler_Init_Bits+N_Reserved_Service_Bits;
    N_Post_Pad = N_Tail_bits+N_PAD;
    N_Data_Bits = N_DATA-N_Pre_Pad-N_Post_Pad;
elseif (strcmp(type,'LDPC'))
    % LDPC matrix initialization
    LDPC(0, false, true, 1/2);
    htErrorCalc = comm.ErrorRate;
    
    code_block = 324; % Our LDPC matricies are defined for 324 fixed input
    N_Pre_Pad = 0;
    N_Data_Bits = code_block;
    N_Post_Pad = 0;
end

%% Create a vector to store the BER computed during each iteration
BERVec = zeros(3,length(SNR_Vec)); % Allocate memory to store results
env_c = length(SNR_Vec);

tic;
%waitbars are invalid in parallel pools.
%h = waitbar(0, 'Initializing data cannon...');

% Run the simulation numIter amount of times
% Note that using a parallel pool will not output graphs. 
% Graphs will be generated and can be saved using print
parfor n=1:env_c
  %reset(hErrorCalc)
  %reset(hConvEnc)
  %reset(hVitDec)
  
  hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)', 'SNR', SNR_Vec(n));
  
  hErrorCalc = htErrorCalc.clone;
  if (debug == 0)
  hDeMod = htDeMod;
  elseif (strcmp(type,'BCC'))
  hDeMod = htDeMod;
  hConvEnc = clone(htConvEnc);
  hVitDec = clone(htVitDec);
  elseif (strcmp(type,'LDPC'))
  hDeMod = htDeMod.clone;
  hDeMod.Variance =  1/10^(hChan.SNR/10);
  end
  for i = 1:numIter
    % Generate binary frames of size specified by the frameLength variable
    bits = [zeros(N_Pre_Pad,1);logical(randi([0 1], N_Data_Bits,1));zeros(N_Post_Pad,1)];

    % Interleave the bits   % Not interleaving because parity bit math mess
    txdata = bits; %randintrlv(bits,sum(double('Keenesus'))); 

    % Encode the data
    if (debug == 0)
        encData = txdata; % for debug no encoding
    elseif (strcmp(type,'BCC'))
        encData = step(hConvEnc, txdata); % Conv Enc
    elseif (strcmp(type,'LDPC'))
        encData = LDPC(txdata, true, false, 0); % LDPC Enc
    end
    
    % Modulate the encoded data
    modData = step(hMod, encData);
    
    % Pass the modulated signal through an AWGN channel
    if (strcmp(modType,'PSK'))
        channelOutput = step(hChan, modData);
    elseif (strcmp(type,'LDPC') && debug)
        channelOutput = step(hChan, modData);
    elseif (strcmp(modType,'QAM'))
        channelOutput = awgn(modData, SNR_Vec(n), 'measured'); 
    end
    % Add AWGN, this accounts for 10*log10(R) modification and additional
    % power to the modulation rate. http://www.mathworks.com/examples/matlab-communications/mw/comm-ex70334664-punctured-convolutional-coding
    
    % Demodulate the signal 
    rxsyms = step(hDeMod, channelOutput);

    % Pass the demodulated channel outputs as input to the decoder
    if (debug == 0)
        decData = rxsyms; % for debug no encoding
    elseif (strcmp(type,'BCC'))
        decData = step(hVitDec, rxsyms); % Viterbi dec
    elseif (strcmp(type,'LDPC'))
        decData = LDPC(rxsyms, false, false, 0); % LDPC dec
    end
    
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

%% Compute the theoretical BER for this scenario
figure
berHypo = berawgn(SNR_Vec - 10*log10(k), modType, msgM, 'nondiff');
semilogy(SNR_Vec,berHypo,'r')
hold on
semilogy(SNR_Vec,BERVec(1,:));
hold off
xlabel('SNR (dB)')
legend('Theoretical BER', 'BER');
title(strcat(num2str(msgM), ' ', modType));

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
