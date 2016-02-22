function [BER_Vec] = BCC(SNR_Vec, MCS)

    switch MCS
        case 0 
            disp('BPSK Rate 1/2')
            R = 1/2;
            msgM = 2;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.BPSKModulator;
            hDeMod = comm.BPSKDemodulator;
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(R); %fix??
        case 1 
            disp('QPSK Rate 1/2')
            R = 1/2;
            msgM = 4;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.QPSKModulator;
            hDeMod = comm.QPSKDemodulator;
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(1/2);
        case 2
            disp('QPSK Rate 3/4')
            R = 3/4;
            msgM = 4;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.QPSKModulator;
            hDeMod = comm.QPSKDemodulator;
            puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(3/4);
        case 3 
            disp('16-QAM Rate 1/2')
            R = 1/2;
            msgM = 16;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
            hDeMod = comm.RectangularQAMDemodulator(msgM);
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(1/2);
        case 4 
            disp('16-QAM Rate 3/4')
            R = 3/4;
            msgM = 16;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
            hDeMod = comm.RectangularQAMDemodulator(msgM);
            puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(3/4);
        case 5
            disp('64-QAM Rate 2/3')
            R = 2/3;
            msgM = 64;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
            hDeMod = comm.RectangularQAMDemodulator(msgM);
            puncpat = [1; 1; 1; 0;]; % Rate 2/3 Figure 18-9
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(2/3);
        case 6
            disp('64-QAM Rate 3/4')
            R = 3/4;
            msgM = 64;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
            hDeMod = comm.RectangularQAMDemodulator(msgM);
            puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(3/4);
        case 7
            disp('64-QAM Rate 5/6')
            R = 5/6;
            msgM = 64;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
            hDeMod = comm.RectangularQAMDemodulator(msgM);
            puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(5/6);
        case 8
            disp('256-QAM Rate 3/4')
            R = 3/4;
            msgM = 256;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
            hDeMod = comm.RectangularQAMDemodulator(msgM);
            puncpat = [1; 1; 1; 0; 0; 1;]; % Rate 3/4  Figure 18-9
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(3/4);
        case 9
            disp('256-QAM Rate 5/6')
            R = 5/6;
            msgM = 256;
            k = log2(msgM);   % # of information bits per symbol
            hMod = comm.RectangularQAMModulator(msgM); % See 22.3.10.9
            hDeMod = comm.RectangularQAMDemodulator(msgM);
            puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;]; % Rate 5/6  Figure 20-11
            EbNoEncoderOutput = SNR_Vec - 10*log10(k/msgM) + 10*log10(5/6);
        otherwise 
            warning('Unexpected MCS.')
    end

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
    N_Punc_Pad = length(puncpat) - N_PAD - mod(N_DATA,length(puncpat)); % Padding for computation

    constlen=7;
    codegen = [171 133]; 
    trellis = poly2trellis(constlen, codegen); % Industry standard 18.3.5.6
    htConvEnc = comm.ConvolutionalEncoder(trellis); 
    htVitDec = comm.ViterbiDecoder(trellis, 'InputFormat', 'Hard'); 
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
    BER_Vec = zeros(3,length(EbNoEncoderOutput)); % Allocate memory to store results
    env_c = length(EbNoEncoderOutput);
    frameLength = N_DATA;
    numIter = 10 %1e6;

    % ScatterPlot! 
    %sPlotFig = scatterplot(qammod(0:msgM-1,msgM,0,'gray'),1,0,'k*');

    tic;
    %waitbars are invalid in parallel pools.
    %h = waitbar(0, 'Initializing data cannon...');

    % Run the simulation numIter amount of times
    % Note that using a parallel pool will not output graphs. 
    % Graphs will be generated and can be saved using print
    %parfor 
    for n=1:env_c
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
        bits = [zeros(N_Scrambler_Init_Bits+N_Reserved_Service_Bits,1); ...
            randi([0 1], frameLength-N_Scrambler_Init_Bits-N_Reserved_Service_Bits, 1); zeros(N_PAD+N_Punc_Pad,1)];
        % Interleave the bits 
        txdata = bits; %randintrlv(bits,sum(double('Keenesus')));
        % Convolutionally encode the data
        encData = step(hConvEnc, txdata);
        % Modulate the encoded data
        modData = step(hMod, encData);
        % Pass the modulated signal through an AWGN channel
        channelOutput = step(hChan, modData);
        % Demodulate the signal 
        demodData = step(hDeMod, channelOutput);
        % Pass the real part of the channel complex outputs as the unquantized
        % input to the Viterbi decoder.
        decData = step(hVitDec, demodData);
        % Deinterleave the bits
        data = decData; %randdeintrlv(decData,sum(double('Keenesus')));
        % Compute and accumulate errors
        BER_Vec(:,n) = step(hErrorCalc, bits, data);
      end

      %waitbar(n/env_c,h,sprintf('Evaluating %d',EbNoEncoderOutput(n)))
    end % End iteration
    toc;
    %close(h);

end
