function [htConvEnc, htVitDec, htErrorCalc, ...
    N_Punc_Pad, N_Pre_Pad, N_Data_Bits, N_Post_Pad] = BCC(k, puncpat)

    % Convolutional Encoding Setup
    constlen = 7;
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
    N_DBPS = k;
    N_Scrambler_Init_Bits = 7;
    N_Reserved_Service_Bits = 9;
    N_Tail_bits = 6;
    N_SYM = ceil((N_Scrambler_Init_Bits + N_Reserved_Service_Bits + ...
        8*length_param+N_Tail_bits)/N_DBPS); % 18.3.5.4
    N_DATA = N_SYM * N_DBPS;
    N_PAD = N_DATA - (N_Scrambler_Init_Bits + N_Reserved_Service_Bits + ...
        8 * length_param + N_Tail_bits);
    
    N_Punc_Pad = mod(length(puncpat) - N_PAD - mod(N_DATA,length(puncpat)), ...
        length(puncpat)); % Padding for computation
    N_Punc_Pad = N_Punc_Pad + mod(N_DATA + N_Punc_Pad + N_PAD + N_Tail_bits, ...
        lcm(k,length(puncpat)));
    N_Pre_Pad = N_Scrambler_Init_Bits + N_Reserved_Service_Bits;
    N_Data_Bits = N_DATA - N_Scrambler_Init_Bits - N_Reserved_Service_Bits;
    N_Post_Pad = N_Tail_bits + N_PAD + N_Punc_Pad;
    
end
    