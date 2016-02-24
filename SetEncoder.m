function [htEnc, htDec, htErrorCalc, ...
    N_Pre_Pad, N_Data_Bits, N_Post_Pad] = SetEncoder(encType, debug, k, R, k_TCB, puncpat)

    if (debug == 0)
        
        htEnc = 0; % default
        htDec = 0; % default
        htErrorCalc = comm.ErrorRate; % for debug no encoding
        
        N_Pre_Pad = 0;
        N_Data_Bits = 1e4;
        N_Post_Pad = 0;
        N_Data_Bits = N_Data_Bits + k - mod(N_Data_Bits,k);

    elseif (strcmp(encType, 'BCC'))
        [htEnc, htDec, htErrorCalc, ...
             N_Pre_Pad, N_Data_Bits, N_Post_Pad] = BCC(k, R, k_TCB, puncpat);

    elseif (strcmp(encType,'LDPC'))
        
        % LDPC matrix initialization
        H = LDPC(R);
        htEnc = comm.LDPCEncoder(H);
        htDec = comm.LDPCDecoder(H);
        htErrorCalc = comm.ErrorRate;
        
        code_block = 648 * R; % Our LDPC matricies are defined for 648 * R fixed input
        N_Pre_Pad = 0;
        N_Data_Bits = code_block;
        N_Post_Pad = 0;

    end
end
