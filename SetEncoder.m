function [htConvEnc, htVitDec, htErrorCalc, ...
    N_Punc_Pad, N_Pre_Pad, N_Data_Bits, N_Post_Pad] = SetEncoder(encType, debug, k, puncpat)

    if (debug == 0)
        htConvEnc = 0; % default
        htVitDec = 0; % default
        htErrorCalc = comm.ErrorRate; % for debug no encoding
        N_Punc_Pad = 0; % default
        N_Pre_Pad = 0;
        N_Data_Bits = 1e4;
        N_Post_Pad = 0;
        N_Data_Bits = N_Data_Bits + k - mod(N_Data_Bits,k);

    elseif (strcmp(encType, 'BCC'))
        [htConvEnc, htVitDec, htErrorCalc, ...
            N_Punc_Pad, N_Pre_Pad, N_Data_Bits, N_Post_Pad] = BCC(k, puncpat);

    elseif (strcmp(encType,'LDPC'))
        LDPC(0, false, true, 1/2); % LDPC matrix initialization
        htConvEnc = 0; % default
        htVitDec = 0; % default
        htErrorCalc = comm.ErrorRate;
        code_block = 324; % Our LDPC matrices are defined for 324 fixed input
        N_Punc_Pad = 0; % default
        N_Pre_Pad = 0;
        N_Data_Bits = code_block;
        N_Post_Pad = 0;

    end
end
