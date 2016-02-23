function [ber, berHypo] = Simulation(numIter, SNR_Vec, encType, debug, ...
    modType, k, M, hMod, hDemod, ...
    htConvEnc, htVitDec, htErrorCalc, ...
    N_Pre_Pad, N_Data_Bits, N_Post_Pad)

    % Create a vector to store the BER computed during each iteration
    BER_Vec = zeros(3, length(SNR_Vec)); % Allocate memory to store results
    env_c = length(SNR_Vec);

    tic;
    %waitbars are invalid in parallel pools.
    %h = waitbar(0, 'Initializing data cannon...');

    % Run the simulation numIter amount of times
    % Note that using a parallel pool will not output graphs. 
    % Graphs will be generated and can be saved using print
    for n = 1:env_c
        
        %reset(hErrorCalc)
        %reset(hConvEnc)
        %reset(hVitDec)
        hErrorCalc = htErrorCalc.clone;
        if(~strcmp(encType, 'LDPC'))
            hConvEnc = htConvEnc.clone;
            hVitDec = htVitDec.clone;
        end
        hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)', 'SNR', SNR_Vec(n));
        
        for i = 1:numIter
            
            % Generate binary frames of size specified by the frameLength variable
            bits = [zeros(N_Pre_Pad,1); logical(randi([0 1], N_Data_Bits,1)); zeros(N_Post_Pad,1)];
            
            % Interleave the bits   % Not interleaving because parity bit math mess
            txdata = bits; %randintrlv(bits,sum(double('Keenesus')));
            
            % Encode the data
            if (debug == 0)
                encData = txdata; % for debug no encoding
            elseif (strcmp(encType, 'BCC'))
                encData = step(hConvEnc, txdata); % Conv Enc
            elseif (strcmp(encType, 'LDPC'))
                encData = LDPC(txdata, true, false, 0); % LDPC Enc
            end
            
            % Modulate the encoded data
            modData = step(hMod, encData);
            
            % Pass the modulated signal through an AWGN channel
            if (strcmp(modType, 'PSK'))
                channelOutput = step(hChan, modData);
            elseif (strcmp(modType, 'QAM'))
                channelOutput = awgn(modData, SNR_Vec(n), 'measured');
            end
            % Add AWGN, this accounts for 10*log10(R) modification and additional
            % power due to the modulation rate. http://www.mathworks.com/examples/matlab-communications/mw/comm-ex70334664-punctured-convolutional-coding
            
            % Demodulate the signal
            rxsyms = step(hDemod, channelOutput);
            
            % Pass the demodulated channel outputs as input to the decoder
            if (debug == 0)
                decData = rxsyms; % for debug no encoding
            elseif (strcmp(encType,'BCC'))
                decData = step(hVitDec, rxsyms); % Viterbi dec
            elseif (strcmp(encType,'LDPC'))
                decData = LDPC(rxsyms, false, false, 0); % LDPC dec
            end
            
            % Deinterleave the bits % Not interleaving because parity bit math mess
            data = decData; %randdeintrlv(decData,sum(double('Keenesus')));
            
            % Compute and accumulate errors
            BER_Vec(:,n) = step(hErrorCalc, bits, data);
        end
        
        %waitbar(n/env_c,h,sprintf('Evaluating %d',EbNoEncoderOutput(n)))
    end % End iteration
    
    toc;
    %close(h);

    % Extract the BERs
    ber = BER_Vec(1,:);
    
    % Compute the theoretical BERs for this scenario
    berHypo = berawgn(SNR_Vec - 10*log10(k), modType, M, 'nondiff');
    
end

    
