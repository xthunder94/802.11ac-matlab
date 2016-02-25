% Helper function to simulate for a provided number of iterations, returns
% both actual and hypothetical BER vectors

function [ber, berHypo] = Simulation(numIter, SNR_Vec, encType, debug, ...
    modType, k, R, M, hMod, htDemod, ...
    htEnc, htDec, htErrorCalc, ...
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
    parfor n = 1:env_c
        
        %reset(hErrorCalc)
        %reset(hEnc)
        %reset(hDec)
        
        hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)', 'SNR', SNR_Vec(n));
    
        hErrorCalc = clone(htErrorCalc);
        hDemod = clone(htDemod);
        
        if(debug ~= 0) % for encoding
            hEnc = clone(htEnc);
            hDec = clone(htDec);                
            if (strcmp(encType,'LDPC'))
                hDemod.DecisionMethod = 'Approximate log-likelihood ratio';
                hDemod.Variance =  1/10^(hChan.SNR/10);
                if(~((M == 2) || (M == 4))) % NormalizationMethod doesn't exist for BPSK or QPSK
                    hDemod.NormalizationMethod = 'Average power';
                    hDemod.AveragePower = 1;
                end
            end
        end
        
        for i = 1:numIter
            
            % Generate binary frames of size specified by the frameLength variable
            bits = [zeros(N_Pre_Pad,1); logical(randi([0 1], N_Data_Bits,1)); zeros(N_Post_Pad,1)];
            
            % Interleave the bits   % Not interleaving because parity bit math mess
            txdata = bits; %randintrlv(bits,sum(double('Keenesus')));
            
            % Encode the data
            if (debug == 0)
                encData = txdata; % for debug no encoding
            else
                encData = step(hEnc, txdata); % for encoding
            end
            
            % Modulate the encoded data
            modData = step(hMod, encData);
            
            % Pass the modulated signal through an AWGN channel
            if (strcmp(modType,'PSK') || (strcmp(encType,'LDPC') && debug))
                channelOutput = step(hChan, modData);
            elseif (strcmp(modType,'QAM'))
                channelOutput = awgn(modData, SNR_Vec(n), 'measured');
            end
            % Add AWGN, this accounts for 10*log10(R) modification and additional
            % power due to the modulation rate. http://www.mathworks.com/examples/matlab-communications/mw/comm-ex70334664-punctured-convolutional-coding
            
            % Demodulate the signal
            rxsyms = step(hDemod, channelOutput);
            
            % Pass the demodulated channel outputs as input to the decoder
            if (debug == 0)
                decData = rxsyms; % for debug no encoding
            else
                decData = step(hDec, rxsyms); % for decoding
            end
            
            % Deinterleave the bits % Not interleaving because parity bit math mess
            data = decData; %randdeintrlv(decData,sum(double('Keenesus')));
            
            % Compute and accumulate errors
            BER_Vec(:,n) = step(hErrorCalc, bits, double(data));
        end
        
        %waitbar(n/env_c,h,sprintf('Evaluating %d',EbNoEncoderOutput(n)))
    end % End iteration
    
    toc;
    %close(h);

    % Extract the BERs
    ber = BER_Vec(1,:);
    
    % Compute the theoretical BERs for this scenario
    berHypo = berawgn(SNR_Vec - 10*log10(k*R), modType, M, 'nondiff');
    
end
