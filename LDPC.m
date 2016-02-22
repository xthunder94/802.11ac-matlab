function [bitsOutput] = LDPC(bitsInput, bEncode)
    % LDPC Encodes LDPC According to 802.11AC Standard
    %    802.11AC Optional LDPC Encoding
    persistent H hEnc hDec;
    if isempty(H)
        H = HTLDPC(1/2);
        hEnc = comm.LDPCEncoder(H);
        hDec = comm.LDPCDecoder(H);
    end
    bitsOutput1 = step(hEnc, bitsInput);
    bitsOutput1 = double(bitsOutput1);
    bitsOutput1(bitsOutput1 == 1) = -1;
    bitsOutput1(bitsOutput1 == 0) = 1;
    bitsOutput = step(hDec, bitsOutput1);
%%%%%%%%%%%%%%%%%%%%%%%% THIS IS THE TESTING CODE PLS CHANGERINO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hMod = comm.RectangularQAMModulator('ModulationOrder',256,'BitInput',true);
    hChan = comm.AWGNChannel(...
            'NoiseMethod','Signal to noise ratio (SNR)','SNR',-5);
    hDemod = comm.RectangularQAMDemodulator('ModulationOrder',256, 'BitOutput',true,...
            'DecisionMethod','Approximate log-likelihood ratio', ...
            'Variance', 1/10^(hChan.SNR/10));
    hError = comm.ErrorRate;
    for counter = 1:10000
      data           = logical(randi([0 1], 324, 1));
      encodedData    = step(hEnc, data);
      modSignal      = step(hMod, encodedData);
      receivedSignal = step(hChan, modSignal);
      demodSignal    = step(hDemod, receivedSignal);
      receivedBits   = step(hDec, demodSignal);
      errorStats     = step(hError, data, receivedBits);
    end
    fprintf('Error rate       = %1.2f\nNumber of errors = %d\n', ...
      errorStats(1), errorStats(2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [C] = B(i)
    % Subblock Size Z=27
    C = circshift(eye(27), i, 2);
end

function [C] = E()
    % Empty Subblock Z=27
    C = zeros(27);
end

function [H] = HTLDPC(r)
    if r == 1/2
        H = [B(00) E     E     E     B(00) B(00) E     E     B(00) E     E     B(00) B(01) B(00) E     E     E     E     E     E     E     E     E     E    ; ...
             B(22) B(00) E     E     B(17) E     B(00) B(00) B(12) E     E     E     E     B(00) B(00) E     E     E     E     E     E     E     E     E    ; ...
             B(06) E     B(00) E     B(10) E     E     E     B(24) E     B(00) E     E     E     B(00) B(00) E     E     E     E     E     E     E     E    ; ...
             B(02) E     E     B(00) B(20) E     E     E     B(25) B(00) E     E     E     E     E     B(00) B(00) E     E     E     E     E     E     E    ; ...
             B(23) E     E     E     B(03) E     E     E     B(00) E     B(09) B(11) E     E     E     E     B(00) B(00) E     E     E     E     E     E    ; ...
             B(24) E     B(23) B(01) B(17) E     B(03) E     B(10) E     E     E     E     E     E     E     E     B(00) B(00) E     E     E     E     E    ; ...
             B(25) E     E     E     B(08) E     E     E     B(07) B(18) E     E     B(00) E     E     E     E     E     B(00) B(00) E     E     E     E    ; ...
             B(13) B(24) E     E     B(00) E     B(08) E     B(06) E     E     E     E     E     E     E     E     E     E     B(00) B(00) E     E     E    ; ...
             B(07) B(20) E     B(16) B(22) B(10) E     E     B(23) E     E     E     E     E     E     E     E     E     E     E     B(00) B(00) E     E    ; ...
             B(11) E     E     E     B(19) E     E     E     B(13) E     B(03) B(17) E     E     E     E     E     E     E     E     E     B(00) B(00) E    ; ...
             B(25) E     B(08) E     B(23) B(18) E     B(14) B(09) E     E     E     E     E     E     E     E     E     E     E     E     E     B(00) B(00); ...
             B(03) E     E     E     B(16) E     E     B(02) B(25) B(05) E     E     B(01) E     E     E     E     E     E     E     E     E     E     B(00); ...
             ];
        H = sparse(H);
    end
end