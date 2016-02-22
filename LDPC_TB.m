clear all; close all; clc;
format compact;
LDPC(0, false, true, 1/2);
for qam = [4 16 32 64 128 256]
    hMod = comm.RectangularQAMModulator('ModulationOrder', qam, 'BitInput', true);
    for snr = 0:16
        hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', snr);
        hDemod = comm.RectangularQAMDemodulator('ModulationOrder', qam, 'BitOutput', true, ...
                'DecisionMethod', 'Approximate log-likelihood ratio', ...
                'Variance', 1/10^(hChan.SNR/10));
        hError = comm.ErrorRate;
        for counter = 1:10000
          data           = logical(randi([0 1], 324, 1));
          encodedData    = LDPC(data, true, false, 0);
          modSignal      = step(hMod, encodedData);
          receivedSignal = step(hChan, modSignal);
          demodSignal    = step(hDemod, receivedSignal);
          receivedBits   = LDPC(demodSignal, false, false, 0);
          errorStats     = step(hError, data, receivedBits);
        end
        fprintf('Error rate       = %1.2f\nNumber of errors = %d\n', ...
          errorStats(1), errorStats(2))
    end
end