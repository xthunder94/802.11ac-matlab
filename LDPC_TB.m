clear all; close all; clc;
format compact;
H = LDPC(1/2);
hEnc = comm.LDPCEncoder(H);
hDec = comm.LDPCDecoder(H);
hMod = comm.RectangularQAMModulator('ModulationOrder', 256, 'BitInput', true);
berawgn(10 - 10*log10(3), 'qam', 256)
for snr = 10
    hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', snr);
    hDemod = comm.RectangularQAMDemodulator('ModulationOrder', 256, 'BitOutput', true);
    hDemod.NormalizationMethod = 'Average power';
    hDemod.AveragePower = 1;
    hDemod.DecisionMethod = 'Approximate log-likelihood ratio';
    hDemod.Variance =  1/10^(hChan.SNR/10);
    hError = comm.ErrorRate;
    for counter = 1:1
      data           = logical(randi([0 1], 324, 1));
      encodedData    = step(hEnc, data);
      modSignal      = step(hMod, encodedData);
      scatterplot(modSignal)
      receivedSignal = step(hChan, modSignal);
      demodSignal    = step(hDemod, receivedSignal);
      receivedBits   = step(hDec, demodSignal);
      errorStats     = step(hError, data, receivedBits);
    end
    fprintf('Error rate       = %1.2f\nNumber of errors = %d\n', ...
      errorStats(1), errorStats(2))
end