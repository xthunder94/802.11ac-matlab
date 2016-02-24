clear all; close all; clc;
format compact;
H = LDPC(1/2);
hEnc = comm.LDPCEncoder(H);
hDec = comm.LDPCDecoder(H);
hMod = comm.RectangularQAMModulator('ModulationOrder', 16, 'BitInput', true);
berawgn(10 - 10*log10(3), 'qam', 16)
for snr = 10
    hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', snr);
    hDemod = comm.RectangularQAMDemodulator('ModulationOrder', 16, 'BitOutput', true);
    hError = comm.ErrorRate;
    for counter = 1:1000
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
end