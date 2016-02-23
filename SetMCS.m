function [modType, M, k, R, puncpat hMod, hDemod] = SetMCS(MCS)

    % Reference Materials
    %{
    M = [2 4 16 64 256];                        % the M-ary number
    k = log2(M);                                % # of information bits per symbol
    R = [1/2 3/4 5/6];                          % the encoding rate
    puncpat = -1;                               % Rate 1/2  Default Rate; no puncture 
    puncpat = [1; 1; 1; 0;];                    % Rate 2/3  Figure 18-9
    puncpat = [1; 1; 1; 0; 0; 1;];              % Rate 3/4  Figure 18-9
    puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;];	% Rate 5/6  Figure 20-11
    %}
    
    % Choosing Modulation and Coding Scheme
    switch MCS
        
        case 0
            disp('BPSK Rate 1/2')
            modType = 'PSK';
            M = 2;
            k = log2(M);
            R = 1/2;
            puncpat = -1;
            hMod = comm.BPSKModulator;
            hDemod = comm.BPSKDemodulator;
            
        case 1
            disp('QPSK Rate 1/2')
            modType = 'PSK';
            M = 4;
            k = log2(M);
            R = 1/2;
            puncpat = -1;
            hMod = comm.QPSKModulator('BitInput', true);
            hDemod = comm.QPSKDemodulator('BitOutput', true);
  
        case 2
            disp('QPSK Rate 3/4')
            modType = 'PSK';
            M = 4;
            k = log2(M);
            R = 3/4;
            puncpat = [1; 1; 1; 0; 0; 1;];
            hMod = comm.QPSKModulator('BitInput', true);
            hDemod = comm.QPSKDemodulator('BitOutput', true);
            
        case 3 
            disp('16-QAM Rate 1/2')
            modType = 'QAM';
            M = 16;
            k = log2(M);
            R = 1/2;
            puncpat = -1;
            hMod = comm.RectangularQAMModulator('ModulationOrder',M, 'BitInput', true); % See 22.3.10.9
            hDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 4 
            disp('16-QAM Rate 3/4')
            modType = 'QAM';
            M = 16;
            k = log2(M);
            R = 3/4;
            puncpat = [1; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            hDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 5
            disp('64-QAM Rate 2/3')
            modType = 'QAM';
            M = 64;
            k = log2(M);
            R = 2/3;
            puncpat = [1; 1; 1; 0;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            hDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 6
            disp('64-QAM Rate 3/4')
            modType = 'QAM';
            M = 64;
            k = log2(M);
            R = 3/4;
            puncpat = [1; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            hDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);

        case 7
            disp('64-QAM Rate 5/6')
            modType = 'QAM';
            M = 64;
            k = log2(M);
            R = 5/6;
            puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            hDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 8
            disp('256-QAM Rate 3/4')
            modType = 'QAM';
            M = 256;
            k = log2(M);
            R = 3/4;
            puncpat = [1; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            hDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 9
            disp('256-QAM Rate 5/6')
            modType = 'QAM';
            M = 256;
            k = log2(M);
            R = 5/6;
            puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            hDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);

        otherwise 
            warning('Unexpected MCS.')
            
    end
end




