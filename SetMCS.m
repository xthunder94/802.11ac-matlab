% Helper function that sets modulation and coding scheme and related
% parameters based on inputs

function [display, modType, lSpec, M, k, R, k_TCB, puncpat, ...
    hMod, htDemod] = SetMCS(MCS, encType, debug)

    % Reference Materials
    %{
    M = [2 4 16 64 256];                        % the M-ary number
    k = log2(M);                                % # of information bits per symbol
    R = [1/2 3/4 5/6];                          % the encoding rate
    k_TCB = [5 7.5 10 15]                       % traceback depth constant
    puncpat = -1;                               % Rate 1/2  Default Rate; no puncture 
    puncpat = [1; 1; 1; 0;];                    % Rate 2/3  Figure 18-9
    puncpat = [1; 1; 1; 0; 0; 1;];              % Rate 3/4  Figure 18-9
    puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;];	% Rate 5/6  Figure 20-11
    %}

    % Choosing Modulation and Coding Scheme
    switch MCS
        
        case 0
            display = 'BPSK Rate 1/2';
            modType = 'PSK';
            lSpec = '*r-';
            M = 2;
            k = log2(M);
            R = 1/2;
            k_TCB = 5;
            puncpat = -1;
            hMod = comm.BPSKModulator;
            htDemod = comm.BPSKDemodulator;
            
        case 1
            display = 'QPSK Rate 1/2';
            modType = 'PSK';
            lSpec = '*g-';
            M = 4;
            k = log2(M);
            R = 1/2;
            k_TCB = 5;
            puncpat = -1;
            hMod = comm.QPSKModulator('BitInput', true);
            htDemod = comm.QPSKDemodulator('BitOutput', true);
            
        case 2
            display = 'QPSK Rate 3/4';
            modType = 'PSK';
            lSpec = 'xg-';
            M = 4;
            k = log2(M);
            R = 3/4;
            k_TCB = 10;
            puncpat = [1; 1; 1; 0; 0; 1;];
            hMod = comm.QPSKModulator('BitInput', true);
            htDemod = comm.QPSKDemodulator('BitOutput', true);
            
        case 3 
            display = '16-QAM Rate 1/2';
            modType = 'QAM';
            lSpec = '*c-';
            M = 16;
            k = log2(M);
            R = 1/2;
            k_TCB = 5;
            puncpat = -1;
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            htDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 4 
            display = '16-QAM Rate 3/4';
            modType = 'QAM';
            lSpec = 'xc-';
            M = 16;
            k = log2(M);
            R = 3/4;
            k_TCB = 10;
            puncpat = [1; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            htDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 5
            display = '64-QAM Rate 2/3';
            modType = 'QAM';
            lSpec = '+b-';
            M = 64;
            k = log2(M);
            R = 2/3;
            k_TCB = 7.5;
            puncpat = [1; 1; 1; 0;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            htDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 6
            display = '64-QAM Rate 3/4';
            modType = 'QAM';
            lSpec = 'xb-';
            M = 64;
            k = log2(M);
            R = 3/4;
            k_TCB = 10;
            puncpat = [1; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            htDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);

        case 7
            display = '64-QAM Rate 5/6';
            modType = 'QAM';
            lSpec = '.b-';
            M = 64;
            k = log2(M);
            R = 5/6;
            k_TCB = 15;
            puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            htDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 8
            display = '256-QAM Rate 3/4';
            modType = 'QAM';
            lSpec = 'xm-';
            M = 256;
            k = log2(M);
            R = 3/4;
            k_TCB = 10;
            puncpat = [1; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            htDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);
            
        case 9
            display = '256-QAM Rate 5/6';
            modType = 'QAM';
            lSpec = '.m-';
            M = 256;
            k = log2(M);
            R = 5/6;
            k_TCB = 15;
            puncpat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1;];
            hMod = comm.RectangularQAMModulator('ModulationOrder', M, 'BitInput', true); % See 22.3.10.9
            htDemod = comm.RectangularQAMDemodulator('ModulationOrder', M, 'BitOutput', true);

        otherwise 
            warning('Unexpected MCS.')
    end
    
    % Configure LDPC moderator to use average power
    if(~((M == 2) || (M == 4)) && strcmp(encType, 'LDPC')) % NormalizationMethod doesn't exist for BPSK or QPSK
        hMod.NormalizationMethod = 'Average power';
        hMod.AveragePower = 1;
    end
    
    if (debug == 0)
        R = 1; % No encoding
    end

end
