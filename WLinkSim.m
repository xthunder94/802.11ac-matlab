% Functions that runs simulations for a single modulation scheme

function WLinkSim(MCS_in, SNR_Vec_in, encType_in, debug_in, numIter_in, h_in)


    %%%%%%%%%%%%%%%%%%%%%%%%%
    %        INPUTS         %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    % Specify modulation and coding scheme (ranges from 0-9)
    MCS = MCS_in;
    %{
    0: BPSK Rate 1/2
    1: QPSK Rate 1/2
    2: QPSK Rate 3/4
    3: 16-QAM Rate 1/2
    4: 16-QAM Rate 3/4
    5: 64-QAM Rate 2/3
    6: 64-QAM Rate 3/4
    7: 64-QAM Rate 5/6
    8: 256-QAM Rate 3/4
    9: 256-QAM Rate 5/6
    %}
    
    SNR_Vec = SNR_Vec_in;
    encType = encType_in;
    debug = debug_in;
    numIter = numIter_in;
    h = h_in;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %       ENCODING        %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    % Set modulation and coding scheme
    [display, modType, lSpec, M, k, R, k_TCB, puncpat, ...
        hMod, htDemod] = SetMCS(MCS, encType, debug);

    % Set encoding method (BCC or LDPC, debug or no debug)
    [htEnc, htDec, htErrorCalc, ...
        N_Pre_Pad, N_Data_Bits, N_Post_Pad] = ...
        SetEncoder(encType, debug, k, R, k_TCB, puncpat);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %  SIMULATION & GRAPHS  %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Run simulation and retrieve BERs
    [ber, berHypo] = Simulation(numIter, SNR_Vec, encType, debug, ...
        modType, k, R, M, hMod, htDemod, ...
        htEnc, htDec, htErrorCalc, ...
        N_Pre_Pad, N_Data_Bits, N_Post_Pad)
    
    
    % Graph theoretical or actual BERs
    
    name = strcat(display,{' (MCS '}, num2str(MCS), ')');
    nameTh = strcat({'Theoretical (MCS '}, num2str(MCS), ')');
    
    if(debug == 0)
        semilogy(SNR_Vec, berHypo, 'k', 'DisplayName', nameTh{1});
        hold on
        semilogy(SNR_Vec, ber, lSpec, 'DisplayName', name{1}); % for debug
        hold on
    else
        semilogy(SNR_Vec, ber, lSpec, 'DisplayName', name{1});
        hold on
    end
    
    legend('off')
    legend('show')
    
    
    % Update waitbar
    waitbar(MCS/9,h,sprintf('Evaluating MCS %d',MCS));

end
