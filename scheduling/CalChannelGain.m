function G = CalChannelGain(x1,y1,x2,y2)
    n0 = 1e-16; % = pow(10,-17.4)/1000;Watt, ref:Wikipwdia Thermal Noise Steven: 1e-19
    pathLoss0_MacroUE = 131.1; % no unit ref:LTE release 9
    pathLossAlpha_MacroUE = 42.68; % no unit ref:LTE release 9
    pathLoss0_MacroUE_Gib = 131.1; % For simulation Adjustment
    pathLossAlpha_MacroUE_Gib = 42.68;
    
    tempDib = (x1-x2)^2 + (y1-y2)^2;
    pathLossDBMacroUE = pathLoss0_MacroUE + pathLossAlpha_MacroUE * (0.5*log10(tempDib) - 3);
    G = 10^( -(pathLossDBMacroUE/10) );
end