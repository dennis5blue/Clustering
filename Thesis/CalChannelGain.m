function G = CalChannelGain(x1,y1,x2,y2)
    % reference from TR36.888    
    tempDib = (((x1-x2)^2 + (y1-y2)^2)^0.5)*1e-3; %km
    pathLossDBMacroUE = 120.9 + 37.6*(log10(tempDib)); % operate at 900MHz
    G = 10^( -(pathLossDBMacroUE/10) );
end