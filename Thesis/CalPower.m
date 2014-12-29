% calculate transmission power given a transmitter and a set of interferers
function [power] = CalPower(Tx,Rx,vec_interferer,vec_interfererPower,mat_Gij,vec_Fi,bandwidth,N0)
    power = 0;
    if length(vec_interferer) == 0
        power = (vec_Fi(Tx)*bandwidth*N0)/mat_Gij(Tx,Rx);
    else
        p = 0;
        for i = 1:length(vec_interferer)
            p = p + vec_interfererPower(i)*mat_Gij(vec_interferer(i),Rx);
        end
        power = (vec_Fi(Tx)*(bandwidth*N0 + p))/mat_Gij(Tx,Rx);
    end 
end