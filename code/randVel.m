function [velocity] = randVel(Temp,v_th,numElec,C)
%RANDVEL Create column vector of random velocities
%   

dummyVar = C.kb*Temp/C.m_0;

velocity = v_th + sqrt(dummyVar).*randn(numElec,1);

% velocity = (2*dummyVar*log(n*(2*pi*dummyVar)^(1/2))).^(1/2);

end

