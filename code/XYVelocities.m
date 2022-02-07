function [Xdir,Ydir] = XYVelocities(currVel,currDir,numElec)
%XYVelocities Calculate X and Y components of velocity vector
%   Take column vector of current electron velocities and directions,
%   calculate the X and Y velocities vectors for all electrons from given
%   magnitude and angle

for n = 1:numElec
    
    Xdir(n,1) = currVel(n,1)*cos(currDir(n,1));
    Ydir(n,1) = currVel(n,1)*sin(currDir(n,1));
    
end

end

