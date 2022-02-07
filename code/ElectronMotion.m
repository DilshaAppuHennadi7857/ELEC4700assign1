function ElectronMotion(L, W, Temp)
%ElectronMotion Model the motion of electrons in the region
%   First it randomly assign a place in the region for each electron. Then
%   it assigns each electron with an initial velocity

global C
global m_n x y Vx Vy numAtoms nAtoms

% Assign an x,y coordinate to each electron such that they are within the
% region.
for i = 1:numAtoms
    
   x(i) = (L - 0)*rand;
   y(i) = (W - 0)*rand;
    
end

x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5);
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5);

% Assign each particle with an initial velocity
if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / m_n);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms));
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms));

nAtoms = nAtoms + numAtoms;

end

