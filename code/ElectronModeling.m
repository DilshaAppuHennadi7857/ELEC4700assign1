%%
% 
% 1. Thermal velocity, assuming a temperature of T = 300K

Temp = 300;
v_th = sqrt(C.kb * Temp / m_n);

%%
%
% Given a temperature of 300K, the thermal velocity is calculated to be
% about 1.3224e5 m/s

%%
%
% 2. If the mean time between collisions is t_mn = 0.2 ps, what is the mean
% free path?

t_mn = 0.2e-12;
L_n = v_th*t_mn; %m

%%
%
% The mean free path is taken as the thermal velosity times the mean time
% between collisions, which gives us a distance of 2.6449e-8 m.

%%
%
%
% 3. Write a program to model the random motion of electrons.

currX = (nomRegionL).*rand(numElec, 1); % set random initial x position
currY = (nomRegionW).*rand(numElec, 1); % set random initial y position
currVel = zeros(numElec,1); % initialize column vector to hold velocity, velocities initialized with zero
currVel(:,1) = v_th; % set all values in the first column as v_th
currDir = (2*pi).*rand(numElec,1); % create column vector with current direction of each electron
currVX = []; % set up x and y velocity vectors as empity column vector
currVY = [];

% Calculate initial vx and vy for each electron given their direction
[currVX, currVY] = XYVelocities(currVel,currDir,numElec);

%%
%
% At this point we should have the initial x and y locations of all
% electron. The locations should be within the bounded region. All
% electrons should have their initial velocity set to the thermal velocity,
% v_th and a random direction between 0 to 2*pi. From the magnitude of the
% velocity and the direction, the x and y vector components of velocity can
% be determined. From here, we can now see how the electrons will move
% about the region. To plot trajectories, we need to save the
% previous locations of each electrons - i.e need a matrix where rows are
% individual electrons and columns are different times - these will be
% saved in a separate matric from the currXXX column vectors.

for n = 0:numTimeStep
    
    currTime(n+1) = n*dt; % determine current time (fs)
    
    % Update according to Newton's laws of motion
    % new position = old position + velocity*time
    
    if n > 0 %update position after t=0
        
        currX(:,1) = currX(:,1) + currVX(:,1)*dt; % calculate new X
        % if e- crosses right bound, set X=0; if e- crosses left bound, set
        % X=upper bound
        crossRight = currX > nomRegionL;
        currX(crossRight) = 0;
        crossLeft = currX < 0;
        currX(crossLeft) = nomRegionL;
        
        newY = currY + currVY*dt; % check if new Y crosses boundary
        bounce = (newY>nomRegionW) | (newY<0); % bounce electron if it hits bounds
        currVY(bounce) = -currVY(bounce);
        currY(:,1) = currY(:,1) + currVY(:,1)*dt;
        
    end
    
    % save previous position to plot trajectory
    prevX(:,n+1) = currX;
    prevY(:,n+1) = currY;
    
    % Calculate average kinetic energy of electrons:
    avgE_k = C.m_0*(sum(sqrt(currVX.^2 + currVY.^2).^2)/numElec)/2;
    % Calculate temperature of material:
    currTemp(n+1) = (2*avgE_k)/(3*C.kb);
    
end

% make plot

figure(1)
% Colour = hsv(nPlottedElec);
for i = 1:nPlottedElec
    
    plot(prevX(i,:),prevY(i,:))
    hold on

end
hold off
axis([0 200e-9 0 100e-9])

figure(2)
plot(currTime,currTemp)
xlabel('Time (fs)')
ylabel('Temperature')


