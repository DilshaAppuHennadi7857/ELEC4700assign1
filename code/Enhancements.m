
global B

%%
%
% In this section, we will add boxes to the system which electrons will not
% be able to pass through. We can define the sides of the boxes as occuring
% at some fraction of the region length/width.
%

% Define boxes
% Box 1 (top)
B.Left1 = 0.35;
B.Right1 = 0.65;
B.Top1 = 1;
B.Bottom1 = 0.6;

% Box 2 (bottom)
B.Left2 = 0.35;
B.Right2 = 0.65;
B.Top2 = 0.5;
B.Bottom2 = 0;

currX = (nomRegionL).*rand(numElec, 1); % set random initial x position
currY = (nomRegionW).*rand(numElec, 1); % set random initial y position

%%
%
% Before moving forward, we must ensure all initial positions lie outisde
% the boxes.
%
for i = 1:numElec
    while ((currX(i)>=(B.Left1*nomRegionL)) && (currX(i)<=(B.Right1*nomRegionL)) && ((currY(i)>=(B.Bottom1*nomRegionW)) || (currY(i)<=(B.Top2*nomRegionW))))
        currX = (nomRegionL).*rand(numElec, 1); % set random initial x position
        currY = (nomRegionW).*rand(numElec, 1); % set random initial y position
    end
end

currVel = zeros(numElec,1); % initialize column vector to hold velocity, velocities initialized with zero
currDir = (2*pi).*rand(numElec,1); % create column vector with current direction of each electron
currVX = []; % set up x and y velocity vectors as empity column vector
currVY = [];
prevX = []; % clear vars
prevY = [];
currTime = [];
randVal = [];
bounce = [];

% Calculate velocity vectors
currVX = v_th + sqrt(C.kb*Temp/C.m_0)*randn(numElec,1);
currVY = v_th + sqrt(C.kb*Temp/C.m_0)*randn(numElec,1);

%%
%
% The probability of scattering is given as: P_scat = 1 - exp(-dt/t_mn)

P_scat = 1 - exp(-(dt)/t_mn);

for n = 0:numTimeStep
    
    currTime(n+1) = n*dt; % determine current time (ms)
    
    % Update according to Newton's laws of motion
    % new position = old position + velocity*time
    
    if n > 0 %update position after t=0
        
        randVal = rand(numElec,1); % assign scatter probability
        scatter = randVal<=P_scat;
        currVX(scatter) = v_th + sqrt(C.kb*Temp/C.m_0)*randn;
        currVY(scatter) = v_th + sqrt(C.kb*Temp/C.m_0)*randn;
        
        % chance to invert direction when scattering
        randVal = rand(numElec,1);
        invertDir = scatter & (randVal<=0.5);
        currVX(invertDir) = -currVX(invertDir);
        currVY(invertDir) = -currVY(invertDir);
        
        checkBounce;
        
        % if e- crosses right bound, set X=0; if e- crosses left bound, set
        % X=upper bound
        crossRight = newX > nomRegionL;
        currX(crossRight) = 0;
        crossLeft = newX < 0;
        currX(crossLeft) = nomRegionL;
        
    end
    
    % save previous position to plot trajectory
    prevX(:,n+1) = currX;
    prevY(:,n+1) = currY;
    
    currVel = sqrt(currVX.^2 + currVY.^2);
    
    % Calculate average kinetic energy of electrons:
    avgE_k = C.m_0*(sum(currVel.^2)/numElec)/2;
    % Calculate temperature of material:
    currTemp(n+1) = (2*avgE_k)/(3*C.kb);
    
end

% make plot

figure(6)
for i = 1:nPlottedElec
    plot(prevX(i,:),prevY(i,:))
    hold on
end
makeBox(B.Left1,B.Right1,B.Top1,B.Bottom1,nomRegionL,nomRegionW);
makeBox(B.Left2,B.Right2,B.Top2,B.Bottom2,nomRegionL,nomRegionW);
hold off


