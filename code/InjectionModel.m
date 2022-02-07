
global B

%%
%
% Continuing on from the boxes, we'll now simulate an injection of
% electrons into the system. The electrons will start out on the left side
% of the region at the mid-way point of the width and will be sent in with
% a positive x velocity and 0 y velocity component. Right away the
% electrons may experience scattering.
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
B.Top2 = 0.45;
B.Bottom2 = 0;

currX = zeros(numElec,1); % start electrons a left side of the screen
currY = 0.5*nomRegionW.*ones(numElec, 1); % middle of the region width
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
currVY = zeros(numElec,1); % no Y component to velocity vector initially

%%
%
% The probability of scattering is given as: P_scat = 1 - exp(-dt/t_mn)

P_scat = 1 - exp(-(dt)/t_mn);
numScat = zeros(1,numTimeStep+1); % just look at scattering of electron 1

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
        
    end
    
    % save previous position to plot trajectory
    prevX(:,n+1) = currX;
    prevY(:,n+1) = currY;
    
    % make plot

    figure(7)
    for i = 1:nPlottedElec
        plot(prevX(i,:),prevY(i,:))
        hold on
    end
    axis([0 nomRegionL 0 nomRegionW])
    makeBox(B.Left1,B.Right1,B.Top1,B.Bottom1,nomRegionL,nomRegionW);
    makeBox(B.Left2,B.Right2,B.Top2,B.Bottom2,nomRegionL,nomRegionW);
    hold off
    
end


