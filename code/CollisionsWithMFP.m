%%
%
% We reset the position and velocity variables from the previous section
% for use here.
%
currX = (nomRegionL).*rand(numElec, 1); % set random initial x position
currY = (nomRegionW).*rand(numElec, 1); % set random initial y position
currVel = zeros(numElec,1); % initialize column vector to hold velocity, velocities initialized with zero
currDir = (2*pi).*rand(numElec,1); % create column vector with current direction of each electron
currVX = []; % set up x and y velocity vectors as empity column vector
currVY = [];
prevX = []; % clear vars
prevY = [];
currTime = [];
randVal = [];
bounce = [];
numScat = [];

%%
% To calculate velocity vectors, the vector components are randomly
% assigned according to the Boltzmann-Maxwell Distribution. For velocity
% components, this distribution is Gaussian with a standard deviation of
% sqrt(k*T/m). The distribution is also adjusted such that the average
% velocity is the thermal velocity.

currVX = v_th + sqrt(C.kb*Temp/C.m_0)*randn(numElec,1);
currVY = v_th + sqrt(C.kb*Temp/C.m_0)*randn(numElec,1);

%%
%
% We check the distribution of the velocities by obtaining the magnitude of
% the velocities and plotting them into a histogram.
%

currVel = sqrt(currVX.^2 + currVY.^2);

figure(3)
hist(currVel);

%%
%
% In this section, we will be studying the effects of scattering. The
% probability of scattering is given by P_scat = 1 - exp(-dt/t_mn). We can
% use this formula to obtain the probability of scattering within our
% system. We can also obtain a rough idea of the scattering of electrons by
% observing a random electron. In this case, we will observe the scattering
% of electron 1.
%

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
        
        % if electron 1 scatters, count it
        if scatter(1,1)
            numScat(n+1) = numScat(n+1) + 1;
        end
        
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
    
    currVel = sqrt(currVX.^2 + currVY.^2);
    
    % Calculate average kinetic energy of electrons:
    avgE_k = C.m_0*(sum(currVel.^2)/numElec)/2;
    % Calculate temperature of material:
    currTemp(n+1) = (2*avgE_k)/(3*C.kb);
    
end

% make plot

figure(4)
for i = 1:nPlottedElec
    plot(prevX(i,:),prevY(i,:))
    hold on
end
hold off

%%
%
% Throughout the iterations, the temperature of the system was taken. With
% the temperature at each time step recorded, the data can be plotted to
% observe the temperature of the system over time.
%

figure(5)
plot(currTime,Temp,'o-')
xlabel('Time (fs)')
ylabel('Temperature')

%%
%
% To calculate the mean time between collisions, we identify the times a
% collision occurs for an electron. Then, we can look at the difference in
% time between scattering events and take the average. Then obtain the mean
% free path as the product of the electron's velocity and mean time between
% collisions.
%

deltaT = 0; % initialized time
sumT = 0;

for t = 0:(numTimeStep-1)
    
    if (numScat(t+1) > 0)
        
        sumT = sumT + deltaT; % every time a scattering event occured take record of how much time had passed
        deltaT = 0; % reset "timer"
        
    end
    
    deltaT = deltaT + dt;
    
end

%%
% 
% calculate mean time between collisions by dividing the sum of time
% between collisions by the number of collision events.
%

calc_t_mn = sumT/sum(numScat);

calc_L_n = v_th*calc_t_mn;

