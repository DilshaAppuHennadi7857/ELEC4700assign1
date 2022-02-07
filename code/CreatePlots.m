t = 0;              % start time
dt = 0.1e-15;         % time step
TStop = 1000*dt;    % stop time

%%
%
% To create all plots for each part in the same while, a different set of
% velocity and position vectors will be made for each plot. They will all
% be initialized with the same value but as electrons are updated with each
% time step, the velocities and positions may be different.

% Part 1
Vx1 = Vx;
Vy1 = Vy;
x1 = x;
y1 = y;

% Part 2
% Vx2 = Vx;
% Vy2 = Vy;
% x2 = x;
% y2 = y;

% Part 3
% Vx3 = Vx;
% Vy3 = Vy;
% x3 = x;
% y3 = y;

while t < TStop
    % Part 1 Plot
    ElectronModelingPlot
    
    % Part 2 Plot
    
    % Part 3 Plot
    
    t = t + dt; % increase time
    
    pause(0.00001)
end