% Check for bouncing

%%
% First obtain new Y value for all electrons. We will also need to check X
% values, we will use the ones calculated after the boundary check. So, we
% are checking newY and currX for bouncing.

newY = currY + currVY*dt;
newX(:,1) = currX(:,1) + currVX(:,1)*dt;

% bounce electron if it hits bounds

% hits roof or floor, adjust Y velocity
bounceBound = (newY>nomRegionW) | (newY<0);
currVY(bounceBound) = -currVY(bounceBound);

% hits left or right side of either box, adjust X velocity
bounceLeft = (((newY<=(B.Top2*nomRegionW)) | (newY>=(B.Bottom1*nomRegionW))) & ((newX>=(B.Left1*nomRegionL)) & (newX<=(B.Right1*nomRegionL))));
currVX(bounceLeft) = -currVX(bounceLeft);

% hits box bottom or top
bounceBox = ((newX>=(B.Left1*nomRegionL)) & (newX<=(B.Right1*nomRegionL)) & ((newY>=(B.Bottom1*nomRegionW)) | (newY<=(B.Top2*nomRegionW))));
currVY(bounceBox) = -currVY(bounceBox);

currY(:,1) = currY(:,1) + currVY(:,1)*dt;
currX(:,1) = currX(:,1) + currVX(:,1)*dt;