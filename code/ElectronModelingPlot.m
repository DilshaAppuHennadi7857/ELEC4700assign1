
MarkerSize = 12e-9;
Limits = [0 200e-9 -10e-9 110e-9];

hold on
plot(x1, y1);
hold off
axis(Limits);
title('Motion of Electrons')
xlabel('X')
ylabel('Y')

% update velocity + position
x1Prev = x1;
y1Prev = y1;

for i = 1:numAtoms
    if ((x1Prev(i) + Vx1(i)*dt) > nomRegionL)
        x1(i) = 0;
    elseif ((x1Prev(i) + Vx1(i)*dt) < 0)
        x1(i) = 200;
    else
        x1(i) = x1Prev(i) + Vx1(i)*dt;
    end

    if ((y1Prev(i) + Vy1(i)*dt) > nomRegionW)
        Vy1(i) = -Vy1(i);
        y1(i) = y1Prev(i) + Vy1(i)*dt;
    elseif ((y1Prev(i) + Vy1(i)*dt) < 0)
        Vy1(i) = -Vy1(i);
        y1(i) = y1Prev(i) + Vy1(i)*dt;
    else
        y1(i) = y1Prev(i) + Vy1(i)*dt;
    end
end
    

