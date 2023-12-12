function wheel = import_processWheelEphys(wheel,smoothMethod,windowSize)

% Set default values
ticksPerRev = 1024; % N
diameter = 17;      % cm


halfMax = max(wheel.pos)/2;
wheel.unwrapped = unwrap(wheel.pos, halfMax);

% not sure I have the right ticks per rev.
wheel.dist = wheel2unitInt(wheel.unwrapped, ticksPerRev, diameter).*-1; % pos, ticks/rev, wheel diam

wheel.rawSpeed = diff(wheel.dist)./diff(wheel.Time);
wheel.rawSpeed = movmean(wheel.rawSpeed, 2);
wheel.rawSpeed = [wheel.rawSpeed(1); wheel.rawSpeed];
wheel.smthSpeed = smoothdata(wheel.rawSpeed,smoothMethod,windowSize);

wheel.acc = diff(wheel.rawSpeed)./diff(wheel.Time);
wheel.acc = movmean(wheel.acc,2);
wheel.acc = [wheel.acc(1); wheel.acc];
wheel.smthAcc = smoothdata(wheel.acc,smoothMethod,windowSize);

wheel.timeStat = numel(find(abs(wheel.rawSpeed)<1)) * median(diff(wheel.Time));
wheel.timeMove = numel(find(abs(wheel.rawSpeed)>1)) * median(diff(wheel.Time));

moveIdx = find(abs(wheel.rawSpeed)>1);
wheel.meanRunSpeed = mean(wheel.rawSpeed(moveIdx));
wheel.stdRunSpeed = std(wheel.rawSpeed(moveIdx));

end

function dist = wheel2unitInt(wheelTicks,ticksPerRev, diameter)
    circum = pi*diameter;
    dist = (wheelTicks/ticksPerRev) * circum; % distance in wheel rotation, * the cms of that
end    