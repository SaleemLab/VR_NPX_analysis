function dist = wheel2unit(wheelTicks,ticksPerRev, diameter)
    circum = pi*diameter;
    dist = (wheelTicks/ticksPerRev) * circum; % distance in wheel rotation, * the cms of that
    
    