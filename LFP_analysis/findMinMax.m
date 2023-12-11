function findMinMax(w, getY, markMaxPoint, markMinPoint)
    points = zeros(1, w-5);  % Since MATLAB indexing starts from 1

    minY = 10000000;
    maxY = 0;
    
    for ix = 3:w-2
        height = getY(ix);
        if isempty(height) 
            points(ix-2) = -1;
            continue;
        end
        points(ix-2) = height;

        if maxY < height
            maxY = height;
        end
        if minY > height
            minY = height;
        end
    end

    lookingFor = -1;
    localMaxY = -inf;  % To ensure that any real number will be higher
    localMinY = inf;   % To ensure that any real number will be lower
    localMaxYx = 0;
    localMinYx = 0;

    minChange = (maxY - minY) * 0.1;
    nextMaxY = points(1) + minChange;
    nextMinY = points(1) - minChange;
    npoints = length(points);
    for x = 1:npoints
        y = points(x);

        if y == -1
            continue;
        end

        if y > nextMaxY
            if lookingFor == 1
                markMinPoint(localMinYx, localMinY);
            end
            nextMinY = y - minChange;
            nextMaxY = y + minChange;
            lookingFor = 0; % look for minimum, but save highest until we find it
            localMaxY = y;  % reset the local highest
            localMaxYx = x;
        end
        if localMaxY <= y
            localMaxY  = y;
            localMaxYx = x;
        end

        if y < nextMinY
            if lookingFor == 0
                markMaxPoint(localMaxYx, localMaxY);
            end
            nextMaxY = y + minChange;
            nextMinY = y - minChange;
            lookingFor = 1; % look for maximum, but save lowest until we find it
            localMinY = y;  % reset the local lowest
            localMinYx = x;
        end
        if localMinY >= y
            localMinY = y;
            localMinYx = x;
        end
    end
end
