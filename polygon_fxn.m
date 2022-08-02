function P = polygon_fxn(n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    points = randi([0, 10], n,2); % choosing random points

    for  i = 1:n-1
        for j = i+1:n
            while points(i,:) == points(j,:) % update points if two are the same
                points(j,:) = randi([0, 10], 1,2);
            end
        end
    end

    % reordering vertices to avoid self-intersection
    xc = points(:,1);
    yc = points(:,2);
    meanx = mean(xc);
    meany = mean(yc);
    angles = atan2( (yc-meany) , (xc-meanx));
    [sortedAngles, sortIndices] = sort(angles);
    reorderedx = xc(sortIndices);
    reorderedy = yc(sortIndices);

    reorderedx(end+1) = reorderedx(1);
    reorderedy(end+1) = reorderedy(1);

    reorder =  [reorderedx reorderedy];

    P = reorder;
end

function draw(P)
% plotting grid + points
    plot(P(:,1), P(:,2), 'bo-', 'LineWidth', 2)
    grid on;
    xlim([0 1]);
    ylim([0 1]);
end
