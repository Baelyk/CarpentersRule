clear
clf
close all

% n -- number of vertices 
n=10
points = 0.1*randi([0 10],n,2) % choosing random points

for i =1:n-1
    for j = 2:n
        if points(i,:) == points(j,:) % update points if two are the same
            points(j,:) = 0.1*randi([0 10],1,2)
        else
            break
        end
    end
end


% reordering vertices to avoid self-intersection
xc = points(:,1);
yc = points(:,2);
meanx = mean(xc)
meany = mean(yc)
angles = atan2( (yc-meany) , (xc-meanx))
[sortedAngles, sortIndices] = sort(angles);
reorderedx = xc(sortIndices)
reorderedy = yc(sortIndices)

reorderedx(end+1) = reorderedx(1)
reorderedy(end+1) = reorderedy(1)


% plotting grid + points
x = linspace(0,10);
y = linspace(0,10);
plot(reorderedx, reorderedy, 'bo-', 'LineWidth', 2);
grid on;
xlim([0 1])
ylim([0 1])

