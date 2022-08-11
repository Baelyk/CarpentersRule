function heatmap(positions, count)
    colors = colormap;

	figure;
    hold on;
    for j = 1 : count
        % Get ith position
        p = positions{j};
        % Shadow p with its plot
        p = plot(p(:, 1), p(:, 2), "-o");
        % Change the color of p to the next color in the heatmap
        color = colors(ceil(j / count * length(colors)), :);
        p.Color = color;
    end
end
