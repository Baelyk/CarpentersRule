function heatmap(positions, count, v_res, filename)
    colors = colormap(winter);

	skip = 2;
	if count < length(positions)
		skip = 1 + round(length(positions) / count);
		count = length(positions);
	end
	% Skip by skip, but also skip the first since it's plotted manually, end
	% with the last
	js = [1 : skip : length(positions), length(positions)];
	js = js(2 : end);

	% Get max and min coordinates
	max_coords = [-Inf -Inf];
	for j = 1 : count
		max_coords = max([max_coords; max(positions{j})]);
	end
	min_coords = [Inf Inf];
	for j = 1 : count
		min_coords = min([min_coords; min(positions{j})]);
    end
    % Create figure with proper aspect ratio
    aspect_ratio = abs((max_coords(1) - min_coords(1)) ...
        / (max_coords(2) - min_coords(2)));
    h_res = round(v_res * aspect_ratio);

	figure("Position", [0 0 h_res v_res], "Resize", false);
    hold on;
	p = positions{1};
	p = plot(p(:, 1), p(:, 2), "-o", LineWidth = 4);
	p.Color = "#c5050c";
    for j = js
        % Get ith position
        p = positions{j};
        % Shadow p with its plot
        p = plot(p(:, 1), p(:, 2), "-o", LineWidth = 3);
        % Change the color of p to the next color in the heatmap
        color = colors(ceil(j / count * length(colors)), :);
        p.Color = color;
    end

	set(gca, Visible = false);
	exportgraphics(gca, filename + ".png");
	%print(filename, "-dpng")
end
