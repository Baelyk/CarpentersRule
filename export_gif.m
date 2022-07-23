function export_gif(filename, positions, count)
	fps = 15;
	padding = 0.25;

	% Get max and min coordinates
	max_coord = 0;
	for j = 1 : count
		max_coord = max([max_coord; max(positions{j}(:))]);
	end
	min_coord = 0;
	for j = 1 : count
		min_coord = min([min_coord; min(positions{j}(:))]);
	end

	fig = figure("Position", [0 0 280 210], "Resize", false);
	for j = 1 : count
		fprintf("\r%4u/%4u", j, count);

		% Plot the jth P
		P = positions{j};
		plot(P(:, 1), P(:, 2), "-o");
		title("i = " + num2str(j));
		% Set the axis
		axis([min_coord max_coord min_coord max_coord] + padding * [-1 1 -1 1]);

		% Capture the figure as a frame and converting it to the image data
		% `imwrite` needs
		frame = getframe(fig);
		[A, map] = rgb2ind(frame2im(frame), 256);

		% Add this frame to the gif, appending after the first frame
		if j == 1
			imwrite(A, map, filename, "gif", "LoopCount", inf, "DelayTime", 1 / fps);
		else
			imwrite(A, map, filename, "gif", "WriteMode", "append", "DelayTime", 1 / fps);
		end
	end
	fprintf("\n");
end
