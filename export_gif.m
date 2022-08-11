function export_gif(filename, positions, count)
	min_fps = 15;
	max_s = 10;
	max_fps = 100;
	padding = 0.25;
    h_res = 200;

	skip = 2;
	fps = max(min_fps, count / max_s);
	if fps > max_fps
		fps = max_fps;
		skip = 1 + ceil(count / fps);
	end
	fprintf("%u fps, skipping by %u\n", fps, skip - 1);
	% Start at one, the increase to skip (by skip - 1), but do include the last
	% frame
	js = [1 : skip : count, count];

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
    v_res = h_res * aspect_ratio;

	fig = figure("Position", [0 0 v_res h_res], "Resize", false);
	for j = js
		fprintf("\r%4u/%4u", j, count);

		% Plot the jth P
		P = positions{j};
		plot(P(:, 1), P(:, 2), "-o");
		title("i = " + num2str(j));
		% Set the axis
		axis([min_coords(1) max_coords(1) min_coords(2) max_coords(2)] ...
            + padding * [-1 1 -1 1]);

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
