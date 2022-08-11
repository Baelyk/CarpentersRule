function [polygons, speeds, chordarcs] = spiral_speedlimit(max_epsilon, min_epsilon, num, do_plot)
	if ~exist("do_plot", "var")
		do_plot = false;
	end

	polygons = cell(1, num);
	speeds = zeros(num, 1);
	chordarcs = zeros(num, 1);
	ignored = false(num, 1);

	resolution = 10^-3;

	epsilons= linspace(max_epsilon, min_epsilon, num);

	for idx = 1 : num
		epsilon = epsilons(idx);
		fprintf("\r%4u/%4u, epsilon = %6.4g", idx, num, epsilon);
		polygon = spiral(epsilon);

		[vel, exitflag] = find_velocity(polygon, false, resolution);

		if exitflag ~= 1
			fprintf(" !! find_velocity returned %d, will not plot\n", exitflag);
			ignored(idx) = true;
		end

		speeds(idx) = max(norm(vel));
		%chordarcs(idx) = find_chordarc(polygon, 10);
	end

	if do_plot
		plot_speeds = speeds(~ignored);
		plot_idx = 1 : num;
		plot_idx = plot_idx(~ignored);
		plot_epsilons = epsilons(~ignored);
		semilogy(plot_epsilons, plot_speeds, "o");
	end

	output = [speeds chordarcs];
	timestamp = datestr(datetime, 30)
	writematrix(output, "speeds_" + timestamp + ".csv");
	writecell(polygons, "polygons_" + timestamp + ".csv");
end
