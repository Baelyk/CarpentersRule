function [output, polygons, ignored] = speedlimit(options)
	arguments
		% The shape to generate
		options.Shape (1, 1) string { mustBeMember(options.Shape, ["random", "spiral", "ngonarc", "w"]) }
		% Number of vertices for random polygons, [start_vertices, end_vertices]
		% range for ngonarcs
		options.Vertices (1, :) double = 10
		% The number of shapes, except for Shapes ngonarc and w, which override
		% this with their arguments
		options.Num (1, 1) double = 10
		% The [start_epsilon, end_epsilon] range for the spirals
		options.Epsilon (1, 2) double = [1/4, 1/8]
		% The [start_num_us, end_num_us] range for the npleus
		options.Us (1, 2) double = [1, 10]

		% The expander resolution
		options.Resolution (1, 1) double = 10^-3
		% Number of samples for the chord arc calculator
		options.ChordArcSamples (1, 1) double = 10

		% Whether to save the output
		options.SaveOutput (1, 1) logical = true
		% Whether or not to plot the results
		options.Plot (1, 1) logical = false
		% Whether or not to plot against the shape generate domain (e.g.
		% epsilon, vertices) instead of chordarc
		options.PlotAgainstDomain (1, 1) logical = false
	end

	is_polygon = false;

	switch options.Shape
		case "random"
			generate_shape = @(num_vertices) generate_polygon(num_vertices);
			domain = options.Vertices(1) * ones(options.Num, 1);
			is_polygon = true;
		case "spiral"
			generate_shape = @(epsilon) spiral(epsilon);
			domain = linspace(options.Epsilon(1), options.Epsilon(2), options.Num);
		case "ngonarc"
			if size(options.Vertices(:), 1) ~= 2
				error("Vertices must be [min_vertices, max_vertices] with shape ngonarc");
			end
			generate_shape = @(num_vertices) ngonarc(num_vertices);
			domain = options.Vertices(1) : options.Vertices(2);
		case "w"
			generate_shape = @(num_us) npleu(num_us);
			domain = options.Us(1) : options.Us(2);
		otherwise
			error("Unexpected shape %s", options.Shape);
	end

	num = length(domain);

	polygons = cell(1, num);
	speeds = zeros(num, 1);
	chordarcs = zeros(num, 1);
	ignored = false(num, 1);

	for idx = 1 : num
		% Generate the shape and log this iteration
		x = domain(idx);
		fprintf("\r%4u/%4u, x = %6.4g ... finding velocity", idx, num, x);
		polygon = generate_shape(x);
		polygons{idx} = polygon;

		% Find the velocity
		[vel, exitflag] = find_velocity(polygon, is_polygon, options.Resolution);

		% Ignore this velocity if something went iffy
		if exitflag ~= 1
			fprintf(" !! find_velocity returned %d, will not plot\n", exitflag);
			ignored(idx) = true;
		end

		% Save the max speed (vecnorm(vel') does norm of each row)
		speeds(idx) = max(vecnorm(vel'));

		% Calculate the chordarc, only if not plotting against domain
		if ~options.PlotAgainstDomain
			fprintf("\r%4u/%4u, x = %6.4g ... calculating chord arc with %d samples", ...
				idx, num, x, options.ChordArcSamples);
			chordarcs(idx) = find_chordarc(polygon, options.ChordArcSamples);
		end
	end

	% Plot the results
	if options.Plot
		plot_speeds = speeds(~ignored);
		if options.PlotAgainstDomain
			plot_X = domain(~ignored);
		else
			plot_X = chordarcs(~ignored);
		end
		semilogy(plot_X, plot_speeds, "o");
	end

	output = [chordarcs speeds];
	% Save the results
	if options.SaveOutput
		timestamp = datestr(datetime, 30);
		if options.Plot
			print("plot_" + timestamp, "-dpng")
		end
		writematrix(output, "speeds_" + timestamp + ".csv");

		filename = "polygons_" + timestamp + ".csv";
		for idx = 1 : num
			polygon = polygons{idx};
			writematrix(polygon, filename, WriteMode = "append");
			writelines("", filename, WriteMode = "append");
		end
	end

	fprintf("\n");
end

function polygon = generate_polygon(num_vertices)
	polygon = polygon_fxn(num_vertices);
	polygon = polygon(1 : end - 1, :);
end
