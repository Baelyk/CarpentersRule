function x = poly_gamma(polygon, t)
	% Get the lengths of each segment of the polygon
	n = size(polygon, 1);
	lengths = vecnorm((polygon(1 : n - 1, :) - polygon(2 : n, :))');
	distance = cumsum(lengths);
	edge_number = find(distance >= t, 1);
	lambda = (distance(edge_number) - t) / distance(edge_number);
	x1 =
end
