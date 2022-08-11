function points = npleu(n)
	% Change n from number of Us to number of points
	n = 2 * n + 1;
	points = linspace(0, 1, n)';
	% Alternate between 1 and 0, by seeing if -1^k < 0
	points = [points, (-1) .^ (1 : n)' < 0];
end
