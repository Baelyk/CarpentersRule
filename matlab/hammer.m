function points = hammer(epsilon, handle_length)
	points = [...
		0, 0
		1, 0
		1, 1 - epsilon
		handle_length, 1 - epsilon
		handle_length, 1 + epsilon
		1, 1 + epsilon
		1, 2
		0, 2
		0, 0 ...
	];
end
