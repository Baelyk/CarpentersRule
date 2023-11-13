% Random Chord-arc Constant Calculation
function chordarc = find_chordarc(polygon, samples)
	n = length(polygon) - 1;
	reorder = polygon;

	a=samples; % (input) number of random points on each edge
	ran = linspace(0, 1, samples)';
	% a a*1 matrix of random scalar
	% drawn from the uniform distribution in the interval (0,1).


	%for k = 1:a-1
		%for l = 2:a
			%if ran(k)== ran(l) % update a to get unique values
				%ran(l) = rand(1);
			%else
				%break
			%end
		%end
	%end

	%% reorder ran
	%ran = sort(ran);

	edgePoints = [];
	for idx = 1:n
		P1=reorder(idx,:);  % every edge has some random points selection
		P2=reorder(idx+1,:);
		ep = bsxfun(@times,ran,P2-P1)+P1; %select random convex comb
		edgePoints = [edgePoints;P1;ep];
		% let v denote vertex, e denote points on the edges between 2 vertices
		% this follows:
		% v1 e_11...e_1a v2 e_21...e_2a v3...vn e_n1...e_na
	end
	edgePoints = [edgePoints;reorder(end,:)];
	% this connects the last point on edge to the first vertex
	%plot(edgePoints(:,1),edgePoints(:,2),'ro-','LineWidth',2)





	% Chord-Arc Constant
		% shortestpath()

	vertex_size = size(edgePoints(:,1),1)-1;
	s =  [];
	t = [];
	for i = 1:(vertex_size-1)
		s = [s;i]; % this connects every points to
		t = [t;i+1];% the next one (last one is the first vertex)
	end
	s = [s;vertex_size];
	t=[t;1];
	% need to calculate weight for every edge (1-2, 2-3, ... , n-1)
	weights = [];
	for i = 1:vertex_size
		dis = norm(edgePoints(i,:)-edgePoints(i+1,:));
		weights = [weights;dis];
	end

	G = graph(s, t, weights);
	%plot(G,'EdgeLabel',G.Edges.Weight)
	%trying to output the result with weights but does not work

	% Caculate Chord-Arc Constant between any two points
	chordarc = [];
	for k = 1:vertex_size-1
		for l = k+1:vertex_size
			[P,m] = shortestpath(G,k,l); % should return the shortest path
			% p will return the path
			% m will return the length
			euclid = norm (edgePoints(k,:)-edgePoints(l,:));
			c = m/euclid;
			chordarc = [chordarc;c];
		end
	end


	% To-do
	%find supremum
	chordarc = max(chordarc);
end

% Vertex chord arc:
%function chord_arc_const_list = chordarc(positions, lengths, count)
	%chord_arc_const_list=[]
	%L = lengths;

	%for j = 1 : count
		%P = positions{j};

		%% add the new chord arc constant to the list
		%chord_arc_const_list(end+1) = chord_arc_constant(P,L);
	%end

	%% plot constant vs iteration
	%% clf;
	%figure
	%plot(chord_arc_const_list)
	%xlabel('nth iteration')
	%ylabel('chord arc constant')
	%print("chord_arc_constant" + datestr(datetime, 31), "-dpng")
%end

%function C = chord_arc_constant(P, L)
    %C=0;
    %for i = 1 : size(P, 1)
        %for j = i+1 : size(P, 1)
            %arc_length = sum(L(i:j-1));
            %euc_dis = norm(P(j,:)-P(i,:));
            %const = arc_length/euc_dis;
            %if const>C
                %C=const;
            %end
        %end
    %end
%end
