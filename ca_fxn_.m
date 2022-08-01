function [outputArg1] = ca_fxn(n,a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% (input) n -- number of vertices
% (input) a -- number of random points on each edge
    
    % reorder ran
    ran = sort(ran);
    reorder = polygon_fxn(n);
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
    % plot(edgePoints(:,1),edgePoints(:,2),'ro-','LineWidth',2)
    % shg
    
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
    
    %find supremum
    constant  = sort(chordarc);
    ca = max(constant);
    
    
    outputArg1 = ca;
end