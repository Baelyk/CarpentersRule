%reorder = [1 0; 0 0 ; 0 1; 1, 0]
function [result] = chordarcfxn(reorder,a)

    ran = linspace(0,1,a)';
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


    reorder_size = size(reorder(:,1),1)-1
    edgePoints = [];
    for idx = 1:reorder_size
        P1=reorder(idx,:);  % every edge has some random points selection
        P2=reorder(idx+1,:);
        ep = bsxfun(@times,ran,P2-P1)+P1; %select random convex comb
        edgePoints = [edgePoints;P1;ep];
        % let v denote vertex, e denote points on the edges between 2 vertices
        % this follows:

        % v1 e_11...e_1a v2 e_21...e_2a v3...vn e_n1...e_na
    end
    edgePoints = [edgePoints;reorder(end,:)];

    plot(edgePoints(:,1),edgePoints(:,2),'-o')
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

    % Caculate Chord-Arc Constant between any two points
    chordarc = [];
    cor1 = [];
    cor2 = [];
    result = [];
    for k = 1:vertex_size-1
        for l = k+1:vertex_size
            [P,m] = shortestpath(G,k,l); % should return the shortest path
            % p will return the path
            % m will return the length
            euclid = norm (edgePoints(k,:)-edgePoints(l,:));
            c = m/euclid;
            chordarc = [chordarc;c];
            cor1 = [cor1;edgePoints(k,:)];
            cor2 = [cor2;edgePoints(l,:)];
        end
    end
    result = [result; cor1 cor2 chordarc];

    mid= sortrows(result,5,"ascend")
    result = mid(end,:)

    plot(reorder(:,1),reorder(:,2),'bo-','LineWidth',2)
    hold on;
    plot(result(1),result(2),'-ro')
    plot(result(3),result(4),'-ro')

end
