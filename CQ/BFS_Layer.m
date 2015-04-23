%% One layer BFS search
function [BFS_Dis, BFS_Sq, BFS_Sqh, BFS_Sqt] = BFS_Layer(G, BFS_Dis, BFS_Sq, BFS_Sqh, BFS_Sqt, Layer)

% BFS Compute breadth first search distances, times, and tree for a graph
%
% [d dt pred] = bfs(A,u) returns the distance (d) and the discover time
% (dt) for each vertex in the graph in a breadth first search 
% starting from vertex u.
%   d = dt(i) = -1 if vertex i is not reachable from u
% pred is the predecessor array.  pred(i) = 0 if vertex (i)  
% is in a component not reachable from u and i != u.
%
% [...] = bfs(A,u,v) stops the bfs when it hits the vertex v
%
% Example:
%   load_gaimc_graph('bfs_example.mat') % use the dfs example from Boost
%   d = bfs(A,1)
%
% David F. Gleich
% Copyright, Stanford University, 2008-20098
%
% This is a modification

%% Initialization
rp = G.rp;
ci = G.ci;

%% Start BFS one layer search
while BFS_Sqt-BFS_Sqh > 0 && BFS_Dis(BFS_Sq(BFS_Sqh+1)) == Layer-1
    
    BFS_Sqh = BFS_Sqh+1;
    PoPNode = BFS_Sq(BFS_Sqh); % Pop a node off the head of the queue
    
    for ri = rp(PoPNode):rp(PoPNode+1)-1
        
        Neighbor = ci(ri);
        
        if BFS_Dis(Neighbor) < 0
            BFS_Sqt = BFS_Sqt+1;
            BFS_Sq(BFS_Sqt) = Neighbor;
            BFS_Dis(Neighbor) = BFS_Dis(PoPNode)+1;
        end
        
    end
    
end

end