%% One step expansion of Dijkstra's algorithm
function [u, Len, H, P, Dis_s] = DijkstraExpansion(rp, ci, vi, H, P, Dis_s, Len)

%%% Input parameters
%
% rp, ci, ai: See sparse_to_csr.m
% H: The heap of node indices
% P: The heap positions of nodes
% Dis_s: The distance vector of each node to the source node
% Len: The length of the heap
%
% This is a modification of the codes of David F. Gleich, Stanford University.

%% Pop the head off the heap
u = H(1);
Tail = H(Len);
H(1) = Tail;
P(Tail) = 1;
Len = Len-1;

%% Maintain the min-heap
[H, P] = MinHeap(H, P, Dis_s, Len);

%% Relax the neighbors of u
[Len, H, P, Dis_s] = Relax(rp, ci, vi, H, P, Dis_s, u, Len);

end

%% Maintaining min-heap
function [H, P] = MinHeap(H, P, Dis_s, Len)

Pos = 1;
FirstVertex = H(Pos);

% Move the first node down the heap
while 1
    
    Idx = 2*Pos;
    
    if Idx > Len % No child
        break;
    elseif Idx == Len % One child
        v = H(Idx);
    else % Two children, pick the smaller one
        Left = H(Idx);
        Right = H(Idx+1);
        v = Left;
        if Dis_s(Right) < Dis_s(Left)
            Idx = Idx+1;
            v = Right;
        end
    end
    
    if Dis_s(FirstVertex) < Dis_s(v)
        break;
    else
        H(Pos) = v;
        P(v) = Pos;
        H(Idx) = FirstVertex;
        P(FirstVertex) = Idx;
        Pos = Idx;
    end
    
end

end

%% Neighbor relax
function [Len, H, P, Dis_s] = Relax(rp, ci, vi, H, P, Dis_s, u, Len)

for EdgeIdx = rp(u):rp(u+1)-1
    
    v = ci(EdgeIdx); % v is a neighbor of u
    EdgeWeight = vi(EdgeIdx);
    
    % Relax edge (u, v)
    if Dis_s(v) > Dis_s(u)+EdgeWeight
        
        Dis_s(v) = Dis_s(u)+EdgeWeight;
        Pos = P(v);
        
        if Pos == 0 % v is not in the heap
            Len = Len+1;
            H(Len) = v;
            P(v) = Len;
            Pos = Len;
        end
        
        % Move v up the heap
        while Pos > 1
            
            ParPos = floor(Pos/2);
            Par = H(ParPos);
            
            if Dis_s(Par) < Dis_s(v)
                break;
            else
                H(ParPos) = v;
                P(v) = ParPos;
                H(Pos) = Par;
                P(Par) = Pos;
                Pos = ParPos;
            end
            
        end
        
    end
    
end

end