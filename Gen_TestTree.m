function [tree] = Gen_TestTree(sidedness,minPart,N)
%function [tree] = Gen_TestTree(sidedness,N)
%This code is for testing purposes only. Artificually generates a tree 
%structure that is not associated with any function.  Function takes in the
%size of the desired matrix and generates a tree by splitting this
%repeatedly in half on the left and right for a 'complete' tree.
%left/right trees are created by cutting asymmetricly and having the 
%upper/lower diagonal block contain only 1 entry. A symmetric tree will have
%log2(N) cuts, and a right/left tree will have N cuts. A tree skewed to the
%center ('center') will have N/2 cuts.
%
%
%INPUT:             minPart     (int) minimum number of partitions 
%                   sidedness   {'left', 'right', 'center', 'complete'} 
%                               type of tree that will be generated
%                   N           (scalar) number of points in desired 
%                               function
%OUTPUT:            tree        (structure) contains split dimensions at   
%                               each level, as well as whether or not the 
%                               node is a leaf
%
% Author: Kristen Lessel, klessel@engineering.ucsb.edu
%
% Copyright (C) 2015 Kristen Lessel, (klessel@engineering.ucsb.edu)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

tree.m = N;
tree.leaf = 0;
tree.root = 1;
if strcmp(sidedness, 'complete')
    
    [tree] = Gen_TestTree_complete_iter(tree,minPart,N);
    
elseif strcmp(sidedness, 'left')
    
    [tree] = Gen_TestTree_left_iter(tree,minPart,N);
    
elseif strcmp(sidedness, 'right')
    
    [tree] = Gen_TestTree_right_iter(tree,minPart,N);

elseif strcmp(sidedness, 'center')
    
    [tree] = Gen_TestTree_center_iter(tree,minPart,N);
    
else
    
    display(sprintf('The variable sidedness can take the value ''left'', ''right'', ''center'', or ''symmetric''. \nPlease enter one of these. Usage:\nfunction [tree] = Gen_TestTree(N,sidedness)'))
    
end

end

function [tree] = Gen_TestTree_complete_iter(tree,minPart,N)
if N/2 < minPart
    tree.m = N;
    tree.leaf = 1;
    tree.root = 0;
else
    
tree.treeL.leaf = 0;
tree.treeL.root = 0;
% If N/2 is not an integer, give the left branch one more row than the
% right.
if N/2 ~= floor(N/2)    
    tree.treeL.m = ceil(N/2);
    treeL = Gen_TestTree_complete_iter(tree.treeL,minPart,ceil(N/2));
else
    tree.treeL.m = N/2;
    treeL = Gen_TestTree_complete_iter(tree.treeL,minPart,N/2);
end
tree.treeL = treeL;


tree.treeR.leaf = 0;
tree.treeR.root = 0;
if N/2 ~= floor(N/2)
    tree.treeR.m = floor(N/2);
    treeR = Gen_TestTree_complete_iter(tree.treeR,minPart,floor(N/2));    
else
    tree.treeR.m = N/2;
    treeR = Gen_TestTree_complete_iter(tree.treeR,minPart,N/2);
end
tree.treeR = treeR;


end
end

function [tree] = Gen_TestTree_left_iter(tree,minPart,N)
if N/2 < minPart
    tree.m = N;
    tree.leaf = 1;
    tree.root = 0;
else
    
tree.treeL.m = N-minPart;
tree.treeL.leaf = 0;
tree.treeL.root = 0;

tree.treeR.m = minPart;
tree.treeR.leaf = 1;
tree.treeR.root = 0;

treeL = Gen_TestTree_left_iter(tree.treeL,minPart,N-minPart);
tree.treeL = treeL;

end
end

function [tree] = Gen_TestTree_right_iter(tree,minPart,N)
if N/2 < minPart
    tree.m = N;
    tree.leaf = 1;
    tree.root = 0;
else
    
tree.treeR.m = N-minPart;
tree.treeR.leaf = 0;
tree.treeR.root = 0;

tree.treeL.m = minPart;
tree.treeL.leaf = 1;
tree.treeL.root = 0;

treeR = Gen_TestTree_right_iter(tree.treeR,minPart,N-minPart);
tree.treeR = treeR;

end
end

function [tree] = Gen_TestTree_center_iter(tree,minPart,N)
    
    tree.treeL.m = N/2;
    tree.treeL.leaf = 0;
    tree.treeL.root = 0;
    treeL = Gen_TestTree_right_iter(tree.treeL,minPart,N/2);
    tree.treeL = treeL;
    
    tree.treeR.m = N/2;
    tree.treeR.leaf = 0;
    tree.treeR.root = 0;
    treeR = Gen_TestTree_left_iter(tree.treeR,minPart,N/2);
    tree.treeR = treeR;

end
