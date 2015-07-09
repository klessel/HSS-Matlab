function [hss] = Gen_HSS_Memory_Efficient(tree,rank)
%This is a memory efficient, two pass algorithm that takes in a tree 
%structure for a matrix for an equation defined in the function, f, below 
%and returns its corresponding HSS structure.  Computation of U,R,V,W is 
%done via deepest first post-ordering.  B matrices are computed by 
%descending from root to child.  HSS matrices are then multiplied and added 
%from child to root in order to compute B.  (Bottom Up Routine to compute 
%B).  Can compute HSS structure for both symmetric and unsymmetric trees.
%INPUT:             rank    (int) largest allowable rank of hankel blocks -
%                           this determines the amount of compression, and
%                           should be compatible with your input tree.
%                           (if the maximum allowable rank is p, then the
%                           tree partitions should be no smaller than 3p)
%                           Corresponding singular values below this rank
%                           will be dropped.
%                   tree    (structure) contains partition dimensions for 
%                           each subdivision of the input matrix, for each
%                           of which there are a left and right child
%                           structure
%OUTPUT:            hss    (structure) contains matrices (Us, Vs, Bs Rs and
%                           Ws) that compose the hss representation of the
%                           input matrix, in addition to the partition
%                           dimensions
% Author: Kristen Lessel, klessel@engineering.ucsb.edu
% v1.0 created March 2015
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


%label tree depth and index at each node
d=1;
[tree, ~]= label_tree(tree,d);

N = tree.m;
[hss] = HSS_Basis_TranslationOp(tree,N,rank);

row_idx = 1;
col_idx = 1;
[hss] = HSS_Expansion_Coeffs_Diag(hss,row_idx,col_idx,N);

end

function [hss, rH, cH] = HSS_Basis_TranslationOp(tree,N,rank)
%Computes U's, V's, R's, W's and D's of HSS structure, and stores these
%heirarchically
%INPUT:             N       (int) grid size
%                   rank    (int) largest allowable rank of hankel blocks -
%                           this determines the amount of compression, and
%                           should be compatible with your input tree.
%                           (if the maximum allowable rank is p, then the
%                           tree partitions should be no smaller than 3p)
%                           Corresponding singular values below this rank
%                           will be dropped.
%                   tree    (structure) contains partition dimensions for
%                           each subdivision of the input matrix, for each
%                           of which there are a left and right child
%                           structure
%OUTPUT:            hss     (structure) contains matrices (Us, Vs, Bs Rs and
%                           Ws) that compose the hss representation of the
%                           input matrix, in addition to the partition
%                           dimensions
%                   rH      (matrix) row hankel block with 'current' U
%                           removed
%                   cH      (matrix) column hankel block with 'current' V'
%                           removed

if tree.leaf %if leaf node
    
    m =tree.m;
    d = tree.d;
    
    %generate indices for currrent diagonal block
    m_idx = d:(d+m)-1;
    %generate indices for current off diagonal hankel block
    butm_idx = [1:m_idx(1)-1 m_idx(end)+1:N];
    
    %function call to generate lowest level row and column hankel blocks
    rH2 = f(m_idx,butm_idx,N);

    cH2= f(butm_idx,m_idx,N);
    
    %take svds of upper row/left column and lower row/right column hankel blocks
    [Ur, Rr, Vr] = svd_ranktol(rH2,rank);
    [Uc, Rc, Vc] = svd_ranktol(cH2,rank);
    
    rH2=[];
    cH2=[];
    
    %store U, V, D
    hss.U = Ur;
    hss.V = Vc;
    hss.D = f(m_idx,m_idx,N);
    
    %store dimension of current diagonal block m and if leaf node
    hss.m = m;
    hss.leaf = tree.leaf;
    hss.root = tree.root;
    
    %return variables, v, d and row and col hankel blocks
    V = hss.V;
    
    %increment counter
    %d = d + m;
    
    %return hankel blocks with U removed from the row hankel block
    %and V' removed from the column hankel block.
    rH = Rr*Vr';
    cH = Uc*Rc;
    
else %If not leaf node
    
    %Descend into the tree, depth first
    if tree.treeL.depth >= tree.treeR.depth %left node deeper than right
        %left recursive call
        [hssL, upperRow, leftCol] = HSS_Basis_TranslationOp(tree.treeL,N,rank);
        
        %right recursive call
        [hssR, lowerRow, rightCol] = HSS_Basis_TranslationOp(tree.treeR,N,rank);
    else  %right node deeper than left
        
        %right recursive call
        [hssR, lowerRow, rightCol] = HSS_Basis_TranslationOp(tree.treeR,N,rank);
        
        %left recursive call
        [hssL, upperRow, leftCol] = HSS_Basis_TranslationOp(tree.treeL,N,rank);
    end
    
    %starting index (row/col value) of its current diagonal block.
    d = tree.d;
    
    %indexes of rol/col hankel blocks excluding corresponding part of diagonal block
    uR_idx = [1:(d-1) (d+hssR.m):(N-hssL.m)];
    lR_idx = [1:(d-1) (d+hssL.m):(N-hssR.m)];
    
    rH_top = upperRow(:, uR_idx);
    rH_bottom = lowerRow(:, lR_idx);
    cH_left = leftCol(uR_idx,:);
    cH_right = rightCol(lR_idx,:);
    
    rH2 = [rH_top; rH_bottom];
    cH2 = [cH_left cH_right];
    
    %take svd of remaining portions of blocks to determine Rs and Ws
    [Ur, Rr, Vr]= svd_ranktol(rH2,rank);
    [Uc, Rc, Vc]= svd_ranktol(cH2,rank);
    
    %partition Us to get R's, partition V's to get Ws
    hss.Rl = Ur(1:size(upperRow,1),:);
    hss.Rr = Ur(size(upperRow,1)+1:end,:);
    hss.Wl = Vc(1:size(leftCol,2),:);
    hss.Wr = Vc(size(leftCol,2)+1:end,:);
    
    %store dimension of current diagonal block and if leaf node
    hss.m = tree.m;
    hss.leaf = tree.leaf;
    hss.root= tree.root;
    
    %store recursive structure
    hss.hssL = hssL;
    hss.hssR = hssR;
    
    %return relavant portions of hankel blocks
    rH = Rr*Vr';
    cH = Uc*Rc;
    
end


end

function [Bu, Bl] = HSS_Expansion_Coeffs_OffDiag(treeL,treeR,row_idx,col_idx,N)
%Computes Expansion Coefficients (B's) of HSS structure for the upper right
%and lower left block of the current node
%   -------
%  |   | x |
%  |-------|
%  | x |   |
%   -------
%INPUT:             treeL   (struct) left branch of current node
%                   treeR   (struct) right branch of current node
%                   row_idx (int) row index of current node
%                   col_idx (int) col index of current node
%                   N       (int) grid size
%OUTPUT:            Bu      (matrix) Expansion coefficient matrix B at the
%                           current level which corresponds to the upper
%                           right block
%                   Bl      (matrix) Expansion coefficient matrix B at the
%                           current level which corresponds to the lower
%                           left block

if treeL.leaf && treeR.leaf %If node is a leaf
    m_idx = (row_idx):(row_idx+treeL.m-1);
    butm_idx = col_idx:(col_idx+treeR.m-1);
    
    Au = f(m_idx,butm_idx,N);
    Bu = treeL.U'*Au*treeR.V;
    Au = [];
    
    Al = f(butm_idx,m_idx,N);
    Bl = treeR.U'*Al*treeL.V;
    Al = [];
    
elseif ~treeL.leaf && ~treeR.leaf %If neither node is a leaf
    
%     treeLL = treeL.hssL;
%     treeLR = treeL.hssR;
%     treeRL = treeR.hssL;
%     treeRR = treeR.hssR;
    
    %COMPUTE EXPANSION COEFFICIENTS 
    
    %upper left corner of both upper right and lower left blocks
    %(these are computed at the same time. One block will be of 
    %dimension MxN, the other will always be of dimension NxM)
    %
    %    ---------------
    %   |       | x |<--|-- MxN
    %   |        -------|  
    %   |       |   |   |
    %   |---------------|
    %   | x |<--|-------|-- NxM
    %   |-------|       |
    %   |   |   |       |
    %    ---------------
    %           (1)
    
    [Bu_ll, Bl_ll] = HSS_Expansion_Coeffs_OffDiag(treeL.hssL,treeR.hssL,row_idx,col_idx,N);
    bu_ll = treeL.Rl'*Bu_ll*treeR.Wl;
    bl_ll = treeR.Rl'*Bl_ll*treeL.Wl;
    Bu_ll = [];
    Bl_ll = [];
    
    %Lower left corner lower left block, and also the upper right corner of
    %the upper right block
    %    ---------------
    %   |       |   | x |
    %   |       |---|---|
    %   |       |   |   |
    %   |---------------|
    %   |   |   |       |
    %   |-------|       |
    %   | x |   |       |
    %    ---------------
    %           (2)
    
    col_idx = col_idx +treeR.hssL.m;
    [Bu_lr, Bl_lr] = HSS_Expansion_Coeffs_OffDiag(treeL.hssL,treeR.hssR,row_idx,col_idx,N);
    bu_lr = treeL.Rl'*Bu_lr*treeR.Wr;
    bl_lr = treeR.Rr'*Bl_lr*treeL.Wl;
    Bu_lr = [];
    Bl_lr = [];
    
    %upper right corner of lower block, lower left corner of upper block
    %    ---------------
    %   |       |   |   |
    %   |       |---|---|
    %   |       | x |   |
    %   |---------------|
    %   |   | x |       |
    %   |-------|       |
    %   |   |   |       |
    %    --------------- 
    %           (3)

    col_idx = col_idx -treeR.hssL.m;
    row_idx = row_idx +treeL.hssL.m;
    [Bu_rl, Bl_rl] = HSS_Expansion_Coeffs_OffDiag(treeL.hssR,treeR.hssL,row_idx,col_idx,N);
    bu_rl = treeL.Rr'*Bu_rl*treeR.Wl;
    bl_rl = treeR.Rl'*Bl_rl*treeL.Wr;
    Bu_rl = [];
    Bl_rl = [];
    
    %lower right corner of both lower left and upper left blocks 
    %    ---------------
    %   |       |   |   |
    %   |       |---|---|
    %   |       |   | x |
    %   |---------------|
    %   |   |   |       |
    %   |-------|       |
    %   |   | x |       |
    %    ---------------    
    %           (4)
    
    col_idx = col_idx +treeR.hssL.m;
    [Bu_rr, Bl_rr] = HSS_Expansion_Coeffs_OffDiag(treeL.hssR,treeR.hssR,row_idx,col_idx,N);
    bu_rr = treeL.Rr'*Bu_rr*treeR.Wr;
    bl_rr = treeR.Rr'*Bl_rr*treeL.Wr;
    Bu_rr = [];
    Bl_rr = [];
    
    %Add b's to form the B's
    Bu = bu_ll+bu_rl+bu_lr+bu_rr;
    Bl = bl_ll+bl_rl+bl_lr+bl_rr;
    
    bu_ll = [];
    bu_rl = [];
    bu_lr = [];
    bu_rr = [];
    bl_ll = [];
    bl_rl = [];
    bl_lr = [];
    bl_rr = [];
    
elseif treeL.leaf && ~treeR.leaf %If left node is a leaf, and right node is not
    
%     treeRL = treeR.hssL;
%     treeRR = treeR.hssR;
    
    %COMPUTE EXPANSION COEFFICIENTS
    
    %left block of upper right corner and upper block of lower left corner 
    %(these are computed in the same pass)
    %    ---------------
    %   |       |   |   |
    %   |       | x |<--|-- MxN  
    %   |       |   |   |
    %   |---------------|
    %   |   x<--|-------|-- NxM
    %   |-------|       |
    %   |       |       |
    %    ---------------
    %           (5)
    
    [Bu_l, Bl_l] = HSS_Expansion_Coeffs_OffDiag(treeL,treeR.hssL,row_idx,col_idx,N); 
    bu_l = Bu_l*treeR.Wl;
    bl_l = treeR.Rl'*Bl_l;
    Bu_l = [];
    Bl_l = [];

    %right block of upper right corner, upper left block of lower left corner.
    %    ---------------
    %   |       |   |   |
    %   |       |   | x	|  
    %   |       |   |   |
    %   |---------------|
    %   |       |       |	
    %   |-------|       |
    %   |   x   |       |
    %    ---------------
    %           (6)

    col_idx = col_idx + treeR.hssL.m;
    [Bu_r, Bl_r] = HSS_Expansion_Coeffs_OffDiag(treeL,treeR.hssR,row_idx,col_idx,N);
    col_idx = col_idx - treeR.hssL.m;
    bu_r = Bu_r*treeR.Wr;
    bl_r = treeR.Rr'*Bl_r;
    Bu_r = [];
    Bl_r = [];
    
    Bu = bu_l + bu_r;
    Bl = bl_l + bl_r;
    
    bu_l = [];
    bu_r = [];
    bl_l = [];
    bl_r = [];
elseif ~treeL.leaf && treeR.leaf %If right node is a leaf, and left node is not
%     treeLL = treeL.hssL;
%     treeLR = treeL.hssR;
    
    %COMPUTE EXPANSION COEFFICIENTS
    
    %upper block of upper right corner and left block of lower left corner 
    %(these are computed in the same pass)
    %    ---------------
    %   |       |   x <-|-- MxN
    %   |       |-------|  
    %   |       |       |
    %   |---------------|
    %   |   |   |       |
    %   | x |<--|-------|-- NxM
    %   |   |   |       |
    %    ---------------
    %           (7)
    
    [Bu_l, Bl_l] = HSS_Expansion_Coeffs_OffDiag(treeL.hssL,treeR,row_idx,col_idx,N); 
    bu_l = treeL.Rl'*Bu_l;
    bl_l = Bl_l*treeL.Wl;
    Bu_l = [];
    Bl_l = [];
    
    %upper block of upper right corner and left block of lower left corner 
    %(these are computed in the same pass)
    %    ---------------
    %   |       |       |
    %   |       |-------|  
    %   |       |   x   |
    %   |---------------|
    %   |   |   |       |
    %   |   | x |       |
    %   |   |   |       |
    %    ---------------
    %           (8)

    row_idx = row_idx + treeL.hssL.m;
    [Bu_r, Bl_r] = HSS_Expansion_Coeffs_OffDiag(treeL.hssR,treeR,row_idx,col_idx,N); 
    row_idx = row_idx - treeL.hssL.m;
    bu_r = treeL.Rr'*Bu_r;
    bl_r = Bl_r*treeL.Wr;
    Bu_r = [];
    Bl_r = [];
    
    Bu = bu_l + bu_r;
    Bl = bl_l + bl_r;
    
    bu_l = [];
    bu_r = [];
    bl_l = [];
    bl_r = [];
end

end

function [hss]= HSS_Expansion_Coeffs_Diag(tree,row_idx,col_idx,N)
%Recursively computes Expansion Coefficients (B's) of HSS structure for the
%upper left and lower right block of the current node
%   -------
%  | x |   |
%  |-------|
%  |   | x |
%   -------
%INPUT:             tree    (struct) branch of current node
%                   row_idx (int) row index of current node
%                   col_idx (int) col index of current node
%                   N       (int) grid size
%OUTPUT:            hss     (structure) contains matrices (Us, Vs, Bs Rs and
%                           Ws) that compose the HSS representation of the
%                           input tree for a corresponding matrix, as well
%                           as corresponding partition dimensions for each
%                           node

%EX:
%HSS_Expansion_Coeffs_OffDiag(A2;1,2,A2;2,1)
%HSS_Expansion_Coeffs_Diag(A2;1,2, A2;1,1)
%HSS_Expansion_Coeffs_Diag(A2;1,2, A2;2,2)


if tree.hssL.leaf  && tree.hssR.leaf %if both nodes are leaves
    col_idx = col_idx + tree.hssL.m; 
    m_idx = (row_idx):(row_idx+tree.hssL.m-1);
    butm_idx = col_idx:(col_idx+tree.hssR.m-1);
    
    Au = f(m_idx,butm_idx,N);
    Bu = tree.hssL.U'*Au*tree.hssR.V;
    Au = [];
    
    Al = f(butm_idx,m_idx,N);
    Bl = tree.hssR.U'*Al*tree.hssL.V;
    Al = [];
    
    tree.Bu = Bu;
    tree.Bl = Bl;
    hss = tree;
    
elseif ~tree.hssL.leaf  && ~tree.hssR.leaf %neither node is a leaf
    
    %Calculate B's for off-diagonal blocks
    col_idx = col_idx + tree.hssL.m;
    [Bu, Bl] = HSS_Expansion_Coeffs_OffDiag(tree.hssL,tree.hssR,row_idx,col_idx,N);
    col_idx = col_idx - tree.hssL.m;
    
    hss = tree;
    hss.Bl = Bl;
    hss.Bu = Bu;
    
    %Calculate B's for upper left diagonal block
    [hssL] = HSS_Expansion_Coeffs_Diag(tree.hssL,row_idx,col_idx,N);
    
    %Calculate B's for lower right diagonal block
    row_idx = row_idx + tree.hssL.m;
    col_idx = col_idx + tree.hssL.m;
    [hssR] = HSS_Expansion_Coeffs_Diag(tree.hssR,row_idx,col_idx,N);
    
    hss.hssL = hssL;
    hss.hssR = hssR;
    
elseif tree.hssL.leaf  && ~tree.hssR.leaf %left node is a leaf, right node is not
    
    %Calculate B's for off-diagonal blocks
    
    %left block of upper right corner and upper block of lower left corner 
    %(these are computed in the same pass)
    %    ---------------
    %   |       |   |   |
    %   |       | x |   |  
    %   |       |   |   |
    %   |---------------|
    %   |   x   |   |   |
    %   |-------|---|---|
    %   |       |   |   |
    %    ---------------
    %           (9)
    
    col_idx = col_idx + tree.hssL.m;
    [Bu_l, Bl_l] = HSS_Expansion_Coeffs_OffDiag(tree.hssL,tree.hssR.hssL,row_idx,col_idx,N);
    bu_l = Bu_l*tree.hssR.Wl;
    bl_l = tree.hssR.Rl'*Bl_l;
    Bu_l = [];
    Bl_l = [];
    
    %right block of upper right corner and lower block of lower left corner 
    %(these are computed in the same pass)
    %    ---------------
    %   |       |   |   |
    %   |       |   | x |  
    %   |       |   |   |
    %   |---------------|
    %   |       |   |   |
    %   |-------|---|---|
    %   |   x   |   |   |
    %    ---------------
    %           (10)

    col_idx = col_idx + tree.hssR.hssL.m;
    [Bu_r, Bl_r] = HSS_Expansion_Coeffs_OffDiag(tree.hssL,tree.hssR.hssR,row_idx,col_idx,N);
    col_idx = col_idx - tree.hssL.m -tree.hssR.hssL.m;
    bu_r = Bu_r*tree.hssR.Wr;
    bl_r = tree.hssR.Rr'*Bl_r;
    Bu_r = [];
    Bl_r = [];
    
    hss = tree;
    hss.Bl = bl_l + bl_r;
    hss.Bu = bu_l + bu_r;
    
    bl_l = [];
    bl_r = [];
    bu_l = [];
    bu_r = [];
    
    %Calculate B's for lower right diagonal block
    %upper right block of lower right corner and lower left block of lower 
    % right corner 
    %    ---------------
    %   |       |   |   |
    %   |       |   |   |  
    %   |       |   |   |
    %   |---------------|
    %   |       |   | x |
    %   |-------|---|---|
    %   |       | x |   |
    %    ---------------
    %           (11)
    
    row_idx = row_idx + tree.hssL.m; 
    col_idx = col_idx + tree.hssL.m; 
    [hssR] = HSS_Expansion_Coeffs_Diag(tree.hssR,row_idx,col_idx,N);
    
    hss.hssR = hssR;
    
elseif ~tree.hssL.leaf  && tree.hssR.leaf %right node is a leaf, left node is not
    
    %Calculate B's for off-diagonal blocks
    
    %upper block of upper right corner and left block of lower left corner 
    %(these are computed in the same pass)
    %    ---------------
    %   |   |   |   x   |
    %   |---|---|-------|  
    %   |   |   |       |
    %   |---------------|
    %   |   |   |       |
    %   | x |   |       |
    %   |   |   |       |
    %    ---------------
    %           (12)
    
    col_idx = col_idx + tree.hssL.m;
    [Bu_l, Bl_l] = HSS_Expansion_Coeffs_OffDiag(tree.hssL.hssL,tree.hssR,row_idx,col_idx,N);
    bu_l = tree.hssL.Rl'*Bu_l;
    bl_l = Bl_l*tree.hssL.Wl;
    Bu_l = [];
    Bl_l = [];
    
    %lower block of upper right corner and right block of lower left corner 
    %(these are computed in the same pass)
    %    ---------------
    %   |   |   |       |
    %   |---|---|-------|  
    %   |   |   |   x   |
    %   |---------------|
    %   |   |   |       |
    %   |   | x |       |
    %   |   |   |       |
    %    ---------------
    %           (13)
    
    row_idx  = row_idx + tree.hssL.hssL.m;
    [Bu_r, Bl_r] = HSS_Expansion_Coeffs_OffDiag(tree.hssL.hssR,tree.hssR,row_idx,col_idx,N);
    row_idx  = row_idx - tree.hssL.hssL.m;
    col_idx = col_idx - tree.hssL.m;
    bu_r = tree.hssL.Rr'*Bu_r;
    bl_r = Bl_r*tree.hssL.Wr;
    Bu_r = [];
    Bl_r = [];
    
    hss = tree;
    hss.Bl = bl_l + bl_r;
    hss.Bu = bu_l + bu_r;
    
    bl_l = [];
    bl_r = [];
    bu_l = [];
    bu_r = [];
    
    %Calculate B's for upper right diagonal block
    
    %lower left block of upper left corner and upper right block of upper 
    % left corner 
    %    ---------------
    %   |   | x |       |
    %   |---|---|-------|  
    %   | x |   |       |
    %   |---------------|
    %   |   |   |       |
    %   |   |   |       |
    %   |   |   |       |
    %    ---------------
    %           (14)
    
    [hssL] = HSS_Expansion_Coeffs_Diag(tree.hssL,row_idx,col_idx,N);
    
    hss.hssL = hssL;
    
end

end


function [Uhat, Rhat, Vhat] = svd_ranktol(A,rank)
%INPUT:             A       matrix
%                   rank    (int) largest allowable rank of hankel blocks -
%                           this determines the amount of compression, and
%                           should be compatible with your input tree.
%                           (if the maximum allowable rank is p, then the
%                           tree partitions should be no smaller than 3p)
%                           Corresponding singular vectors below this rank
%                           will be dropped.
%OUTPUT:            Uhat    (array) left singular vectors of A
%                           corresponding to singular values higher than
%                           the chosen tolerance
%                   Vhat    (array) right singular vectors of A
%                           corresponding to singular values higher than
%                           the chosen tolerance
%                   Rhat    (array) diagonal matrix containing the singular
%                           values of A higher than the chosen tolerance

[U, R, V] = svd(A,'econ');
if isempty(U)
    Uhat = U;
    Rhat = R;
    Vhat = V;
else
    Uhat = U(:,1:rank);
    Rhat = R(1:rank,1:rank);
    Vhat = V(:,1:rank);
end

end

function [tree,d] = label_tree(tree,d)
%This function determines and labels the maximum depth of each node in a 
%tree, as well as label each node with the starting index (row/col value)
% of its corresponding diagonal block.
%
%INPUT:             tree    (structure) contains partition dimensions for 
%                           each subdivision of the input matrix, for each
%                           of which there are a left and right child
%                           structure
%                   d       (int) row and col index which corresponds to 
%                           the first element in the current diagonal block 
%                           at every node.
%OUTPUT:            tree    (structure) same as input, but with an added
%                           field that contains tree depth at each node
% Author: Kristen Lessel - Sept 2014

if tree.leaf == 0 % if not a leaf node
    
    tree.d = d;
    [tree.treeL,d] = label_tree(tree.treeL,d);
    depth_l = tree.treeL.depth;
    
    [tree.treeR,d] = label_tree(tree.treeR,d);
    depth_r = tree.treeR.depth;
    
    tree.depth = 1 + max(depth_l,depth_r);

else
    tree.depth = 0;
    tree.d = d;
    d = d + tree.m;
end
end

function H = f(rowIdx,colIdx,N)
a = 0;
b = 1;

%linearly spaced x vector
x = linspace(a,b,N);

%logarithmicly spaced vector
% x = logspace(a,b,N);
% x = (x-1)/9;

ml = length(rowIdx);
nl = length(colIdx);
H = sqrt(abs(repmat(x(rowIdx)',1,nl)-repmat(x(colIdx),ml,1)));

end

% function [H] = f(rowIdx,colIdx,N)
% %generates blocks of a matrix defined by input parameters. This function
% is not efficient in matlab
% %INPUT          sr: (scalar) start row
% %               sc: (scalar) start column
% %               ml: (scalar) number of rows
% %               nl: (scalar) number of columns
% %               sInt: (scalar) start of interval
% %               eInt: (scalar) end of interval
% %               N: (scalar) grid size
% 
% %OUTPUT         H: (matrix) matrix defined by given inputs
% sInt = 0;
% eInt = 1;
% 
% % %start of interval
% % sr = rowIdx(1);
% % %end of interval
% % sc = colIdx(1);
% %number of rows
% ml = length(rowIdx);
% %number of columns
% nl = length(colIdx);
% 
% if ml == 0 || nl == 0
%     H = [];
% else
%     x = linspace(sInt,eInt,N);
%     for ii = 1:ml
%         for jj = 1:nl
%             i_idx = rowIdx(ii);
%             j_idx = colIdx(jj);
%             H(ii,jj) = sqrt(abs(x(i_idx)-x(j_idx)));
%             %H(ii,jj) = log(1+abs(x(ii+sr)-x(jj+sc)));
%         end
%     end
% end
% 
% 
% end
