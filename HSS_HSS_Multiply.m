function [hssC] = HSS_HSS_Multiply(hssA,hssB)
%function [hssC] = HSS_HSS_Multiply(hssA,hssB)
%Takes in an hss representation of a matrix A (hssA) and multiplies this with 
% the hss representation of another matrix B (hssB).  The output is the hss 
% representation of a matrix C (hssC), such that A*B=C.  
%
%INPUT:             hssA   (struct) a matrix in HSS form as given by the
%                           function Gen_HSS_Matrix().
%                   hssB   (struct) a matrix in HSS form as given by the
%                           function Gen_HSS_Matrix().
%
%OUTPUT:            hssC    (struct) a matrix in HSS form; a result of the
%                          multiplication of hss1*hss2.
%
% Author: Kristen Lessel, klessel@engineering.ucsb.edu
%
% Copyright (C) 2013 Kristen Lessel, (klessel@engineering.ucsb.edu)
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


ftree.f = zeros(size(hssA.Rr,2),size(hssB.Wl,2)); %size(hss.Rr,2) = size(hss.Rl,2)

%downsweep to calculate all of the g's
[gtree, ~] = get_gtree(hssA,hssB);

%Main call - contains an upsweep to calculate f's
[hssC,ftree] = HSS_HSS_Multiply_iter(hssA,hssB,ftree,gtree);
end

function [gtree,g] = get_gtree(hssA,hssB)

if hssA.leaf %if leaf node
    g = hssA.V'*hssB.U; %You must ensure the 2nd dimension of U from hssA and V from hssB match.
    gtree= struct;
    
else %if not leaf node

    [gtreeL, gl] = get_gtree(hssA.hssL,hssB.hssL);
    
    [gtreeR, gr] = get_gtree(hssA.hssR,hssB.hssR);
    
    %Compute new g
    g = hssA.Wl'*gl*hssB.Rl +hssA.Wr'*gr*hssB.Rr;
    
    gtreeL.g = gl;
    gtreeR.g = gr;
    
    %store g's in tree
    gtree.gtreeL = gtreeL;
    gtree.gtreeR = gtreeR;
end

end

%gtree is hss structure with g's stored at each  node as well
function [hssC, ftree] = HSS_HSS_Multiply_iter(hssA,hssB,ftree,gtree)

%if leaf node
if hssA.leaf
    
    %Calculate U, V and D at the leaf nodes
    hssC.U = [hssA.U  hssA.D*hssB.U];
    hssC.V = [hssB.D*hssA.V hssB.V];
    hssC.D = hssA.D*hssB.D + hssA.U*ftree.f*hssB.V';
    hssC.leaf = 1;
    
else
    %compute f
    ftreeL.f = hssA.Br*gtree.gtreeR.g*hssB.Bl + hssA.Rl*ftree.f*hssB.Wl';
    ftreeR.f = hssA.Bl*gtree.gtreeL.g*hssB.Br + hssA.Rr*ftree.f*hssB.Wr';
    
    %left call
    [hssL, ftreeL] = HSS_HSS_Multiply_iter(hssA.hssL,hssB.hssL,ftreeL,gtree.gtreeL); 
    
    %right call
    [hssR, ftreeR] = HSS_HSS_Multiply_iter(hssA.hssR,hssB.hssR,ftreeR,gtree.gtreeR);
    
    %store f's in tree
    ftree.ftreeL = ftreeL;
    ftree.ftreeR = ftreeR;
    
    %store recursive HSS structure
    hssC.hssL = hssL;
    hssC.hssR = hssR;
    
    %Store new B's
    hssC.Bl = [hssA.Bl hssA.Rr*ftree.f*hssB.Wl'; zeros(size(hssB.Bl,1),size(hssA.Bl,2)) hssB.Bl];
    hssC.Br = [hssA.Br hssA.Rl*ftree.f*hssB.Wr'; zeros(size(hssB.Br,1), size(hssA.Br,2)) hssB.Br];
    
    %Store new W's
    hssC.Wl = [hssA.Wl zeros(size(hssA.Wl,1),size(hssB.Wl,2)); hssB.Bl'*gtree.gtreeR.g'*hssA.Wr hssB.Wl];
    hssC.Wr = [hssA.Wr zeros(size(hssA.Wr,1),size(hssB.Wr,2)); hssB.Br'*gtree.gtreeL.g'*hssA.Wl hssB.Wr];
    
    %Store new R's
    hssC.Rl = [hssA.Rl hssA.Br*gtree.gtreeR.g*hssB.Rr; zeros(size(hssA.Rl,1),size(hssB.Rl,2)) hssB.Rl];
    hssC.Rr = [hssA.Rr hssA.Bl*gtree.gtreeL.g*hssB.Rl; zeros(size(hssA.Rr,1),size(hssB.Rr,2)) hssB.Rr];
    
    hssC.leaf = 0;
    
end

end
