function [b] = HSS_Multiply(hss,x)
%function [b] = HSS_Multiply(hss,x)
%Takes in an hss matrix A_hss, and multiplies this with an input 
%array/matrix, x.  Gives an array/matrix b as output.  A_hss*x = b.  
%
%INPUT:             hss   (struct) a matrix in HSS form as given by the
%                           function Gen_HSS_Matrix().
%                   x       (array/matrix) input array/matrix to multiply with A_hss.
%
%OUTPUT:            b       (array/matrix) output array/matrix from the
%                           resulting multiplication
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


% d = 1;
f= zeros(size(hss.Rr,2),size(x,2)); %size(hss.Rr,2) = size(hss.Rl,2)

gtree = get_gtree(hss,x);

[b] = HSS_Multiply_iter(hss,x,f,gtree);
end

function [gtree,g] = get_gtree(hss,x)

if hss.leaf
    g = hss.V'*x;
    gtree= struct;
    
else
    
    xl = x(1:hss.hssL.m,:);
    [gtreeL, gl] = get_gtree(hss.hssL,xl);
    
    xr = x(hss.hssL.m+1:end,:);
    [gtreeR, gr] = get_gtree(hss.hssR,xr);
    
    %new g
    g = hss.Wl'*gl+hss.Wr'*gr;
    
    gtreeL.g = gl;
    gtreeR.g = gr;
    
    gtree.gtreeL = gtreeL;
    gtree.gtreeR = gtreeR;
end

end

%gtree is hss structure with g's stored at each  node as well
function [b] = HSS_Multiply_iter(hss,x,f,gtree)

%if leaf node
if hss.leaf
    
    b = hss.D*x + hss.U*f;

else
    
    xl = x(1:hss.hssL.m,:);
    fl = hss.Br*gtree.gtreeR.g + hss.Rl*f;
    [bl] = HSS_Multiply_iter(hss.hssL,xl,fl,gtree.gtreeL);
    
    
    xr = x(hss.hssL.m+1:end,:);
    fr = hss.Bl*gtree.gtreeL.g + hss.Rr*f; 
    [br] = HSS_Multiply_iter(hss.hssR,xr,fr,gtree.gtreeR);
    
    
    b = [bl; br];
    
end

end
