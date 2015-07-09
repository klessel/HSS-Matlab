function [A] = Build_Matrix_From_HSS2(hss)
%Takes in an HSS matrix (structure) and returns the dense matrix A 
%(nxn double). This code should be use for error checking/debugging only, 
%in practice you should not construct the full matrix A since the amount of
%memory required to do so will be large. Uses feild labled Bu instead of Br
%
%INPUT:             hss    (structure) contains matrices (Us, Vs, Bs Rs and
%                           Ws) that compose the hss representation of the
%                           input matrix, as well as partition dimensions
%                           and if branch is a leaf.
%
%OUTPUT:            A       dense matrix A 
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

[A] = Build_Matrix_From_HSS_iter(hss);

end

function [D, U, V] = Build_Matrix_From_HSS_iter(hss)

if hss.leaf %if leaf node
    U = hss.U;
    V = hss.V;
    D = hss.D;
    
else
    
    %left recursive call
    [Dl, Ul, Vl] = Build_Matrix_From_HSS_iter(hss.hssL);
    
    %right recursive call
    [Dr, Ur, Vr] = Build_Matrix_From_HSS_iter(hss.hssR);
    
    %construct upper right and lower left block of current level
    urb = Ul*hss.Bu*Vr';
    llb = Ur*hss.Bl*Vl';
    
    %create D, U and V to pass to parent
    D = [Dl urb; llb Dr];
    U = [Ul*hss.Rl; Ur*hss.Rr];
    V = [Vl*hss.Wl; Vr*hss.Wr];
    
end

end
