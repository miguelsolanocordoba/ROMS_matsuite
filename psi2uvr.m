function [ufield,vfield,rfield,pfield,pfield2]=psi2uvr(field);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   compute the values at u,v and psi points...
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Lp,Mp]=size(field);
M1=Mp-1;
L1=Lp-1;
%
pfield=field(2:end-1,2:end-1);
rfield=0.5*(field(1:L1,1:M1)+field(2:Lp,2:Mp)); 

%% From RHO to U/V
[Mr,Lr]=size(rfield');
M=Mr-1;
L=Lr-1;

vfield=0.5*(rfield(:,1:M)+rfield(:,2:Mr));
ufield=0.5*(rfield(1:L,:)+rfield(2:Lr,:));
pfield2=0.5*(ufield(:,1:M)+ufield(:,2:Mr));
%vfield=0.5*(rfield(1:M,:)+rfield(2:Mp,:));
%ufield=0.5*(rfield(:,1:L)+rfield(:,2:Lp));
%pfield2=0.5*(ufield(1:M,:)+ufield(2:Mp,:));

return













