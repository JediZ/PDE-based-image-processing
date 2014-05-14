function [PHIout] = GAC(Iin, PHIin, varargin)
%function [PHIout] = GAC(Iin, PHIin, varargin)
%
%Geodesic Active Contour model: Active Contour model based on Level-sets.
%This function implements the geodesic active contour model of Caselles et al. 1993.
%
%Model:
%$\Phi_t = |\nabla \Phi| DIV( \dfrac{g(I)}{|\nabla \Phi|} \nabla \Phi ) +  g(I) |\nabla \Phi| c$
%where c is the 'balloon force'
%
%Author: Jarno Ralli
%E-mail: jarno@ralli.fi
%Web: www.jarnoralli.fi (or www.jarnoralli.com)
%
%If you use this code, please reference (some) of my papers available at http://www.jarnoralli.fi
%
%The program is delivered as it is and the author accepts no responsibility what so ever of its use and/or results.
%Errors and suggestions are kindly to be communicated to the author.
%
%(C) Copyright 2011-2014, Jarno Ralli
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU Lesser General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU Lesser General Public License for more details.
%
%You should have received a copy of the GNU Lesser General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.

param.tau = 0.25;		%Time step.
param.PHI = [];			%Level-set function as input.
param.c = -0.1;			%'Balloon' force (negative => shrinks, positive => expands)
param.lambda = -1;		%Parameter for defining edge strengths in the 'stopping function' g
				%	-an image gradient smaller than lambda has g>0.5
				%	-an image gradient bigger than lambda has g<0.5
				%	-negative value => approximate automatically (see the code)
param.ITER = 100;		%Number of iterations
param.SMOOTH = 100;		%Smoothness of the curve
param = setParameters(param, varargin{:});

%Reinitialise the level-set funcion to a signed distance function before starting
PHIout = Reinit(single(PHIin), single(10));

%Operators
O_dx = [-1 0 1]*0.5;
O_dy = O_dx';
O_dx_fd = [0 -1 1];		%Forward differencing
O_dx_bd = [-1 1 0];		%Backward differencing
O_dy_fd = [0 -1 1]';		%Forward differencing
O_dy_bd = [-1 1 0]';		%Backward differencing
G = fspecial('gaussian',[7 7],2.5);

%Smoothen input image and approximate derivatives
I = imfilter(Iin, G, 'replicate');
Idx = imfilter( I, O_dx, 'replicate' );
Idy = imfilter( I, O_dy, 'replicate' );
%In order to improve 'edge detection', use the biggest derivatives
Idx = max( Idx, [], 3);
Idy = max( Idy, [], 3);

%----------------------
%Stopping function 'g'
%----------------------
Igrad = Idx.^2 + Idy.^2;
%Evaluate lambda automatically
if param.lambda<0
    [Y I] = sort( Igrad(:) );
    param.lambda = Y( round( 0.7*length(Y) ) );
end
g = 1./( 1 + Igrad./param.lambda );

figure
iter = 0;
while iter<param.ITER

	%----------------------------
	%PHI => Level-set function
	%----------------------------
	%Calculate derivatives
	PHIdx = imfilter( PHIout, O_dx, 'replicate' );
	PHIdy = imfilter( PHIout, O_dy, 'replicate' );
	%'Upwind' derivatives
	PHIx_fd = imfilter( PHIout, O_dx_fd, 'replicate' );
	PHIx_bd = imfilter( PHIout, O_dx_bd, 'replicate' );
	PHIy_fd = imfilter( PHIout, O_dy_fd, 'replicate' );
	PHIy_bd = imfilter( PHIout, O_dy_bd, 'replicate' );

	if param.c<=0
		gradPHIUW = sqrt(	max(PHIx_bd,0).^2 + min(PHIx_fd,0).^2 + ...
					max(PHIy_bd,0).^2 + min(PHIy_fd,0).^2 );
	else
		gradPHIUW = sqrt(	min(PHIx_bd,0).^2 + max(PHIx_fd,0).^2 + ...
					min(PHIy_bd,0).^2 + max(PHIy_fd,0).^2 );
	end
	
	DATA = param.c*g.*gradPHIUW;
	
	%------------------------------------------
	% Update the level-set function Omega (PHI)
	%------------------------------------------
	%Gradient of PHI (for stability eps is added)
	gradPHI = sqrt(PHIdx.^2+PHIdy.^2+eps);
	Diff = gradPHI./g;			%gradPHI/g instead of g/gradPHI because of 'harmonic averaging'

	PHIout = AC_solver_2d(	single( PHIout ), ...
				single( DATA ), ...
				single( gradPHI ), ...
				single( Diff ), ...
				single( param.tau ), ...
				single( param.SMOOTH ) );
	
	imagesc(Iin),hold on, contour(PHIout>=0,'r-'), title([int2str(iter),'/',int2str(param.ITER)]), hold off
	drawnow	

	iter = iter+1;

end
