function [PHIout] = GAC(Iin, PHIin, varargin)
%function [PHIout] = GAC(Iin, PHIin, varargin)
%
%Geodesic Active Contour model: Active Contour model based on Level-sets.
%This function implements the geodesic active contour model of Caselles et al. 1997.
%
%Model:
%$\Phi_t = |\nabla \Phi| DIV( \dfrac{g(I)}{|\nabla \Phi|} \nabla \Phi ) + \nabla g(I) \cdot \nabla \Phi$
%
%In other words, the model has a convection term.
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
G = fspecial('gaussian', [7 7], 2.5);

%Smoothen input image
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
gdx = imfilter( g, O_dx, 'replicate' );
gdy = imfilter( g, O_dy, 'replicate' );

fig = figure;
iter = 0;

while iter<param.ITER

	%----------------------------
	%PHIEGA => Level-set function
	%----------------------------
	%Calculate derivatives
	PHIdx = imfilter( PHIout, O_dx, 'replicate' );
	PHIdy = imfilter( PHIout, O_dy, 'replicate' );
	%Upwind derivatives for the convection term: (nabla g) dot (nabla PHI)
	%where (nabla g) is the 'speed'
	DATA =	max(0, gdx).*( circshift(PHIout,[0 -1])-PHIout ) + ...
		min(0, gdx).*( PHIout-circshift(PHIout,[0 1])  ) + ...
		max(0, gdy).*( circshift(PHIout,[-1 0])-PHIout ) + ...
		min(0, gdy).*( PHIout-circshift(PHIout,[1 0])  );
	
	%------------------------------------------
	% Update the level-set function Omega (PHI)
	%------------------------------------------
	%Gradient of PHI (for stability eps is added)
	gradPHI = sqrt(PHIdx.^2+PHIdy.^2+eps);
	%Diff = (gradPHI);
	Diff = gradPHI./g;			%gradPHI/g instead of g/gradPHI because of 'harmonic averaging'

	PHIout = AC_solver_2d(	single( PHIout ), ...
				single( DATA ), ...
				single( gradPHI ), ...
				single( Diff ), ...
				single( param.tau ), ...
				single( param.SMOOTH ) );

	imagesc(Iin),hold on, contour(PHIout>=0,'r-'), title([int2str(iter),'/',int2str(param.ITER)]), drawnow,
	hold off
	iter = iter+1;

end

close(fig)
