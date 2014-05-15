function [Iout] = TVdenoise(I_in, varargin)
%function [Iout] = TVdenoise(I_in, varargin)
%
%Total Variation based image restoration with 4 neighbours (non-linear diffusion).
%
%INPUT
%Iin
%
%PARAMETERS
%see from the code
%
%OUTPUT
%Iout
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

param.alpha = 5;		%Smoothness
param.omega = 1.75;		%Relaxation parameter (SOR).
param.outer_iter = 10;		%Number of outer iterations (lagged diffusivity iterations)
param.inner_iter = 5;		%Number of inner iterations (solver iterations).
param.solver = 2.0;		%Type of solver used:
					%1 = Normal point-wise Gauss-Seidel relaxation.
					%2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.
param.scl = 0.5;		%Downscale to this size
param.scl_factor = 0.75;	%Downscale factor

param = setParameters(param, varargin{:});

I_in = single(I_in);

%Current size and downscale to size
[rows cols frames] = size(I_in);
ds_rows = ceil(rows*param.scl);
ds_cols = ceil(cols*param.scl);
scales = double(intmax);

%Smoothening mask
G = fspecial('gaussian',[7 7],2.0);
%G = fspecial('gaussian',[5 5],1.25);
Iin{1} = I_in;
isizes{1} = size(Iin{1});
for scl=2:scales
	Iin{scl} = imresize( Iin{scl-1}, param.scl_factor, 'bilinear' );
	isizes{scl} = size(Iin{scl});
	%----------
	% Smoothen
	%----------
	Iin{scl-1} = imfilter( Iin{scl-1}, G, 'replicate' );
	%If scaled image size is less than 10, break
	if isizes{scl}(1)<=ds_rows | isizes{scl}(2)<=ds_cols
		scales = scl;
		%---------------------
		% Smoothen last scale
		%---------------------
		Iin{scl} = imfilter( Iin{scl}, G, 'replicate' );
		break
	end
end
Iout = Iin{scales};

for scl=scales:-1:1

	for iter=0:param.outer_iter

		PsiData = 1./realsqrt( (Iout-Iin{scl}).^2 + eps );
		[wW wN wE wS] = DiffWeights(Iout); 

		TRACE = PsiData + param.alpha*(wW+wN+wE+wS);
		B = PsiData.*Iin{scl};
	  
		[Iout] =	PDEsolver4(	single(Iout), ...
						single(TRACE), ...
						single(B), ...
						single(param.alpha*wW), ...
						single(param.alpha*wN), ...
						single(param.alpha*wE), ...
						single(param.alpha*wS), ...
						single(param.inner_iter), ...		%solver iterations
						single(param.omega), ...		%omega
						single(param.solver) );			%solver
		%mesh( double(Iout) );
		%drawnow
		
	end

	%---------
	% Upscale
	%---------
	if (scl-1)>0
		Iout = imresize( Iout, 'bilinear' ,'OutputSize', isizes{scl-1}(1:2)  );
	end
end

%Iout(Iout<0)=0;
%Iout(Iout>1)=1;

function [wW wN wE wS] = DiffWeights(D)

	  [rows cols frames] = size( D );

	  %Spatial differences
	  Dver = imfilter(D,[0.25 0 -0.25]','replicate');
	  Dhor = imfilter(D,[0.25 0 -0.25],'replicate');
	  
	  wW = (circshift(D,[0 1 0])-D).^2 + (Dver+circshift(Dver,[0 1 0])).^2;
	  wE = (circshift(D,[0 -1 0])-D).^2 + (Dver+circshift(Dver,[0 -1 0])).^2;
	  wN = (circshift(D,[1 0 0])-D).^2 + (Dhor+circshift(Dhor,[1 0 0])).^2;
	  wS = (circshift(D,[-1 0 0])-D).^2 + (Dhor+circshift(Dhor,[-1 0 0])).^2;
	  
	  %Find the maximum derivative...is much more effective at stopping "bleeding" at the edges
	  wW = max(wW,[],3);
	  wE = max(wE,[],3);
	  wN = max(wN,[],3);
	  wS = max(wS,[],3);
  
	  wW = 1./realsqrt( wW + 0.00001 );
	  wE = 1./realsqrt( wE + 0.00001 );
	  wN = 1./realsqrt( wN + 0.00001 );
	  wS = 1./realsqrt( wS + 0.00001 );
	  
	  %wW = 1./( wW + 0.001 ).^0.3;
	  %wE = 1./( wE + 0.001 ).^0.3;
	  %wN = 1./( wN + 0.001 ).^0.3;
	  %wS = 1./( wS + 0.001 ).^0.3;

	  wW(:,1,:) = 0;
	  wE(:,end,:) = 0;
	  wN(1,:,:) = 0;
	  wS(end,:,:) = 0;

	  wW = repmat( wW, [1 1 frames] );
	  wN = repmat( wN, [1 1 frames] );
	  wE = repmat( wE, [1 1 frames] );
	  wS = repmat( wS, [1 1 frames] );

