function [Iout] = TVdenoise(I_in, varargin)
%function [Iout] = TVdenoise(I_in, varargin)
%
%Total Variation based image restoration with 8 neighbours (anisotropic diffusion).
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

param.alpha = 500;		%Smoothness
param.omega = 1.75;		%Relaxation parameter (SOR).
param.outer_iter = 20;		%Number of outer iterations (lagged diffusivity iterations)
param.inner_iter = 4;		%Number of inner iterations (solver iterations).
param.solver = 2.0;		%Type of solver used:
					%1 = Normal point-wise Gauss-Seidel relaxation.
					%2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.
param.scl = 0.75;		%Downscale to this size
param.scl_factor = 0.75;	%Downscale factor

param = setParameters(param, varargin{:});

%Current size and downscale to size
[rows cols frames] = size(I_in);
ds_rows = ceil(rows*param.scl);
ds_cols = ceil(cols*param.scl);
scales = double(intmax);

%Smoothening mask
%G = fspecial('gaussian',[7 7],2.0);
G = fspecial('gaussian',[5 5],1.25);
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
		Itin{scl} = imfilter( Iin{scl}, G, 'replicate' );
		break
	end
end
Iout = Iin{scales};

for scl=scales:-1:1
	for iter=0:param.outer_iter
	
		[wW wNW wN wNE wE wSE wS wSW] = ADdiffWeights(Iout);

		PsiData = 1./sqrt( (Iout-Iin{scl}).^2 + eps );
		TRACE = PsiData + param.alpha*(wW+wNW+wN+wNE+wE+wSE+wS+wSW);
		B = PsiData.*Iin{scl};

		[Iout] =	PDEsolver8(	single(Iout), ...
						single(TRACE), ...
						single(B), ...
						single(param.alpha*wW), ...
						single(param.alpha*wNW), ...
						single(param.alpha*wN), ...
						single(param.alpha*wNE), ...
						single(param.alpha*wE), ...
						single(param.alpha*wSE), ...
						single(param.alpha*wS), ...
						single(param.alpha*wSW), ...
						single(param.inner_iter), ...		%solver iteration
						single(param.omega), ...		%omega
						single(param.solver) );			%solver

		%Idisp=Iout;Idisp(Idisp>1)=1;Idisp(Idisp<0)=0;
		%Idisp = Iout;
		%imagesc(Idisp),colormap(gray(256)),drawnow
		
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

function [W NW N NE E SE S SW] = ADdiffWeights(D,varargin)
%function [W NW N NE E SE S SW]=ADdiffWeights(D,varargin)
%
%%Calculates anisotropic diffusion weights between pixels.
%
%Diffusion tensor is based on the following equation:
%T(dI) = 1/(||dI||^2+2*lambda^2)    *   [ ( (dI/dy)^2 + lambda^2 )  -( dI/dx*dI/dy );
%                                        -( dI/dx*dI/dy )            ( (dI/dx)^2 + lambda^2 )]
%where the output is marked as:
%T(dI) = [dyy -dxy;-dxy dxx]
%
%INPUT
%
%OUTPUT
%
%PARAMETERS
%
%

%------------------------
% Set and get parameters
%------------------------
param.Ddx = [];
param.Ddy = [];
param.lambda = -1;
param.derivator='alvarez';
param.quantile = 0.5;
%param = setParameters_v10( param,varargin{:} );

%--------------------
% Derivation kernels 
%--------------------
if strcmp(param.derivator,'sobel')
	%Sobel derivative operators
	O_dx = [1 0 -1;2 0 -2;1 0 -1]/8;
	O_dy = [1 2 1;0 0 0;-1 -2 -1]/8;
elseif strcmp(param.derivator,'alvarez')
	%Alvarez derivative operators
	O_dx =[ 1           0   -1;
		sqrt(2)     0   -sqrt(2);
		1           0   -1]./(4+sqrt(8));
	O_dy =[	1   sqrt(2)     1;
		0   0           0;
	-1   -sqrt(2)     -1]./(4+sqrt(8));
elseif strcmp(param.derivator,'simoncelli')
	
else
    error('diffusionTensor: non-existent derivator operator!!')
end

D = double(D);

[rows cols frames] = size(D);

if isempty(param.Ddx)
	param.Ddx = imfilter( D, O_dx, 'replicate', 'conv' );
end
if isempty(param.Ddy)
	param.Ddy = imfilter( D, O_dy, 'replicate', 'conv' );
end

if frames>1
	%Calculate gradient norm
	[Mval Mind] = max( (param.Ddx.^2+param.Ddy.^2),[],3 );
	[X Y] = meshgrid( 1:cols, 1:rows );
	indices_from = sub2ind( [rows cols frames],Y,X,Mind );
		
	maxDdx = param.Ddx(indices_from);
	maxDdy = param.Ddy(indices_from);
else
	maxDdx = param.Ddx;
	maxDdy = param.Ddy;
end
	
normDdxDdy = maxDdx.^2 + maxDdy.^2;

%If parameter lambda is not given, choose one adaptively
if param.lambda < 0
	sorted = sort( reshape( normDdxDdy,1,numel(normDdxDdy) ) );
	sorted(sorted==0) = [];
	if ~isempty(sorted)
		param.lambda = sorted( round( numel(sorted)*param.quantile+eps ) );
 	else
		param.lambda = 1;
	end
end

%Multiplicator
multip = 1./( (normDdxDdy)+2*param.lambda );

dyy = multip.*( maxDdy.^2 + param.lambda );
dxx = multip.*( maxDdx.^2 + param.lambda );
dxy = -multip.*( maxDdx.*maxDdy );

W =	0.5*(	dyy + circshift(dyy,[0 1]) );		W(:,1,:) = 0;
NW =	0.25*(	dxy + circshift(dxy,[1 1]) );		NW(:,1,:) = 0;		NW(1,:,:) = 0;
N =	0.5*(	dxx + circshift(dxx,[1 0]) );		N(1,:,:) = 0;
NE =	-0.25*(	dxy + circshift(dxy,[1 -1]) );		NE(:,end,:) = 0;	NE(1,:,:) = 0;
E =	0.5*(	dyy + circshift(dyy,[0 -1]) );		E(:,end,:) = 0;
SE =	0.25*(	dxy + circshift(dxy,[-1 -1]) ); 	SE(:,end,:) = 0;	SE(end,:,:) = 0;
S =	0.5*(	dxx + circshift(dxx,[-1 0]) ); 		S(end,:,:) = 0;
SW =	-0.25*(	dxy + circshift(dxy,[-1 1]) );		SW(end,:,:) = 0;	SW(:,1,:) = 0;

if frames>1
    W = repmat(W, [1 1 frames]);
    NW = repmat(NW, [1 1 frames]);
    N = repmat(N, [1 1 frames]);
    NE = repmat(NE, [1 1 frames]);
    E = repmat(E, [1 1 frames]);
    SE = repmat(SE, [1 1 frames]);
    S = repmat(S, [1 1 frames]);
    SW = repmat(SW, [1 1 frames]);
end

