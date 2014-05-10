function [U varargout] = DispEminND_llin_sym_2D(Il, Ir, varargin)
%function [U ] = DispEminND_llin_sym_2D(Il, Ir, varargin)
%
%DispEminX_Y_Z
%	X	=	ND/AD (nonlinear diffusion/anisotropic diffusion)
%	Y	=	elin/llin/sym (early/late linearization/symmetry imposed)
%	Z	=	2D/3D (2D or 3D smoothness constraint)
%
%Calculates disparity based on energy minimization with "brightness" and/or gradient constancy as
%data terms and spatial regularization. Uses late linearization (aka non-linear constancy terms) and warping
%over several scales with symmetry constraint.
%
%INPUT
%Il		=		Left stereo-image.
%Ir		=		Right stereo-image.
%
%OUTPUT
%U		=		Disparity (in this case contains two disparity maps).
%
%PARAMETERS
%See in the code
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

%------------------------
% Get and set parameters
%------------------------
param.alpha = 0.035;			%Smoothness factor (diffusion coefficient)
param.beta = 0.4;			%Symmetry factor
param.omega = 1.9;			%Relaxation parameter for SOR
param.firstLoop = 3;			%Number of iterations, outer fixed-point loop
param.secondLoop = 4;			%Number of iterations, inner fixed-point loop
param.iter = 4;				%Number of SOR iterations
param.b1 = 0.25;			%Brightness constancy coefficient
param.b2 = 0.72;			%Gradient constancy coefficient
param.scales = double(intmax);		%Used scales, defaults to 'realmax' but downscaling is 
					%stopped when image-size is less than 10 pixels.
param.scl_factor = 0.75;			%Scale factor
param.solver = 2.0;			%Type of solver used:
						%1 = Normal point-wise Gauss-Seidel relaxation.
						%2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.

param = setParameters(param, varargin{:});

%--------------------
% Derivation kernels 
%--------------------
%Spatial and temporal prefilters
prefilter_spa = [0.037659 0.249724 0.439911 0.249724 0.037659];
%Derivation kernels
O_dx = [0.104550 0.292315 0.0 -0.292315 -0.104550];
O_dy = O_dx';
O_dxx = [0.232905 0.002668 -0.471147 0.002668 0.232905];
O_dyy = O_dxx';

%------------------
% Smoothening mask 
%------------------
G = fspecial('gaussian',[3 3],1);

%--------------------
% Scale input images 
%--------------------
It0{1} = Il;
It1{1} = Ir;
isizes{1} = size(It0{1});
for scl=2:param.scales
	It0{scl} = imresize( It0{scl-1}, param.scl_factor, 'bilinear' );
	It1{scl} = imresize( It1{scl-1}, param.scl_factor, 'bilinear' );
	isizes{scl} = size(It0{scl});
	%----------
	% Smoothen
	%----------
	It0{scl-1} = imfilter( It0{scl-1}, G, 'replicate' );
	It1{scl-1} = imfilter( It1{scl-1}, G, 'replicate' );
	%If scaled image size is less than 10, break
	if isizes{scl}(1)<=10 | isizes{scl}(2)<=10
		param.scales = scl;
		break
	end
end

U = [];
iter = 1;

%----------------------------
% Multiscale, coarse to fine
%----------------------------
for scl=param.scales:-1:1
	
	[rows cols channels] = size( It1{scl} );
	[X Y] = meshgrid(1:cols,1:rows);
	It1w = [];
	It0w = [];

	%--------------------------------
	% Initial guess is "zero flow"
	%--------------------------------
	if isempty(U)
		U = zeros( rows, cols, 2 );
	end

	%----------------------------
	% Scaled symmetry difference
	%----------------------------
	srDiff = 2*(1/param.scl_factor)^-(scl-1);
	ASdiff = 2*(1/param.scl_factor)^-(scl-1); %Scaled spatial difference

	%-----------------------
	% First iteration loop
	%-----------------------
	for firstLoop=1:param.firstLoop

		%-------------
		% Warp images 
		%-------------	
		It0w = BilinInterp_2d( single(It0{scl}), single(X+U(:,:,2)), single(Y) );
		It1w = BilinInterp_2d( single(It1{scl}), single(X+U(:,:,1)), single(Y) );		

		%-----------
		% Warp flow
		%-----------
		U0w = interp2( X,Y, U(:,:,1), X+U(:,:,2), Y );
		U1w = interp2( X,Y, U(:,:,2), X+U(:,:,1), Y );

		%---------------------------
		% Calculate image derivates
		%---------------------------
		[Idt0 Idx1 Idy1]		= FstDerivatives5( single(It0{scl}), single(It1w) );
		[Idxt0 Idyt0 Idxx1 Idyy1 Idxy1]	= SndDerivatives5( single(It0{scl}), single(It1w) );
		
		[Idt1 Idx0 Idy0]		= FstDerivatives5( single(It1{scl}), single(It0w) );
		[Idxt1 Idyt1 Idxx0 Idyy0 Idxy0]	= SndDerivatives5( single(It1{scl}), single(It0w) );
		
		%----------------------------
		% Calculate flow derivatives
		%----------------------------
		Udt0 = (U(:,:,1)+U1w)*0.5;
		Udx1 = imfilter( imfilter(U1w,prefilter_spa','replicate','conv'), O_dx,'replicate','conv' );
		Udy1 = imfilter( imfilter(U1w,prefilter_spa,'replicate','conv'), O_dy,'replicate','conv' );

		Udt1 = (U(:,:,2)+U0w)*0.5;
		Udx0 = imfilter( imfilter(U0w,prefilter_spa','replicate','conv'), O_dx,'replicate','conv' );
		Udy0 = imfilter( imfilter(U0w,prefilter_spa,'replicate','conv'), O_dy,'replicate','conv' );

		%------------------------------------------------------
		% Calculate terms that are constant in the second loop
		%------------------------------------------------------
		%From data term
		CuD0 =	param.b1*Idt0.*Idx1 + param.b2*(Idxt0.*Idxx1 + Idyt0.*Idxy1);
		DuD0 =	param.b1*Idx1.*Idx1 + param.b2*(Idxx1.*Idxx1 + Idxy1.*Idxy1);
		CuD1 =	param.b1*Idt1.*Idx0 + param.b2*(Idxt1.*Idxx0 + Idyt1.*Idxy0);
		DuD1 =	param.b1*Idx0.*Idx0 + param.b2*(Idxx0.*Idxx0 + Idxy0.*Idxy0);
		%From symmetry term
		CuS0 =	Udt0.*(1+Udx1);
		DuS0 =	1 + Udx1 + Udx1 + Udx1.*Udx1;
		CuS1 =	Udt1.*(1+Udx0);
		DuS1 =	1 + Udx0 + Udx0 + Udx0.*Udx0;

		dU0 = zeros( rows, cols );
		dU1 = zeros( rows, cols );

		%-----------------------
		% Second iteration loop
		%-----------------------
		for secondLoop=1:param.secondLoop

			dU0rep = repmat(dU0,[1 1 channels]);
			dU1rep = repmat(dU1,[1 1 channels]);

			%----------------------------------
			% Non-quadratic "robust" functions
			%----------------------------------
			%Data term(s)
			OPnorm0 =	param.b1*	(Idt0 - Idx1.*dU0rep).^2 + ...
					param.b2*(	(Idxt0 - Idxx1.*dU0rep).^2 + ...
							(Idyt0 - Idxy1.*dU0rep).^2 );
			gD0 = 1./(param.alpha*sqrt(OPnorm0+0.00001));
			OPnorm1 =	param.b1*	(Idt1 - Idx0.*dU1rep).^2 + ...
					param.b2*(	(Idxt1 - Idxx0.*dU1rep).^2 + ...
							(Idyt1 - Idxy0.*dU1rep).^2 );
			gD1 = 1./(param.alpha*sqrt(OPnorm1+0.00001));
			%Symmetry term(s)
			Snorm0 = ( dU0 + Udt0 + Udx1.*dU0 ).^2;
			Snorm1 = ( dU1 + Udt1 + Udx0.*dU1 ).^2;
			
			gSYM0 = (channels*param.beta/param.alpha) ./ (1 + Snorm0/srDiff^2 );
			gSYM1 = (channels*param.beta/param.alpha) ./ (1 + Snorm1/srDiff^2 );
			%gSYM0 = (channels*param.beta/param.alpha) * exp(1-(Snorm0./5)).*(1/srDiff - Snorm0./25);
			%gSYM1 = (channels*param.beta/param.alpha) * exp(1-(Snorm1./5)).*(1/srDiff - Snorm1./25);
			%gSYM0 = (channels*param.beta/param.alpha) * exp(-(Snorm0./srDiff));
			%gSYM1 = (channels*param.beta/param.alpha) * exp(-(Snorm1./srDiff));

			%-------------------
			% Diffusion weights
			%-------------------
			[wW0 wN0 wE0 wS0] = DdiffWeights( single(U(:,:,1)+dU0),single(0.00001) );
			[wW1 wN1 wE1 wS1] = DdiffWeights( single(U(:,:,2)+dU1),single(0.00001) );

			CuG0 = sum( cat(3, gD0.*CuD0, -gSYM0.*CuS0),3 );
			DuG0 = sum( cat(3, gD0.*DuD0, gSYM0.*DuS0),3 );
			CuG1 = sum( cat(3, gD1.*CuD1, -gSYM1.*CuS1),3 );
			DuG1 = sum( cat(3, gD1.*DuD1, gSYM1.*DuS1),3 );

			[dU0 dU1] =	 Disp_sor_llin_sym4_2d(	single(U(:,:,1)), ...
								single(dU0), ...
								single(CuG0), ...
								single(DuG0), ...
								single(wW0), ...
								single(wN0), ...
								single(wE0), ...
								single(wS0), ...
								single(U(:,:,2)), ...
								single(dU1), ...
								single(CuG1), ...
								single(DuG1), ...
								single(wW1), ...
								single(wN1), ...
								single(wE1), ...
								single(wS1), ...
								single(param.iter), ...
								single(param.omega), ...
								single(param.solver) ...
						);

			iter = iter + 1;

		end
		%-------------------
		% End of second loop
		%-------------------
		U(:,:,1) = medfilt2( U(:,:,1) + dU0, [3 3], 'symmetric' );	
		U(:,:,2) = medfilt2( U(:,:,2) + dU1, [3 3], 'symmetric' );
	end 
	%-------------------
	% End of first loop
	%-------------------

	%--------------
	% Upscale flow
	%--------------
	if (scl-1)>0
		U = imresize( U.*(1/param.scl_factor), 'bilinear' ,'OutputSize', isizes{scl-1}(1:2) );
	end

end

