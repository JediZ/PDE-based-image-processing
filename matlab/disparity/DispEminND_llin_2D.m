function [U varargout] = DispEminND_llin_2D( Il, Ir, fstTerm, sndTerm, varargin)
%function [U] = DispEminND_llin_2D(Il, Ir, fstTerm, sndTerm, varargin)
%
%DispEminX_Y_Z
%	X	=	ND/AD (nonlinear diffusion/anisotropic diffusion)
%	Y	=	elin/llin/sym (early/late linearization/symmetry imposed)
%	Z	=	2D/3D (2D or 3D smoothness constraint)
%
%EXAMPLE: DispEminND_llin_2D( I_l, I_r, 'grad', 'gradmag', 'alpha', 0.15 )
%
%Uses functions medfilt2 and imresize from Image Processing Toolbox
%
%INPUT
%Il		=		Left stereo-image.
%Ir		=		Right stereo-image.
%fstTerm	=		First constancy term: 'rgb','grad'
%sndTerm	=		Second constancy term: 'none', 'rgb', 'gradmag'
%
%OUTPUT
%U		=		Horizontal component of movement (disparity).
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

%---------------------------------------------------------------
% Get and set parameters (basic parameters work with grad+none)
%---------------------------------------------------------------
param.alpha 		= 0.042;		%Smoothness weight (diffusion coefficient)
param.gammaS 		= 0.005;		%Apriori weight, spatial.
param.omega 		= 1.9;			%Relaxation parameter for SOR
param.firstLoop 	= 4;			%Number of iterations, outer fixed-point loop
param.secondLoop 	= 6;			%Number of iterations, inner fixed-point loop
param.iter 		= 4;			%Number of SOR iterations
param.b1 		= 1.48;			%First constancy term coefficient
param.b2 		= 0.29;			%Second constancy term coefficient
param.scales = double(intmax);			%Used scales, defaults to 'realmax' but downscaling is 
						    %stopped when image-size is less than 10 pixels.
param.scl_factor 	= 0.75;			%Scale factor
param.solver 		= 2.0;			%Type of solver used:
						    %1 = Normal point-wise Gauss-Seidel relaxation.
						    %2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.
param.Us = [];					%Spatial constraint for U.
param = setParameters(param, varargin{:});

fstTerm = upper( fstTerm );
sndTerm = upper( sndTerm );

%Input is probably uint8 (0-255). For calculation we need 'single' (float in C).
%Given parameters work 'properly' for images scaled between 0-1.
Il = single(Il)./255;
Ir = single(Ir)./255;

%----------------------
% Variable definitions
%----------------------

Cu1 = [];
Du1 = [];
gD1 = [];
Cu2 = [];
Du2 = [];
gD2 = [];

ASCu = [];
ASDu = [];
gS = [];
USap = [];

U = [];

%-------------------
% Smoothening mask 
%-------------------
G = fspecial('gaussian',[5 5],1.25);

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
		%---------------------
		% Smoothen last scale
		%---------------------
		It0{scl} = imfilter( It0{scl}, G, 'replicate' );
		It1{scl} = imfilter( It1{scl}, G, 'replicate' );
		break
	end
end

%----------------------
% First constancy term
%----------------------
for scl=1:param.scales

	switch (fstTerm)

		case 'RGB' 		%RGB constancy
			I1t0{scl} = It0{scl};
			I1t1{scl} = It1{scl};
		case 'GRAD' 		%Gradient constancy
			I1t0{scl} = rgb2grad(It0{scl});
			I1t1{scl} = rgb2grad(It1{scl});
		otherwise
			error('No such fstTerm');
	end
end

%-----------------------
% Second constancy term
%-----------------------
for scl=1:param.scales

	switch (sndTerm)

		case 'NONE'		%No second term
			I2t0{scl} = [];
			I2t1{scl} = [];
		case 'RGB'		%RGB
			I2t0{scl} = It0{scl};
			I2t1{scl} = It1{scl};
		case 'GRADMAG'		%Gradient magnitude constancy (based on the first term)
			I2t0{scl} = It0{scl};
			I2t1{scl} = It1{scl};
		otherwise
			error('No such sndTerm');
	end
end

%-------------------------------------
% Scale spatial apriori disparity map
%-------------------------------------
if ~isempty(param.Us)

	param.Us( isnan(param.Us) ) = 0;
	USap{1} = param.Us;

	for scl=2:param.scales
		USap{scl} = imresize( USap{scl-1}.*param.scl_factor, param.scl_factor, 'bilinear' );
	end

	U = USap{param.scales};
end

%----------------------------
% Multiscale, coarse to fine
%----------------------------
for scl=param.scales:-1:1
	
	ASdiff = 1.75*(1/param.scl_factor)^-(scl-1);
	%ASdiff = 0.9*(1/param.scl_factor)^-(scl-1);
	
	[rows cols channels1]	= size( I1t1{scl} );
	[rows2 cols2 channels2] = size( I2t1{scl} );

	[X Y] = meshgrid(1:cols,1:rows);
	I1t1w = [];
	I2t1w = [];

	%--------------------------------
	% Initial guess is "zero flow"
	%--------------------------------
	if isempty(U)
		U = zeros( rows, cols );
	end

	%-----------------------
	% First iteration loop
	%-----------------------
	for firstLoop=1:param.firstLoop

		%-------------------------
		% Warp 1st constancy term
		%-------------------------
		I1t1w = BilinInterp_2d( single(I1t1{scl}), single(X+U), single(Y) );

		%-------------------------
		% Warp 2nd constancy term
		%-------------------------
		if( ~isempty(I2t1{scl}) )
			I2t1w = BilinInterp_2d( single(I2t1{scl}), single(X+U), single(Y) );
		end

		%------------------------------
		% Calculate 1st constancy term 
		%------------------------------
		[I1dt I1dx I1dy] = FstDerivatives5( single(I1t0{scl}), single(I1t1w) );
		Cu1 = I1dt.*I1dx;
		Du1 = I1dx.*I1dx;

		%------------------------------
		% Calculate 2nd constancy term 
		%------------------------------
		if( ~isempty(I2t1{scl}) )
			if( ~strcmp(sndTerm,'GRADMAG') )
				[I2dt I2dx I2dy] = FstDerivatives5( single(I2t0{scl}), single(I2t1w) );
				Cu2 = I2dt.*I2dx;
				Du2 = I2dx.*I2dx;
			else
				[I2dxt I2dyt I2dxx I2dyy I2dxy] = SndDerivatives5( single(I2t0{scl}), single(I2t1w) );
				Cu2 = I2dxt.*I2dxx + I2dyt.*I2dxy;
				Du2 = I2dxx.*I2dxx + I2dxy.*I2dxy;
			end
		end

		%-------------------------
		% Apriori constancy terms
		%-------------------------
		%Term from spatial apriori
		if ~isempty(USap)		
			ASCu = (USap{scl}-U);
			ASDu = 1;
		end

		dU = zeros( rows, cols );

		%-----------------------
		% Second iteration loop
		%-----------------------
		for secondLoop=1:param.secondLoop

			dU1rep = repmat(dU,[1 1 channels1]);
			dU2rep = repmat(dU,[1 1 channels2]);

			%----------------------------------
			% Non-quadratic "robust" functions
			%----------------------------------
			%First constancy term
			OPnorm1 = (I1dt - I1dx.*dU1rep ).^2;
			gD1 = param.b1./(param.alpha.*realsqrt(OPnorm1+0.00001));
			%Second constancy term
			if( ~isempty(I2t1{scl}) )
				if( ~strcmp(sndTerm,'GRADMAG') )
					OPnorm2 = (I2dt - I2dx.*dU2rep ).^2;
				else
					OPnorm2 = (I2dxt - I2dxx.*dU2rep).^2 + (I2dyt - I2dxy.*dU2rep).^2;
				end
				gD2 = param.b2./(param.alpha.*realsqrt(OPnorm2+0.00001));
			end
			%Spatial apriori constancy
			if ~isempty(USap)
				%--------------------------------------------------
				% gS => spatial apriori term: "influence function"
				%--------------------------------------------------
				APnorm = (USap{scl}-U-dU).^2;
				%gS = (param.gammaS)./(param.alpha.*(1 + APnorm/ASdiff.^2));
				gS = (param.gammaS)/(param.alpha)*exp(-APnorm/(ASdiff.^2));
			end

			%-------------------
			% Diffusion weights
			%-------------------
			[wW wN wE wS] = DdiffWeights( single(U+dU),single(0.00001) );

			CuGd = sum( cat(3, Cu1.*gD1, Cu2.*gD2, ASCu.*gS ),3 );
			DuGd = sum( cat(3, Du1.*gD1, Du2.*gD2, ASDu.*gS ),3 );

			%--------
			% Solver
			%--------
			[dU] = Disp_sor_llin4_2d(	single(U), ...
							single(dU), ...
							single(CuGd), ...
							single(DuGd), ...
							single(wW), ...
							single(wN), ...
							single(wE), ...
							single(wS), ...
							single(param.iter), ...
							single(param.omega), ...
							single(param.solver) ...
						);


		end
		%-------------------
		% End of second loop
		%-------------------
		U = medfilt2(U + dU,[3 3], 'symmetric');

	end 
	%-------------------
	% End of first loop
	%-------------------

	%--------------
	% Upscale flow
	%--------------
	if (scl-1)>0
		U = imresize( U.*(1/param.scl_factor), 'bilinear' ,'OutputSize', isizes{scl-1}(1:2)  );
	end
	
end

if nargout>1
	varargout{1} = I1t1w(:,:,1);
end

%---------------------------------------
%--- rgb2grad: rgb to gradient
%---------------------------------------
function OUT=rgb2grad( IN )

	[rows cols frames] = size(IN);

	Odx = [1 0 -1];
	Ody = Odx';

	for i=1:frames
	    OUT(:,:,i*2-1) = imfilter( IN(:,:,i), Odx, 'replicate' );
	    OUT(:,:,i*2) = imfilter( IN(:,:,i), Ody, 'replicate' );
	end
