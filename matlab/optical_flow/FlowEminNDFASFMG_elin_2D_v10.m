function [U V] = FlowEminND_elin_2D(Iin, channels, varargin)
%function [U V] = FlowEminND_elin_2D(Iin, channels, varargin)
%
%FlowEminX_Y_Z
%	X	=	ND/AD (nonlinear diffusion/anisotropic diffusion)
%	Y	=	elin/llin/sym (early/late linearization, symmetry used)
%	Z	=	2D/3D (2D or 3D smoothness constraint)
%
%Calculates optical flow based on energy minimization with "brightness" and/or gradient constancy as
%data terms and spatial regularization. Uses early linearization and several scales.
%
%Characteristics:
%-FMG = Full Multigrid (DCA = Discretization Coarse grid Approximation)
%-FAS = Full Approximation Scheme
%-Early linearization = Small displacements
%
%INPUT
%Iin		=		Image secuence.
%channels	=		Number of channels per image (1 for intensity, 3 for RGB e.g.)
%
%OUTPUT
%U		=		Horizontal component of movement.
%V		=		Vertical component of movement.
%
%PARAMETERS
%
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

%------------------------
% Get and set parameters
%------------------------
%FOR YOSEMITE
param.alpha = 0.035;			%Smoothness factor (diffusion coefficient)
param.omega = 1.9;			%Relaxation parameter for SOR
param.firstLoop = 4;			%Number of iterations, outer fixed-point loop
param.iter = 4;				%Number of SOR iterations
param.b1 = 0.03;			%Brightness constancy coefficient
param.b2 = 0.97;			%Gradient constancy coefficient
param.scl_factor = 0.5;			%Scale factor
param.solver = 2.0;			%Type of solver used:
					%	1 = Normal point-wise Gauss-Seidel relaxation.
					%	2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.
param.cycle_index = 1;			%Defines the type of multi-grid cycle to use used:
					%	1 = Vcycle.
					%	2 = Wcycle.
param.scales = intmax;

param = setParameters(param, varargin{:});

%Input is probably uint8 (0-255). For calculation we need 'single' (float in C).
%Parameters in this function have been set so that input images are expected to be scaled between 0-255.
Iin = single( Iin );

tic
%----------------------
% Normalized RGB space
%----------------------
%Iin = double(Iin);
%Iin = Iin/max(max(max(Iin)));

%--------------------
% Derivation kernels 
%--------------------
%Spatial and temporal prefilters
prefilter_spa = [0.037659 0.249724 0.439911 0.249724 0.037659];
%Derivation kernels
O_dx = [0.104550 0.292315 0.0 -0.292315 -0.104550];
O_dx_scl = [0.104550 0.292315 0.0 -0.292315 -0.104550]/255;
O_dy = O_dx';
O_dy_scl = O_dx_scl';
O_dxx = [0.232905 0.002668 -0.471147 0.002668 0.232905];
O_dyy = O_dxx';

%------------------------
%--- Smoothening mask ---
%------------------------
G = fspecial('gaussian',[5 5],1);
lpf = [1 4 6 4 1]/16;

%--------------------
% Scale input images 
%--------------------
It0{1} = imfilter( Iin(:,:,1:channels),G,'replicate','conv' );
It1{1} = imfilter( Iin(:,:,channels+1:2*channels),G,'replicate','conv' );
isizes{1} = size(It0{1});
for scl=2:param.scales
	It0{scl} = imfilter( imfilter(It0{scl-1},lpf,'replicate','conv'),lpf','replicate','conv' );
	It1{scl} = imfilter( imfilter(It1{scl-1},lpf,'replicate','conv'),lpf','replicate','conv' );
	It0{scl} = It0{scl}(1:2:end,1:2:end,:);
	It1{scl} = It1{scl}(1:2:end,1:2:end,:);
	isizes{scl} = size(It0{scl});
	
	%If scaled image size is less than 10, break
	if isizes{scl}(1)<=10 | isizes{scl}(2)<=10
		param.scales = scl;
		break
	end
end

%---------------------------
% Calculate image derivates
%---------------------------
for scl=1:param.scales
	
	Ist = (It0{scl}+It1{scl}).*0.55/255;
	Idt{scl} = (It0{scl}-It1{scl})/255;
	
	Idx{scl} = imfilter( imfilter(Ist,prefilter_spa','replicate','conv'), O_dx,'replicate','conv' );
	Idy{scl} = imfilter( imfilter(Ist,prefilter_spa,'replicate','conv'), O_dy,'replicate','conv' );
	Idxx{scl} = imfilter( imfilter(Ist,prefilter_spa','replicate','conv'), O_dxx,'replicate','conv' );
	Idyy{scl} = imfilter( imfilter(Ist,prefilter_spa,'replicate','conv'), O_dyy,'replicate','conv' );
	Idxy{scl} = imfilter( imfilter(Ist,O_dx,'replicate','conv'), O_dy,'replicate','conv' );
	
	Idxt0 = imfilter( imfilter(It0{scl},prefilter_spa','replicate','conv'), O_dx_scl,'replicate','conv' );
	Idxt1 = imfilter( imfilter(It1{scl},prefilter_spa','replicate','conv'), O_dx_scl,'replicate','conv' );
	Idxt{scl} = (Idxt0-Idxt1);
	
	Idyt0 = imfilter( imfilter(It0{scl},prefilter_spa,'replicate','conv'), O_dy_scl,'replicate','conv' );
	Idyt1 = imfilter( imfilter(It1{scl},prefilter_spa,'replicate','conv'), O_dy_scl,'replicate','conv' );
	Idyt{scl} = (Idyt0-Idyt1);

	%------------------------------------------------------
	% Calculate terms that are constant in the second loop
	%------------------------------------------------------
	M{scl} = param.b1*Idy{scl}.*Idx{scl} + param.b2*Idxy{scl}.*( Idxx{scl} + Idyy{scl} );
	Cu{scl} = param.b1*Idt{scl}.*Idx{scl} + param.b2*(Idxt{scl}.*Idxx{scl} + Idyt{scl}.*Idxy{scl});
	Cv{scl} = param.b1*Idt{scl}.*Idy{scl} + param.b2*(Idxt{scl}.*Idxy{scl} + Idyt{scl}.*Idyy{scl});
	Du{scl} = param.b1*Idx{scl}.*Idx{scl} + param.b2*(Idxx{scl}.*Idxx{scl} + Idxy{scl}.*Idxy{scl});
	Dv{scl} = param.b1*Idy{scl}.*Idy{scl} + param.b2*(Idxy{scl}.*Idxy{scl} + Idyy{scl}.*Idyy{scl});
end

U = [];
V = [];

%----------------------------
% Multiscale, coarse to fine
%----------------------------
for scl=param.scales:-1:1

	[rows cols frames] = size( It0{scl} );
	
	%--------------------------------
	% Initial guess is "zero flow"
	%--------------------------------
	if isempty(U)
		U = zeros( rows, cols );
		V = zeros( rows, cols );
	end

	%----------------------------------------------------------
	% If the parameter cycle_index is 1 FAS_VCYCLE can be used 
	%----------------------------------------------------------
	%[U V] = FAS_VCYCLE(U, V, Idx, Idy, Idt, Idxx,	Idyy, Idxy, Idxt, Idyt, M, Cu, Cv, Du, Dv, scl, param.scales, param );
	[U V] = FAS_CYCLE(U, V, Idx, Idy, Idt, Idxx, Idyy, Idxy, Idxt, Idyt, M, Cu{scl}, Cv{scl}, Du, Dv, scl, param.scales, param );

	%--------------
	% Upscale flow
	%--------------
	if (scl-1)>0
		U = imresize( U.*(1/param.scl_factor), isizes{scl-1}(1:2) );
		V = imresize( V.*(1/param.scl_factor), isizes{scl-1}(1:2) );
	end
end

toc
%----------------------------
% !!!END OF MAIN FUNCTION!!!
%----------------------------

%-------------------------------------------------------
%--- FULL-APPROXIMATION-SCHEME (FAS) MULTIGRID CYCLE ---
%-------------------------------------------------------
function [U V] = FAS_CYCLE(U, V, Idx, Idy, Idt, Idxx, Idyy, Idxy, Idxt, Idyt, M, Cu, Cv, Du, Dv, scl, scales, param )
%This function implements both V- and W-cycles. Type of cycle used is based on the parameter param.cycle_index

%--------------------------------------
% Restriction operator (for multigrid)
%--------------------------------------
fw = 1/16 *[1 2 1;2 4 2;1 2 1];		%aka full-weighting

Uc = [];
Vc = [];

if scl<scales
	for ci=1:param.cycle_index
		%Presmoothing
		[U V RU RV] = smooth(	U, V, Idx{scl}, Idy{scl}, Idt{scl}, Idxx{scl}, ...
					Idyy{scl}, Idxy{scl}, Idxt{scl}, Idyt{scl}, ...
					M{scl}, Cu, Cv, Du{scl}, Dv{scl}, param );

		%Restrict residual
		RUres = imfilter( RU*param.scl_factor, fw, 'replicate', 'conv' ); RUres=RUres(1:2:end,1:2:end,:);
		RVres = imfilter( RV*param.scl_factor, fw, 'replicate', 'conv' ); RVres=RVres(1:2:end,1:2:end,:);
		
		%Restrict solution
		Ures = imfilter( U*param.scl_factor, fw, 'replicate', 'conv' ); Ures=Ures(1:2:end,1:2:end,:);
		Vres = imfilter( V*param.scl_factor, fw, 'replicate', 'conv' ); Vres=Vres(1:2:end,1:2:end,:);
		
		%RHS (right hand size of Ax=b) => "residual" terms
		%------------------------------------------
		% Update gD => data "regularizer" function 
		%------------------------------------------
		%Brightness + gradient constancy
		[rows cols channels] = size( Idxx{scl} );
		Urep = repmat( Ures,[1 1 channels] );
		Vrep = repmat( Vres,[1 1 channels] );
	
		OPnorm =	param.b1*	(Idt{scl+1} - Idx{scl+1}.*Urep - Idy{scl+1}.*Vrep).^2 + ...
				param.b2*( 	(Idxt{scl+1} - Idxx{scl+1}.*Urep - Idxy{scl+1}.*Vrep).^2 + ...
						(Idyt{scl+1} - Idxy{scl+1}.*Urep - Idyy{scl+1}.*Vrep).^2 );
	
		gd = 1./(param.alpha*sqrt(OPnorm+0.00001));
		[wW wN wS wE] = OPdiffWeights( Ures,Vres );
	
		MGd = M{scl+1}.*gd;
		DuGd = Du{scl+1}.*gd;
		DvGd = Dv{scl+1}.*gd;
		%LHS (Left-Hand-Side) (based on restricted terms)
		[Au Av] = Oflow_lhs_elin4_2d( 		single(Ures), ...
							single(Vres), ...
							single(MGd), ...
							single(DuGd), ...
							single(DvGd), ...
							single(wW), ...
							single(wN), ...
							single(wE), ...
							single(wS) ...
					);
	
		fu = (RUres + Au)./gd;
		fv = (RVres + Av)./gd;

		[Uc Vc] = FAS_CYCLE(Ures, Vres, Idx, Idy, Idt, Idxx, Idyy, Idxy, Idxt, Idyt, M, fu, fv, Du, Dv, scl+1, scales, param );

		%New estimation
		U = U + imresize( (Uc-Ures)*(1/param.scl_factor),  size(U), 'bilinear' );
		V = V + imresize( (Vc-Vres)*(1/param.scl_factor),  size(V), 'bilinear' );
	end
else
	%Presmoothing
	[U V] = smooth(	U, V, Idx{scl}, Idy{scl}, Idt{scl}, Idxx{scl}, ...
			Idyy{scl}, Idxy{scl}, Idxt{scl}, Idyt{scl}, ...
			M{scl}, Cu, Cv, Du{scl}, Dv{scl}, param );
end

if ~isempty( Uc )

	%Postsmoothing
	[U V] = smooth(	U, V, ...
			Idx{scl}, Idy{scl}, Idt{scl}, Idxx{scl}, ...
			Idyy{scl}, Idxy{scl}, Idxt{scl}, Idyt{scl}, ...
			M{scl}, Cu, Cv, Du{scl}, Dv{scl}, param );
end

%-----------------------------------------------
%--- FULL-APPROXIMATION-SCHEME (FAS) V-CYCLE ---
%-----------------------------------------------
function [U V] = FAS_VCYCLE(U, V, Idx, Idy, Idt, Idxx, Idyy, Idxy, Idxt, Idyt, M, Cu, Cv, Du, Dv, scl, scales, param )
%This function implements a V-cycle

%--------------------------------------
% Restriction operator (for multigrid)
%--------------------------------------
fw = 1/16 *[1 2 1;2 4 2;1 2 1];		%aka full-weighting

fu{scl} = Cu{scl};
fv{scl} = Cv{scl};
Uin{scl} = U;
Vin{scl} = V;

[Uout{scl} Vout{scl} RU RV] = smooth(	Uin{scl}, Vin{scl}, Idx{scl}, Idy{scl}, Idt{scl}, Idxx{scl}, ...
					Idyy{scl}, Idxy{scl}, Idxt{scl}, Idyt{scl}, ...
					M{scl}, fu{scl}, fv{scl}, Du{scl}, Dv{scl}, param );
if scales>scl
	%Fine-to-coarse
	for scli=scl+1:scales
	%scli=scl+1;
		%Restrict residual
		RU = imfilter( RU*param.scl_factor, fw, 'replicate','conv' ); RU = RU(1:2:end,1:2:end);
		RV = imfilter( RV*param.scl_factor, fw, 'replicate','conv' ); RV = RV(1:2:end,1:2:end);
		
		%Restrict solution
		Uin{scli} = imfilter( Uout{scli-1}*param.scl_factor, fw, 'replicate','conv' ); Uin{scli} = Uin{scli}(1:2:end,1:2:end);
		Vin{scli} = imfilter( Vout{scli-1}*param.scl_factor, fw, 'replicate','conv' ); Vin{scli} = Vin{scli}(1:2:end,1:2:end);
		
		%RHS (right hand size of Ax=b) => "residual" terms
		%------------------------------------------
		% Update gD => data "regularizer" function 
		%------------------------------------------
		%Brightness + gradient constancy
		[rows cols channels] = size( Idxx{scli} );
		Urep = repmat( Uin{scli}, [1 1 channels] );
		Vrep = repmat( Vin{scli}, [1 1 channels] );
		OPnorm =	param.b1*	(Idt{scli} - Idx{scli}.*Urep - Idy{scli}.*Vrep).^2 + ...
				param.b2*( 	(Idxt{scli} - Idxx{scli}.*Urep - Idxy{scli}.*Vrep).^2 + ...
						(Idyt{scli} - Idxy{scli}.*Urep - Idyy{scli}.*Vrep).^2 );
		gd = 1./(param.alpha*sqrt(OPnorm+0.00001));
		[wW wN wS wE] = OPdiffWeights( Uin{scli},Vin{scli} );

		MGd = sum(M{scli}.*gd,3);
		DuGd = sum(Du{scli}.*gd,3);
		DvGd = sum(Dv{scli}.*gd,3);
		%LHS (Left-Hand-Side) (based on restricted terms)
		[Au Av] = Oflow_lhs_elin4_2d( 		single(Uin{scli}), ...
							single(Vin{scli}), ...
							single(MGd), ...
							single(DuGd), ...
							single(DvGd), ...
							single(wW), ...
							single(wN), ...
							single(wE), ...
							single(wS) ...
					);


		fu{scli} = (RU + Au)./gd;
		fv{scli} = (RV + Av)./gd;
		
		[Uout{scli} Vout{scli} RU RV] = smooth(	Uin{scli}, Vin{scli}, Idx{scli}, Idy{scli}, Idt{scli}, Idxx{scli}, ...
							Idyy{scli}, Idxy{scli}, Idxt{scli}, Idyt{scli}, ...
							M{scli}, fu{scli}, fv{scli}, Du{scli}, Dv{scli}, param );
	
	end
	
	%Coarse-to-fine
	for scli=scales-1:-1:scl
	%scli=scl;

		Uout{scli} = Uout{scli} + imresize( (Uout{scli+1}-Uin{scli+1})*(1/param.scl_factor), size( Uout{scli} ), 'bilinear' );
		Vout{scli} = Vout{scli} + imresize( (Vout{scli+1}-Vin{scli+1})*(1/param.scl_factor), size( Vout{scli} ), 'bilinear' );

		[Uout{scli} Vout{scli}] = smooth(	Uout{scli}, Vout{scli}, ...
							Idx{scli}, Idy{scli}, Idt{scli}, Idxx{scli}, ...
							Idyy{scli}, Idxy{scli}, Idxt{scli}, Idyt{scli}, ...
							M{scli}, fu{scli}, fv{scli}, Du{scli}, Dv{scli}, param );
	
	end
end

U = Uout{scl};
V = Vout{scl};

%--------------------------------------------------
%--- SMOOTHER 					---
%--- This is the solver, also called 'smoother' ---
%--------------------------------------------------
function [U V varargout] = smooth( U, V, Idx, Idy, Idt, Idxx, Idyy, Idxy, Idxt, Idyt, M, Cu, Cv, Du, Dv, param )
	
	[rows cols channels] = size( Idx );
	%-----------------------
	% First iteration loop
	%-----------------------
	for firstLoop=1:param.firstLoop

		Urep = repmat(U,[1 1 channels]);
		Vrep = repmat(V,[1 1 channels]);

		%------------------------------------------
		% Update gD => data "regularizer" function 
		%------------------------------------------
		%Brightness + gradient constancy
		OPnorm =	param.b1*	(Idt - Idx.*Urep - Idy.*Vrep).^2 + ...
				param.b2*( 	(Idxt - Idxx.*Urep - Idxy.*Vrep).^2 + ...
						(Idyt - Idxy.*Urep - Idyy.*Vrep).^2 );
		gd = 1./(channels*param.alpha*sqrt(OPnorm+0.00001));

		%----------------------------------------------
		% Update gR => smoothness regularizer function 
		%----------------------------------------------
		[wW wN wS wE] = OPdiffWeights(U,V);

		MGd = sum(M.*gd,3);
		CuGd = sum(Cu.*gd,3);
		CvGd = sum(Cv.*gd,3);
		DuGd = sum(Du.*gd,3);
		DvGd = sum(Dv.*gd,3);
		
		[U V] = Oflow_sor_elin4_2d(		single(U), ...
							single(V), ...
							single(MGd), ...
							single(CuGd), ...
							single(CvGd), ...
							single(DuGd), ...
							single(DvGd), ...
							single(wW), ...
							single(wN), ...
							single(wE), ...
							single(wS), ...
							single(param.iter), ...
							single(param.omega), ...
							single(param.solver) ...
						);

	end
	%-----------------------------
	% End of first iteration loop
	%-----------------------------
	
	%---------------------
	% Calculate residuals
	%---------------------
	if nargout>2
		Urep = repmat(U,[1 1 channels]);
		Vrep = repmat(V,[1 1 channels]);

		%------------------------------------------
		% Update gD => data "regularizer" function 
		%------------------------------------------
		%Brightness + gradient constancy
		OPnorm =	param.b1*	(Idt - Idx.*Urep - Idy.*Vrep).^2 + ...
				param.b2*( 	(Idxt - Idxx.*Urep - Idxy.*Vrep).^2 + ...
						(Idyt - Idxy.*Urep - Idyy.*Vrep).^2 );
		gd = 1./(param.alpha*sqrt(OPnorm+0.00001));

		%----------------------------------------------
		% Update gR => smoothness regularizer function 
		%----------------------------------------------
		[wW wN wS wE] = OPdiffWeights(U,V);
		
		MGd = M.*gd;
		CuGd = Cu.*gd;
		CvGd = Cv.*gd;
		DuGd = Du.*gd;
		DvGd = Dv.*gd;
		
		[A B RU RV] = Oflow_sor_elin4_2d(	single(U), ...
							single(V), ...
							single(MGd), ...
							single(CuGd), ...
							single(CvGd), ...
							single(DuGd), ...
							single(DvGd), ...
							single(wW), ...
							single(wN), ...
							single(wE), ...
							single(wS), ...
							single(0), ...		%0 = no iterations, calculate residuals only!
							single(param.omega), ...
							single(param.solver) ...
						);

		varargout{1} = RU;
		varargout{2} = RV;
	end

%-------------------------
%--- DIFFUSION WEIGHTS ---
%-------------------------
function [wW wN wS wE varargout] = OPdiffWeights(U,V)
%function [wW wN wS wE varargout] = OPdiffWeights(U,V)
%
%Approximates "diffusion" weights between pixels. 
%
%INPUT
%
%OUTPUT
%wW			=		Western diffusitivity weights.
%wN			=		Northern diffusitivity weights.
%wS			=		Southern diffusitivity weights.
%wE			=		Eastern diffusitivity weights.
%
%PARAMETERS
%
%
%Implementation of 6-point discretization based on PhD thesis:
%"From Pixels to Regions", Thomas Brox, 2005.
%
%Author: Jarno Ralli
%E-mail: jarno@ralli.fi
%
%The program is delivered as it is and the author accepts no responsibility what so ever of its use and/or results.
%Errors and suggestions are kindly to be communicated to the author.
% 
%(C) Copyright Jarno Ralli 2010

%Cast to double. Depending on the version of Matlab, single-float can cause problems.
U = double(U);
V = double(V);

%Spatial differences
Uver = imfilter(U,[0.25 0 -0.25]','replicate');
Vver = imfilter(V,[0.25 0 -0.25]','replicate');
Uhor = imfilter(U,[0.25 0 -0.25],'replicate');
Vhor = imfilter(V,[0.25 0 -0.25],'replicate');

wW = (circshift(U,[0 1 0])-U).^2 + (Uver+circshift(Uver,[0 1 0])).^2 + (circshift(V,[0 1 0])-V).^2 + (Vver+circshift(Vver,[0 1 0])).^2;
wE = (circshift(U,[0 -1 0])-U).^2 + (Uver+circshift(Uver,[0 -1 0])).^2 + (circshift(V,[0 -1 0])-V).^2 + (Vver+circshift(Vver,[0 -1 0])).^2;
wN = (circshift(U,[1 0 0])-U).^2 + (Uhor+circshift(Uhor,[1 0 0])).^2 + (circshift(V,[1 0 0])-V).^2 + (Vhor+circshift(Vhor,[1 0 0])).^2;
wS = (circshift(U,[-1 0 0])-U).^2 + (Uhor+circshift(Uhor,[-1 0 0])).^2 + (circshift(V,[-1 0 0])-V).^2 + (Vhor+circshift(Vhor,[-1 0 0])).^2;

wW = 1./sqrt( wW + 0.00001 );
wE = 1./sqrt( wE + 0.00001 );
wN = 1./sqrt( wN + 0.00001 );
wS = 1./sqrt( wS + 0.00001 );
