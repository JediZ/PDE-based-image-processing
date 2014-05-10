function [U V]=FlowEminHS_elin_2D(Iin, channels, varargin)
%function [U V]=FlowEminHS_elin_2D(Iin, channels, varargin)
%
%FlowEminX_Y_Z
%	X	=	HS/ND/AD (Horn&Schunk/nonlinear diffusion/anisotropic diffusion)
%	Y	=	elin/llin/ (early/late linearization) sym (symmetrical)
%	Z	=	2D/3D (2D or 3D smoothness constraint)
%
%Calculates optical flow based on energy minimization with "brightness" and/or gradient constancy as
%data terms and spatial regularization. Uses early linearization and several scales. Horn&Schunk "version".
%
%Uses functions medfilt2 and imresize from Matlab Image Processing Toolbox
%
%INPUT
%Iin		=		Image secuence (concatenated i.e. cat(3,It0, It1))
%Channels	=		Number of channels per image (1 for intensity, 3 for RGB )
%scales		=		Number of used scales
%
%OUTPUT
%U		=		Horizontal component of movement
%V		=		Vertical component of movement
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
param.alpha = 0.2;			%Smoothness factor (diffusion coefficient)
param.omega = 1.9;			%Relaxation parameter for SOR
param.iter = 20;			%Number of SOR iterations (solver)
param.b1 = 0.25;			%Brightness constancy coefficient
param.b2 = 0.75;			%Gradient constancy coefficient
param.scl_factor = 0.75;		%Scale factor
param.solver = 2.0;			%Type of solver used:
						%1 = Normal point-wise Gauss-Seidel relaxation.
						%2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.
param.scales = double(intmax);
param = setParameters(param, varargin{:});

%Input is probably uint8 (0-255). For calculation we need 'single' (float in C).
%Given parameters work 'properly' for images scaled between 0-1.
Iin = single( Iin )./255;

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

%------------------------
%--- Smoothening mask ---
%------------------------
G = fspecial('gaussian',[5 5],1.25);

%--------------------
% Scale input images 
%--------------------
It0{1} = Iin(:,:,1:channels);
It1{1} = Iin(:,:,channels+1:2*channels);
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
	if isizes{scl}(1)<=20 | isizes{scl}(2)<=20
		param.scales = scl;
		%---------------------
		% Smoothen last scale
		%---------------------
		It0{scl} = imfilter( It0{scl}, G, 'replicate' );
		It1{scl} = imfilter( It1{scl}, G, 'replicate' );
		break
	end
end

U = [];
V = [];

%----------------------------
% Multiscale, coarse to fine
%----------------------------
for scl=param.scales:-1:1

	[rows cols frames] = size( It0{scl} );
	W = param.alpha*channels*ones( rows, cols );
	
	%--------------------------------
	% Initial guess is "zero flow"
	%--------------------------------
	if isempty(U)
		U = zeros( rows, cols );
		V = zeros( rows, cols );
	end

	%---------------------------
	% Calculate image derivates
	%---------------------------
	Ist = (It0{scl}+It1{scl}).*0.55;
	Idt = (It0{scl}-It1{scl}); 
	
	Idx = imfilter( imfilter(Ist,prefilter_spa','replicate','conv'), O_dx,'replicate','conv' );
	Idy = imfilter( imfilter(Ist,prefilter_spa,'replicate','conv'), O_dy,'replicate','conv' );
	Idxx = imfilter( imfilter(Ist,prefilter_spa','replicate','conv'), O_dxx,'replicate','conv' );
	Idyy = imfilter( imfilter(Ist,prefilter_spa,'replicate','conv'), O_dyy,'replicate','conv' );
	Idxy = imfilter( imfilter(Ist,O_dx,'replicate','conv'), O_dy,'replicate','conv' );
	
	Idxt0 = imfilter( imfilter(It0{scl},prefilter_spa','replicate','conv'), O_dx,'replicate','conv' );
	Idxt1 = imfilter( imfilter(It1{scl},prefilter_spa','replicate','conv'), O_dx,'replicate','conv' );
	Idxt = (Idxt0-Idxt1);
	
	Idyt0 = imfilter( imfilter(It0{scl},prefilter_spa,'replicate','conv'), O_dy,'replicate','conv' );
	Idyt1 = imfilter( imfilter(It1{scl},prefilter_spa,'replicate','conv'), O_dy,'replicate','conv' );
	Idyt = (Idyt0-Idyt1);

	%------------------------------------------------------
	% Calculate terms that are constant in the second loop
	%------------------------------------------------------
	M = param.b1*Idy.*Idx + param.b2*Idxy.*( Idxx + Idyy );
	Cu = param.b1*Idt.*Idx + param.b2*(Idxt.*Idxx + Idyt.*Idxy);
	Cv = param.b1*Idt.*Idy + param.b2*(Idxt.*Idxy + Idyt.*Idyy);
	Du = param.b1*Idx.*Idx + param.b2*(Idxx.*Idxx + Idxy.*Idxy);
	Dv = param.b1*Idy.*Idy + param.b2*(Idxy.*Idxy + Idyy.*Idyy);

	MGd = sum(M,3);
	CuGd = sum(Cu,3);
	CvGd = sum(Cv,3);
	DuGd = sum(Du,3);
	DvGd = sum(Dv,3);

	%--------
	% Solver
	%--------
	[U V] = Oflow_sor_elin4_2d(	single(U), ...
					single(V), ...
					single(MGd), ...
					single(CuGd), ...
					single(CvGd), ...
					single(DuGd), ...
					single(DvGd), ...
					single(W), ...
					single(W), ...
					single(W), ...
					single(W), ...
					single(param.iter), ...
					single(param.omega), ...
					single(param.solver) ...
				);

	%--------------
	% Upscale flow
	%--------------
	if (scl-1)>0
		U = imresize( medfilt2(U.*(1/param.scl_factor), [3 3], 'symmetric' ) ,'OutputSize', isizes{scl-1}(1:2) );
		V = imresize( medfilt2(V.*(1/param.scl_factor), [3 3], 'symmetric' ) ,'OutputSize', isizes{scl-1}(1:2) );
	end
	
end
