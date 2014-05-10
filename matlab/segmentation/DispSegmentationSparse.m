function [PHI SEG SParam] = DispSegmentation(Din, varargin)
%function [PHI SEG SParam] = DispSegmentation(Din, varargin)
%
%Segments a stereo-disparity map. Segmentation is based on Chan&Vese's Active Regions.
%Similarity measure based on distance only. Implicit formulation with AOS (Additive Operator Splitting).
%Works with 'sparse' disparity maps (NaN:s encode values without approximations).
%!! This code is still in BETA testing phase !!
%
%INPUT
%Din		=		Input disparity map.
%
%OUTPUT
%PHIout		=		Level-set functions of the segments.
%SEG		=		Segments numbered.
%SParam		=		Surface parameters.
%
%PARAMETERS
%
%If you use this code, please reference (some) of my papers http://www.jarnoralli.fi
%
%Author: Jarno Ralli
%E-mail: jarno@ralli.fi
%
%The program is delivered as it is and the author accepts no responsibility what so ever of its use and/or results.
%Errors and suggestions are kindly to be communicated to the author.
% 
%Copyright 2014, Jarno Ralli
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

param.tau = 1.0;		%Time step.
param.PHI = [];			%Level-set function as input.
param.AA = [];			%Allowed area for segmentation.
param.srem_thr = 0.002;		%Segment size remove threshold.
param.varLim = 0.7;		%Segment minimum variance limit
param.polyorder = 2;		%Order of the polynomial used for modeling the surface (1 or 2).
param.seeds = 15;		%Number of seeds.
param.scl_factor = 0.75;	%Downscale factor
param.gen_scl = 0.55;		%Downscaled size (percentage) of the orignal image when generating seeds
param.rc_scl = 0.55;		%Downscaled size (percentage) of the original image for region competition
%RANSAC related parameters
param.ransac_min_cset = 0.1;	%Minimum consensus set size.
param.ransac_max_cset = 0.7;	%Maximum consensus set size.
param.ransac_cset_cycles = 10;	%How many steps it takes to go from minimum to maximum consensus set size.
param = setParameters(param, varargin{:});

%RANSAC consensus set vector
cset_vect = [param.ransac_min_cset+(param.ransac_max_cset - param.ransac_min_cset)/param.ransac_cset_cycles * [0:param.ransac_cset_cycles] ];

PHI = param.PHI;
SEG = [];
D{1} = double(Din);
D{1} = nanmedfilt2( D{1}, [5 5] );

%-----------------------------------------------------------------------------------
%Downscale 'D' and generate the 'scales' to be used in the multipyramid calculation
%-----------------------------------------------------------------------------------
%Coarse-to-fine pyramids
seedPyramid = [ 1 ];		%Add the first scale
competitionPyramid = [ 1 ];	%Add the first scale
min_scl = min( [param.gen_scl param.rc_scl] );
%Downscale disparity map
isizes{1} = size(D{1});
scales = double(intmax);
for scl=2:scales
	%D{scl} = imresize( nanmedfilt2(D{scl-1}), param.scl_factor );
	D{scl} = nanmedfilt2(imresize( nanmedfilt2(D{scl-1}, [5 5]), param.scl_factor ), [5 5] );
	isizes{scl} = size(D{scl});
	%Add the current scale into the vector
	if isizes{scl}(1) >= isizes{1}(1)*param.gen_scl & isizes{scl}(2) >= isizes{1}(2)*param.gen_scl
	    seedPyramid = [seedPyramid scl];
	end
	%Add the current scale into the vector
	if isizes{scl}(1) >= isizes{1}(1)*param.rc_scl & isizes{scl}(2) >= isizes{1}(2)*param.rc_scl
	    competitionPyramid = [competitionPyramid scl];
	end
	%Stop downscaling
	if isizes{scl}(1)<isizes{1}(1)*min_scl | isizes{1}(2)<isizes{1}(2)*min_scl
	     break;
	end
end
%Complete the pyramids
seedPyramid = [seedPyramid seedPyramid(end):-1:1 -1];
competitionPyramid = [competitionPyramid competitionPyramid(end):-1:1 -1];

fig = figure;
%-----------------------------------------------
% If PHI is empty, start by generating the seeds
%-----------------------------------------------
if isempty(PHI)
		%1st run
		[PHI SParam] = generateSeeds(	D, ...				%Data in
						seedPyramid, ...		%Scale pyramid
						param.polyorder, ...		%Order of the polynomial used for modeling the surface
						0.7, ...			%Minimum limit of variance (segment model)
						cset_vect, ...			%RANSAC consensus set size vector
						20, ...				%iterations
						param.srem_thr, ...		%Segment remove size threshold (small segments will be removed)
						param.AA==1, ...		%Allowed segmentation area
						param.seeds ...			%Number of seeds to be generated
						);

	if param.seeds ~= 1
		
		[PHI SParam] = regionCompetition(D, ...				%Data in
						competitionPyramid, ...		%Scale pyramid
						param.polyorder, ...		%Order of the polynomial used for modeling the surface
						1.5, ...			%Minimum limit of variance (segment model)
						param.ransac_max_cset, ...	%RANSAC concensus set size
						30, ...				%Number of competition iterations per scale
						param.srem_thr, ...		%Segment remove size threshold (small segments will be removed)
						PHI, ...			%Level-set functions describing the segments
						'inverse' ...			%Competition style (see below).
						);
		%2nd run
		PHI = cat(3, PHI,	generateSeeds(	D, ...			%Data in
						competitionPyramid, ...		%Scale pyramid
						param.polyorder, ...		%Order of the polynomial used for modeling the surface
						1.2, ...			%Minimum limit of variance (segment model)
						cset_vect, ...			%RANSAC consensus set size vector
						20, ...				%iterations
						param.srem_thr, ...		%Segment remove size threshold (small segments will be removed)
						sum(PHI>0,3)==0, ...		%Allowed segmentation area
						param.seeds ...			%Number of seeds to be generated
						));
		[PHI SParam] = regionCompetition(D, ...				%Data in
						competitionPyramid, ...		%Scale pyramid
						param.polyorder, ...		%Order of the polynomial used for modeling the surface
						1.5, ...			%Minimum limit of variance (segment model)
						param.ransac_max_cset, ...	%RANSAC concensus set size
						20, ...				%Number of competition iterations per scale
						param.srem_thr, ...		%Segment remove size threshold (small segments will be removed)
						PHI, ...			%Level-set functions describing the segments
						'inverse' ...			%Competition style (see below).
						);
	end
%---------------------------------------------------------------------------------------------------
% Otherwise make the existing segments compete, generate new seeds and then make them compete again
%---------------------------------------------------------------------------------------------------
else
		%Region competition
		[PHI SParam] = regionCompetition(D, ...				%Data in
						competitionPyramid, ...		%Scale pyramid
						param.polyorder, ...		%Order of the polynomial used for modeling the surface
						1.0, ...			%Minimum limit of variance (segment model)
						param.ransac_max_cset, ...	%RANSAC concensus set size
						20, ...				%Number of competition iterations per scale
						param.srem_thr, ...		%Segment remove size threshold (small segments will be removed)
						PHI, ...			%Level-set functions describing the segments
						'inverse' ...			%Competition style (see below).
						);
		%Generate seeds
		PHI = cat(3, PHI,	generateSeeds(	D, ...			%Data in
						competitionPyramid, ...		%Scale pyramid
						param.polyorder, ...		%Order of the polynomial used for modeling the surface
						1.2, ...			%Minimum limit of variance (segment model)
						cset_vect, ...			%RANSAC consensus set size vector
						20, ...				%iterations
						param.srem_thr, ...		%Segment remove size threshold (small segments will be removed)
						sum(PHI>0,3)==0, ...		%Allowed segmentation area
						1 ...				%Number of seeds to be generated
						));
		%Region competition
		[PHI SParam] = regionCompetition(D, ...				%Data in
						competitionPyramid, ...		%Scale pyramid
						param.polyorder, ...		%Order of the polynomial used for modeling the surface
						2.0, ...			%Minimum limit of variance (segment model)
						param.ransac_max_cset, ...	%RANSAC concensus set size
						20, ...				%Number of competition iterations per scale
						param.srem_thr, ...		%Segment remove size threshold (small segments will be removed)
						PHI, ...			%Level-set functions describing the segments
						'inverse' ...			%Competition style (see below).
						);
%------------
%...and end
%------------
end
close(fig);

%-------------------------------------
% Return individual segments numbered
%-------------------------------------
[rows cols segments] = size(PHI);
H1 = PHI>0.0;
SEG = H1(:,:,1);
for i=2:segments
	SEG = SEG +  H1(:,:,i)*i;
end
H1sum = sum(H1,3);
SEG(H1sum>=2) = 0;
SEG(SEG>segments) = segments+1;

%---------------------
%--- GENERATESEEDS ---
%---------------------
function [PHIout SParam] = generateSeeds( Din, pyramid, polyorder, sigmaLim, ransac_cset_vect, iterations, srem_thr, AAin, seeds )
%INPUT
%Din			=		Disparity map.
%pyramid		=		Pyramid scales.
%polyorder		=		Order of the polynomial used for modeling the surface.
%sigmaLim		=		Minimum limit for standard deviation (Bayes probability).
%ransac_cset_vect	=		RANSAC consensus set size vector.
%iterations		=		Number of iterations.
%srem_thr		=		Segment size remove threshold.
%AAin			=		Allowed area for segmentation.
%seeds			=		Number of seeds to be generated.
%
%OUTPUT
%PHI			=		Level-set function(s), one per segment.
%SParam			=		Surface parameters.

PHIout = [];
SParam = [];

gamma = 0.005;

%Derivator operators
O_dx = [-1 0 1]*0.5;
O_dy = O_dx';

%Area to be segmented (allowed area)
if isempty(AAin)
	AA{1} = ones( size(Din{1}) );
else
	AA{1} = double(AAin);
end
%-------------------
%Level-set function
%-------------------
%Initial level-set containing the 'seeds'
PHIinitial = ones( size(Din{1}) )*-1;
PHIinitial(2:5:end-1,2:5:end-1) = 1;
%Reserve space
PHI{1} = PHIinitial;
for i=2:max(pyramid)
	PHI{i} = ones( size(Din{i}) )*-1;
end

%Signal for emptysegment
SIG_emptysegment = false;

%---------
%1st loop 
%---------
for seed=1:seeds

	%-----------------------------
	%Area allowed to be segmented
	%-----------------------------
	if SIG_emptysegment == false
		%Downscale
		for i=2:max(pyramid)
			AA{i} = imresize( AA{i-1}, size(Din{i}) );
		end
	end
	
	%Reset emptysegment 
	SIG_emptysegment = false;
	minCOV = [sigmaLim];
	
	%---------
	%2nd loop
	%---------
	for cscl=1:numel(pyramid)-1
	
		scl = pyramid(cscl);
		[rows cols frames] = size( Din{scl} );
		gammaScl = gamma*(rows*cols)^(0.7);
		includeFilter = AA{scl}>0.05;
		if cscl==1 PHI{scl}(~includeFilter) = -1; end
		[X Y] = meshgrid(1:cols, 1:rows); X2 = X.^2; Y2 = Y.^2; XY = X.*Y;
		H1eq = [];
		DinNoNaN = Din{scl}; DinNoNaN( isnan( DinNoNaN ) ) = 1000;

		%Once per segment do a 'sanity' check...only the biggest connected element survives
		if cscl == round( numel(pyramid)/2 )
		    %---BWLABEL BASED---%
		    [CC num] = bwlabel( PHI{scl}>0 );
		    CCprop = regionprops( CC, 'Area' );
		    %Find the cell with most elements
		    [CCY CCI] = max([CCprop(:).Area]);
		    Ot = ones(rows, cols)*-5;
		    Ot( CC==CCI ) = 5;
		    PHI{scl} = Ot;
%		    %---BWCONNCOMP BASED---%
%		    CC = bwconncomp( PHI{scl}>0 );
%		    %Find the cell with most elements
%		    [maxY maxI] = max( cellfun(@length, CC.PixelIdxList) );
%		    Ot = ones(rows, cols)*-5;
%		    Ot( CC.PixelIdxList{maxI} ) = 5;
%		    PHI{scl} = Ot;
		end
		
		%---------
		%3rd loop
		%---------
		for iter=1:iterations
			%------------------
			%RANSAC PARAMETERS
			%------------------
			%Iterations
			if iter<=1 & cscl==1
				RITER = 2000;
			else
				RITER = 100;
			end
			%Consensus set size (function of 'iter' in the first scale)
			if cscl==1
				if iter<=length(ransac_cset_vect)
				    RCONS = ransac_cset_vect(iter);
				else
				    RCONS = ransac_cset_vect(end);
				end

			else
				RCONS = ransac_cset_vect(end);
			end
			
			%---------------------------------
			%SURFACE PARAMETERS & LIKELIHOOD
			%---------------------------------
			%Segment defined by H(PHI) => Heaviside function
			H1 = (PHI{scl}>=0.0); Xs = X(H1); Ys = Y(H1);
 
			%Empty segment
			if numel(Xs)<20
				SIG_emptysegment = true;
				break;
			end
			%DISPARITY SIMILARITY MEASURE
			%Surface parameters
			switch polyorder
				case 1
					%H1eq = LSplaneEquation(Xs, Ys, Din{scl}(H1), 0.7,  'iterations', RITER, 'bestModel', H1eq, 'c_set', RCONS );
					[H1eq Err] = SurfaceEquation( single( [Xs Ys ones(numel(Xs),1)] ), ...
								      single(DinNoNaN(H1)), ...
								      single(H1eq), ...
								      single(0.7),  ...
								      single(RCONS), ...
								      single(RITER) );

					%Distance to the fitted surface
					distD = ( H1eq(1)*X + H1eq(2)*Y + H1eq(3) - DinNoNaN ).^2;
				case 2
					%H1eq = LSquadraticEquation(Xs, Ys, Din{scl}(H1), 0.7,  'iterations', RITER, 'bestModel', H1eq, 'c_set', RCONS );
					[H1eq Err] = SurfaceEquation( single( [Xs.^2 Ys.^2 Xs.*Ys Xs Ys ones(numel(Xs),1)] ), ...
								      single(DinNoNaN(H1)), ...
								      single(H1eq), ...
								      single(0.7),  ...
								      single(RCONS), ...
								      single(RITER) );

					%Distance to the fitted surface
					distD = ( H1eq(1)*X2 + H1eq(2)*Y2 + H1eq(3)*XY + H1eq(4)*X + H1eq(5)*Y + H1eq(6) - DinNoNaN ).^2;
				otherwise
				  error('Only polynomials of 1st and 2nd order have been implemented')
			end
			
			temp = distD( H1 );
			temp = temp( temp<100 );
			covect = mean( temp );
			covect( covect<minCOV ) = minCOV( covect<minCOV );

			%Likelihood of belonging to the surface
			P1 = ( 1./realsqrt(2*pi*covect)*exp( -distD./(2*covect) ) );
			%Likelihood of not belonging to the surface
			P0 = 1./realsqrt(2*pi*covect)-P1;
			%Data that 'drives' the level-set function
			DATA = reallog( (P1+eps)./(P0+eps) );
			DATA(~includeFilter) = -2;

			%----------------------------
			%PHI => Level-set function
			%----------------------------
			%Gradient of PHI

			PHIdx = imfilter( PHI{scl}, O_dx, 'replicate' );
			PHIdy = imfilter( PHI{scl}, O_dy, 'replicate' );
			gradPHI = realsqrt(PHIdx.^2+PHIdy.^2);

			%Approximation of derivative of Heaviside function
			%Heaviside function H(z) = \frac{1}{2} \left( 1 + \frac{2}{\pi}atan( \frac{z}{\epsilon} ) \right)
			%Heaviside deltafuntion (DH) \frac{\partial H(z)}{\partial z} = \dfrac{1}{\epsilon \pi (1 + (\frac{z}{\epsilon})^2 )}
			%DH = 1./(0.5*pi*(1+ (PHI{scl}/0.5).^2) );	%e = 0.5
			DH = 1./(pi*(1+(PHI{scl}.^2)./1));		%e = 1.0

			PHI{scl} = CV_solver_2d(	single( PHI{scl} ), ...
							single( DATA ), ...
							single( DH ), ...
							single( gradPHI ), ...
							single( 1.0 ), ...		%Time step
							single( gammaScl ) );
			imagesc(PHI{scl}>0),drawnow
			
		end
		%-------------
		%End 3rd loop
		%-------------
		%Auto-adjust gamma
		if SIG_emptysegment==true
			gamma = gamma*0.8;
			break;
		end

		%Auto-adjust the minimum variance
		if cscl==round( numel(pyramid)/2 )
			temp = distD( H1 );
			temp = temp( temp<100 );
			covect = mean( temp );
			if covect>0.5
			    minCOV = covect;
			end
		end

		if pyramid(cscl+1)~=-1
			PHI{pyramid(cscl+1)} = imresize( PHI{pyramid(cscl)}, size(PHI{pyramid(cscl+1)}) );
		end

	end

	%-------------
	%End 2nd loop
	%-------------
	if SIG_emptysegment~=true
		PHIout = cat(3,PHIout,PHI{1});
		SParam = cat(2,SParam, H1eq);
		AA{1} = double(any(PHIout(:,:,end)<0,3)&AA{1});
	end

	PHI{1} = PHIinitial;

end
%-------------
%End 1st loop
%-------------

%Remove segments below the segment size threshold
[srows scols segments] = size(PHI{1});
H1 = PHI{1}>=0.0;
H1size = sum(sum(H1));
H1empty = find( H1size<srem_thr*srows*scols );
if ~isempty( H1empty )
	%Remove empty/small sets
	PHI{1}(:,:,H1empty) = [];
end

%-------------------------
%--- REGIONCOMPETITION ---
%-------------------------
function [PHIout SParam] = regionCompetition( Din, pyramid, polyorder, sigmaLim, ransac_cset, iterations, srem_thr, PHIin, competition );
%INPUT
%Din			=		Disparity map.
%pyramid		=		Pyramid scales.
%polyorder		=		Order of the polynomial used for modeling the surface.
%sigmaLim		=		Minimum limit for standard deviation (Bayes probability).
%ransac_cset_vect	=		RANSAC consensus set size vector.
%iterations		=		Number of iterations.
%srem_thr		=		Segment size remove threshold.
%PHIin			=		Level-set functions containing the segments.
%competition		=		Competition style (see below).
%
%OUTPUT
%PHIout			=		Levet-set functions after region competition.
%SParam			=		Surface parameters.

%Derivator operators
O_dx = [-1 0 1]*0.5;
O_dy = O_dx';
%Recalculate signal
SIG_recalc = false;

PHI{1} = PHIin;
for i=2:max(pyramid)
	PHI{i} = imresize( PHI{i-1}, size(Din{i}) );
end

minCOV = [sigmaLim];

for cscl=1:numel(pyramid)-1

	scl = pyramid(cscl);
	[irows icols channels] = size( Din{scl} );
	[srows scols segments] = size( PHI{scl} );
	gamma = 0.005*(irows*icols)^(0.7);
	[X Y] = meshgrid(1:scols, 1:srows); XX = X.^2; YY = Y.^2; XY = X.*Y;
	DATA = zeros(size(PHI{scl}));
	P = zeros(size(PHI{scl}));
	selectWC = 1:segments;
	DinNoNaN = Din{scl}; DinNoNaN( isnan( DinNoNaN ) ) = 1000;

	switch polyorder
	    case 1
		surface1 = zeros(3,segments);
	    case 2
		surface1 = zeros(6,segments);
	    otherwise
		error('Only polynomials of 1st and 2nd order have been implemented')
	end

	for iter=1:iterations

		%----------------------------------------------
		%AREA (defined by H(PHI) => Heaviside function)
		%----------------------------------------------
		%H(1) (inside segment)
		H1 = PHI{scl}>=0.0;
		%Check for 'small' segments
		H1size = sum(sum(H1));
		H1empty = find( H1size<srem_thr*srows*scols );
		%Remove small segments and do house holding 
		if ~isempty( H1empty )
			%Remove empty/small sets
			PHI{scl}(:,:,H1empty) = [];
			H1 = PHI{scl}>=0.0;

			%New size & reinitialize matrices
			[srows scols segments] = size( PHI{scl} );
			DATA = zeros(size(PHI{scl}));
			P = zeros(size(PHI{scl}));
			selectWC = 1:segments;
			switch polyorder
			    case 1
				surface1 = zeros(3,segments);
			    case 2
				surface1 = zeros(6,segments);
			    otherwise
				error('Only polynomials of 1st and 2nd order have been implemented')
			end
			
			%Recalculate all
			SIG_recalc = true;

			if segments<1
			    disp('no more segments left!')
			    keyboard
			end
		end

		if mod(iter, 2) || SIG_recalc==true
			%Heaviside function H(z) = \frac{1}{2} \left( 1 + \frac{2}{\pi}atan( \frac{z}{\epsilon} ) \right)
			%Heaviside deltafuntion (DH) \frac{\partial H(z)}{\partial z} = \dfrac{1}{\epsilon \pi (1 + (\frac{z}{\epsilon})^2 )}
			%DH = 1./(0.5*pi*(1+ (PHI{scl}/0.5).^2) );	%e = 0.5
			%DH = 1./(pi*(1+(PHI{scl}.^2)./1));		%e = 1.0
			DH = 1./(pi*(2+(PHI{scl}.^2)./4));		%e = 2.0
			DH( DH<0.04 ) = 0.04;				%To keep things 'active'
			%Calculate derivatives of the level-set functions
			PHIdx = imfilter( PHI{scl}, O_dx, 'replicate' );
			PHIdy = imfilter( PHI{scl}, O_dy, 'replicate' );
			gradPHI = realsqrt( PHIdx.^2 + PHIdy.^2 );
			covect = [];	
			
			%Distance to each segment
			for s=1:segments

			    switch polyorder
				case 1
					%Quadratic surface
					%surface1(:,s) = LSplaneEquation(	X(H1(:,:,s)), Y(H1(:,:,s)), Din{scl}(H1(:,:,s)), ...
					%					1.0, 'bestModel', surface1(:,s), 'c_set', ransac_cset, 'iterations', 10 );
					Xs = X(H1(:,:,s)); Ys = Y(H1(:,:,s));
					[surface1(:,s) Err] = SurfaceEquation(	single( [Xs Ys ones(numel(Xs),1)] ), ...
										single( DinNoNaN(H1(:,:,s)) ), ...
										single( surface1(:,s) ), ...
										single( 1.2 ),  ...
										single( ransac_cset ), ...
										single( 10 ) );
					%Distance from the plane
					distD = ((surface1(1,s)*X + surface1(2,s)*Y + surface1(3,s)) - DinNoNaN).^2;
				case 2
					%Quadratic surface
					%surface1(:,s) = LSquadraticEquation(	X(H1(:,:,s)), Y(H1(:,:,s)), Din{scl}(H1(:,:,s)), ...
					%					1.0, 'bestModel', surface1(:,s), 'c_set', ransac_cset, 'iterations', 10 );
					Xs = X(H1(:,:,s)); Ys = Y(H1(:,:,s));
					[surface1(:,s) Err] = SurfaceEquation(	single( [Xs.^2 Ys.^2 Xs.*Ys Xs Ys ones(numel(Xs),1)] ), ...
										single( DinNoNaN(H1(:,:,s)) ), ...
										single( surface1(:,s) ), ...
										single( 1.2 ),  ...
										single( ransac_cset ), ...
										single( 10 ) );
					%Distance from the plane
					distD = (surface1(1,s)*XX + surface1(2,s)*YY + surface1(3,s)*XY + surface1(4,s)*X + surface1(5,s)*Y + surface1(6,s) - DinNoNaN).^2;
				otherwise
					error('Only polynomials of 1st and 2nd order have been implemented')
			    end

			    %Calculate variance/covariance
			    covtemp = distD(H1(:,:,s)); covtemp = covtemp( covtemp<100 );
			    covtemp = sum(covtemp)/numel(covtemp);
			    covtemp( covtemp<minCOV ) = minCOV( covtemp<minCOV );
			    covect = cat(3,covect,covtemp);

			    %LIKELIHOOD
			    P(:,:,s) = (1./realsqrt(2*pi*covtemp)*exp( -distD./(2*covtemp) ));

			end

			WC = zeros(size(PHI{scl}));
			%REGION COMPETITION STRATEGY	
			switch competition
				case {'surface'}
				%All the segments compete equally for all the pixels. Can generate non-segmented regions.
					for s=1:segments
						%Worst competitors probability
						WC(:,:,s) = max( P(:,:,selectWC(selectWC~=s)),[],3 );
					end
				case {'greedy'}
				%Greedily segments all the pixels. Typically segments the whole domain.
					Hnotany = ~any(H1,3);
					for s=1:segments
						WCtemp = max( P(:,:,selectWC(selectWC~=s)),[],3 );
						%Empty 'segment' near to current segments (s) border => set competing probability to zero
						WCtemp( Hnotany & DH(:,:,s)>0.02 ) = 0;
						%Worst competitors probability
						WC(:,:,s) = WCtemp;
					end
				case {'inverse'}
				%Uses 'inverse likelyhood' based on the segments surface likehood. If a segment intrudes on other segments
				%region, competing likehood (for the segment being intruded) is maximum of 'inverse likelyhood' or the 
				%intruding segments likelyhood. Can generate non-segmented regions.
					Ptemp = P;Ptemp(~H1) = 0;
					for s=1:segments
						%inverse likelyhood vs. intruding segments likelyhood
						WC(:,:,s) = max( cat(3, 1./realsqrt(2*pi*covect(1,1,s)) - P(:,:,s), Ptemp(:,:,selectWC~=s)), [],3);
					end
				otherwise
					error('No such competition strategy!')
			end
			DATA = reallog( (P+eps)./(WC+eps) );
			SIG_recalc = false;
		end

		%-----------------------------------
		%Update PHI => Level-set function
		%-----------------------------------
		PHI{scl} = CV_solver_2d(	single( PHI{scl} ), ...
						single( DATA ), ...
						single( DH ), ...
						single( gradPHI ), ...
						single( 1.0 ), ...
						single( gamma ) );

		%Plot something for the user
		%-------------------------------
		% Individual segments numbered
		%-------------------------------
		[rows cols segments] = size(PHI{scl});
		H1 = PHI{scl}>0.0;
		Sout = H1(:,:,1);
		for s=2:segments
			Sout = Sout +  H1(:,:,s)*s;
		end
		Sout(Sout>segments) = 0;
		imagesc(Sout),title([int2str(iter),'/',int2str(iterations)]),colorbar
		drawnow
	end

	if pyramid(cscl+1)~=-1
			PHI{pyramid(cscl+1)} = imresize( PHI{pyramid(cscl)}, [size(PHI{pyramid(cscl+1)},1) size(PHI{pyramid(cscl+1)},2)] );
	end
end

PHIout = PHI{1};
SParam = surface1;

%---------------
%- NANMEDFILT2 -
%---------------
function D = nanmedfilt2( D, fsize )
%function D = nanmedfilt2( D, fsize )
%
%Applies a medianfilter (2D) to D, ignoring NaN:s

	D = colfilt( D,[3 3], 'sliding', @nanmedian );

