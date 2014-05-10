function img = flow2color( flow, varargin )
%function img = flow2color( flow, maxvalue, varargin )
%
%Color coding of 2D-velocity information: hue codifies direction while intensity codifies magnitude.
%
%Usage:	I = flow2color( cat(3,U,V) )
%	I = flow2color( cat(3,U,V), 'border', 20 )
%
%INPUT
%flow			=	Velocity matrix (optical flow).
%
%OUTPUT
%img			=	color codified velocity.
%
%PARAMETERS
%param.maxvalue		=	maximum magnitud
%param.border		=	frame thickness
%
%Original code by Karl Pauwels (thanx for sharing Karl!)

%Parameters
param.maxvalue = [];
param.border = 0;
param = setParameters(param, varargin{:});

if param.border>0
	[rows cols frames] = size(flow);	
	%calculate image with borders
	brows = rows+2*param.border;
	bcols = cols+2*param.border;
	[X Y] = meshgrid(1:bcols, 1:brows);
	X = (X./bcols-0.5)*10;
	Y = (Y./brows-0.5)*10;
	img_border = flow2color( cat(3, X, Y) );
end

dir = atan2(-flow(:,:,2),-flow(:,:,1));
dir(dir<0) = dir(dir<0) + 2*pi;
dir = dir ./ (2*pi);
mag = sqrt(sum(flow.^2,3));

if isempty(param.maxvalue)
	param.maxvalue = max(mag(:));
	disp(['Max = ' num2str(param.maxvalue)]);
end

mag = mag ./ param.maxvalue;
mag(mag>1) = 1;

valid = (isfinite(flow(:,:,1)) & (mag<=1));

hue = ones(size(flow(:,:,1)));
sat = zeros(size(flow(:,:,1)));
val = ones(size(flow(:,:,1)));

hue(valid) = dir(valid);
sat(valid) = 1;
val(valid) = mag(valid);

img = hsv2rgb(cat(3,hue,sat,val));

%Add directional information to image border if param.border is greater than zero 0
if param.border>0
	imgOut = img_border;
	imgOut( param.border:param.border+rows-1, param.border:param.border+cols-1, : ) = img;
	img = imgOut;
end
