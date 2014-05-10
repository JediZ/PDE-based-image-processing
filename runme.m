addpath( genpath('./images') );
addpath( genpath('./matlab') );
addpath( genpath('./mex/build') );

%Fist compile all...if you have already compiled everything, just comment this line!
disp('Compiling MEX-functions...')
run('./mex/buildAll.m');

%------------
%- Tsukuba -
%------------
disp('Processing Tsukuba...')
I_l_tsukuba = imread('tsukuba_left.png');
I_r_tsukuba = imread('tsukuba_right.png');
%Typical disparity
D_tsukuba = DispEminND_llin_2D( I_l_tsukuba, I_r_tsukuba, 'grad', 'gradmag' );
figure
subplot(2,2,1), imagesc( I_l_tsukuba ), axis off, title('Tsukuba, left')
subplot(2,2,2), imagesc( I_r_tsukuba ), axis off, title('Tsukuba, right')
subplot(2,2,3), imagesc( D_tsukuba ), axis off, title('Disparity')
drawnow
%Symmetric disparity
D_tsukuba_sym = DispEminND_llin_sym_2D( I_l_tsukuba, I_r_tsukuba, 'grad', 'gradmag' );
figure
subplot(2,2,1), imagesc( I_l_tsukuba ), axis off, title('Tsukuba, left')
subplot(2,2,2), imagesc( I_r_tsukuba ), axis off, title('Tsukuba, right')
subplot(2,2,3), imagesc( D_tsukuba_sym(:,:,1) ), axis off, title('Symmetric disparity 1')
subplot(2,2,3), imagesc( D_tsukuba_sym(:,:,2) ), axis off, title('Symmetric disparity 2')
drawnow

%---------
% Urban3 -
%---------
disp('Processing Urban3...')
I_7_urban3 = imread('Urban3_frame07.png');
I_8_urban3 = imread('Urban3_frame08.png');
%Optical flow, late linearization with warping
[U_urban3 V_urban3] = FlowEminND_llin_2D_v10( cat(3,I_7_urban3,I_8_urban3), 3, 'grad', 'gradmag' );
OFC_urban3 = flow2color( cat(3, U_urban3, V_urban3 ), 'border', 10 );
figure
subplot(2,2,1), imagesc( I_7_urban3 ), axis off, title('Urban 3, frame 7')
subplot(2,2,2), imagesc( I_8_urban3 ), axis off, title('Urban 3, frame 8')
subplot(2,2,3), imagesc( sqrt(U_urban3.^2 + V_urban3.^2) ), axis off, title('Late linearization, velocity')
subplot(2,2,4), imagesc( OFC_urban3 ), axis off, title('Late linearization, color codified flow')
drawnow
%Horn&Schunck
[U_urban3_HS V_urban3_HS] = FlowEminHS_elin_2D_v10( cat(3,I_7_urban3,I_8_urban3), 3 );
OFC_urban3_HS = flow2color( cat(3, U_urban3_HS, V_urban3_HS ), 'border', 10 );
figure
subplot(2,2,1), imagesc( I_7_urban3 ), axis off, title('Urban 3, frame 7')
subplot(2,2,2), imagesc( I_8_urban3 ), axis off, title('Urban 3, frame 8')
subplot(2,2,3), imagesc( sqrt(U_urban3_HS.^2 + V_urban3_HS.^2) ), axis off, title('H&S, velocity')
subplot(2,2,4), imagesc( OFC_urban3_HS ), axis off, title('H&S, color codified flow')
drawnow

%------------
%- Beanbags -
%-------------
disp('Processing Beanbags...')
I_10_beanbags = imread('beanbags_frame10.png');
I_11_beanbags = imread('beanbags_frame11.png');
[U_beanbags V_beanbags] = FlowEminND_llin_2D_v10( cat(3,I_10_beanbags,I_11_beanbags), 3, 'rgb', 'none' );
OFC_beanbags = flow2color( cat(3, U_beanbags, V_beanbags ), 'border', 10 );
figure
subplot(2,2,1), imagesc( I_10_beanbags ), axis off, title('Beanbags, frame 10')
subplot(2,2,2), imagesc( I_11_beanbags ), axis off, title('Beanbags, frame 11')
subplot(2,2,3), imagesc( sqrt(U_beanbags.^2 + V_beanbags.^2) ), axis off, title('Late linearization, velocity')
subplot(2,2,4), imagesc( OFC_beanbags ), axis off, title('Late linearization, color codified flow')
drawnow

%----------
% Segment -
%----------
disp('Segmenting...')
%Load pre-calculated disparity maps
load('disparity_maps')
%Segment the dense disparity map
[Ad Bd Cd] = DispSegmentation(Dd);

%Segment the sparse disparity map
[As Bs Cs] = DispSegmentationSparse(Ds);

%And plot the results for the user
figure
subplot(2,2,1), imagesc( Dd ), title('Dense disparity map')
subplot(2,2,2), imagesc( Ds ), title('Sparse disparity map')
subplot(2,2,3), imagesc( Bd ), title('Segmentation (dense)')
subplot(2,2,4), imagesc( Bs ), title('Segmentation (sparse)')
