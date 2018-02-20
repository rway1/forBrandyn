% These PSFs are created using Seidel Formulism.

tic;

%% user defined values which can define the system and sampling

ApertureDiameter = 10; % mm, diameter of the pupil

ApertureSampleDimension = 0.1; % mm

Fnum = 2;

lambda = 0.55; % wavelength in um

PixelPitch = 1; % 1 um

DiffLimSpotSize = 1.22*lambda*Fnum; % valid of circular aperture only

ZeroPaddingFactor = 30; % defines the amount of padding of zeros to be placed around aperture function

% Zeropaddingfactor = 30 gives very reasonable PSF results

% ApertureSampleDimension = TBD for reasonable MTF results

WavesOfDefocus = 0.25;

WavesOfTiltX = 0;

WavesOfTiltY = 0;

WavesOfSpherical = 0;

WavesOfComaX = 0;

WavesOfComaY = 0;

WavesOfAstigmatismX = 0;

WavesOfAstigmatismY = 0;

WavesOfFieldCurvatureX = 0;

WavesOfFieldCurvatureY = 0;

WavesOfDistortionX = 0;

WavesOfDistortionY = 0;

X_FieldExtent = 1; % half field of view in x direction

Y_FieldExtent = 1; % half field of view in y direction

% these values calculated from above numbers not to be edited by user

FocalLength = Fnum.*ApertureDiameter;

PupilSampling = ApertureDiameter./ApertureSampleDimension; % number sample points across aperture

NormalizedFieldX = X_FieldExtent/X_FieldExtent;

NormalizedFieldY = Y_FieldExtent/Y_FieldExtent;

disp('System Parameters Initialized')

%% creates array for pupil function which is padded on all sides

ZeroPaddingSize = (ZeroPaddingFactor.*PupilSampling);

PupilFunctionArray = zeros(PupilSampling+(2*ZeroPaddingSize),PupilSampling+(2*ZeroPaddingSize));

NormalizedPupilFunctionArray = zeros(PupilSampling+(2*ZeroPaddingSize),PupilSampling+(2*ZeroPaddingSize));

NormalizedPupilFunctionArrayX = zeros(PupilSampling+(2*ZeroPaddingSize),PupilSampling+(2*ZeroPaddingSize));

NormalizedPupilFunctionArrayY = zeros(PupilSampling+(2*ZeroPaddingSize),PupilSampling+(2*ZeroPaddingSize));

arraysize = size(PupilFunctionArray);

arrayHmiddle = arraysize(1,2)./2;

arrayVmiddle = arraysize(1,1)./2;

%NormalizedPupilFunctionArray = zeros(

% creates circular aperture of 100% tranmission

for i = 1:max(arraysize(1,1))

for j = 1:max(arraysize(1,1))

RadiusIndex = sqrt(((i-arrayHmiddle).^2) + ((j-arrayVmiddle).^2));

if RadiusIndex <= (PupilSampling/2)

PupilFunctionArray(i,j) = 1;

else

PupilFunctionArray(i,j) = 0;

end

NormalizedPupilFunctionArray(i,j) = sqrt((((i-arrayHmiddle)/(PupilSampling/2)).^2) + (((j-arrayVmiddle)/(PupilSampling/2)).^2));

NormalizedPupilFunctionArrayX(i,j) = (j-arrayVmiddle)/(PupilSampling/2);

NormalizedPupilFunctionArrayY(i,j) = (i-arrayHmiddle)/(PupilSampling/2);

end

end

%% creates dimensional array for Pupil Function and plots Pupil Function

PupilDimensionsArray = linspace(-arraysize(1,2)/2,arraysize(1,2)/2,arraysize(1,2));

PupilDimensions = PupilDimensionsArray.*(ApertureSampleDimension);

PupilDimensionsHorizontal = PupilDimensions;

PupilDimensionsVertical = (PupilDimensions);

disp('Pupil Array generated')

%% fourier transform of the pupil function using a 2D cartesian FFT

n = arraysize(1,1);

m = arraysize(1,2);

Xs = ApertureSampleDimension; % mm

% calculation of the amplitude PSF

AmplitudePSF = fftshift(ifft2(ifftshift(PupilFunctionArray))).*n.*m;

% calculation of the intensity PSF

IntensityPSF = AmplitudePSF.*conj(AmplitudePSF);

% calculation of spatial frequencies in both image dimensions

SpatialFrequenciesX = (-((n/2)-0.5):((n/2)-0.5))./n./Xs;

SpatialFrequenciesY = (-((m/2)-0.5):((m/2)-0.5))./m./Xs;

% calculation of the dimensions in the image plane (in micron) converted

% from spatial frequencies

%PSFdimensionX = SpatialFrequenciesX.*FocalLength.*lambda./ApertureDiameter; %image space dimensions in um

%PSFdimensionY = SpatialFrequenciesX.*FocalLength.*lambda./ApertureDiameter; %image space dimensions in um

PSFdimensionX = SpatialFrequenciesX.*FocalLength.*lambda;

PSFdimensionY = SpatialFrequenciesY.*FocalLength.*lambda;

disp('Diff. Limited Intensity PSF calculated')

%% calculates the diffraction limited OTF and MTF of the system

% the spatial frequencies are not correctly calculated for this MTF. The

% cut off frequency should be 1/lambda/Fnum = 0.79 cycles/mm

MTF_SpatialFrequenciesX = PSFdimensionX./FocalLength./lambda;

MTF_SpatialFrequenciesY = PSFdimensionY./FocalLength./lambda;

OTF = ifftshift(fft2(fftshift(IntensityPSF)))./n./m;

disp('Diff. Limited OTF Calculated')

%% calculation of an aberrated PSF

% aberrated transmission function

% waves of aberration

% Wavefront Error (WFE) in the pupil

% Optical Path Difference (OPD) in the pupil

% phase function as a function of aperture

% coded in switch/case in order to facilitate faster processing

switch WavesOfDefocus

case 0

PhaseFunction_020 = 1;

otherwise

W_020 = WavesOfDefocus*lambda;

WFE_020 = W_020.*((NormalizedPupilFunctionArray).^2);

OPD_020 = (2.*pi./lambda).*WFE_020;

PhaseFunction_020 = exp(1i*OPD_020);

end

switch WavesOfTiltX

case 0

PhaseFunction_111x = 1;

otherwise

W_111x = WavesOfTiltX*lambda;

WFE_111x = W_111x.*(NormalizedFieldX).*(NormalizedPupilFunctionArrayX);

OPD_111x = (2.*pi./lambda).*WFE_111x;

PhaseFunction_111x = exp(1i*OPD_111x);

end

switch WavesOfTiltY

case 0

PhaseFunction_111y = 1;

otherwise

W_111y = WavesOfTiltY*lambda;

WFE_111y = W_111y.*(NormalizedFieldY).*(NormalizedPupilFunctionArrayY);

OPD_111y = (2.*pi./lambda).*WFE_111y;

PhaseFunction_111y = exp(1i*OPD_111y);

end

switch WavesOfSpherical

case 0

PhaseFunction_040 = 1;

otherwise

W_040 = WavesOfSpherical*lambda;

WFE_040 = W_040.*((NormalizedPupilFunctionArray).^4);

OPD_040 = (2.*pi./lambda).*WFE_040;

PhaseFunction_040 = exp(1i*OPD_040);

end

switch WavesOfComaX

case 0

PhaseFunction_131x = 1;

otherwise

W_131x = WavesOfComaX*lambda;

WFE_131x = W_131x.*NormalizedFieldX.*NormalizedPupilFunctionArrayX.*((NormalizedPupilFunctionArray).^2);

OPD_131x = (2.*pi./lambda).*WFE_131x;

PhaseFunction_131x = exp(1i*OPD_131x);

end

switch WavesOfComaY

case 0

PhaseFunction_131y = 1;

otherwise

W_131y = WavesOfComaY*lambda;

WFE_131y = W_131y.*NormalizedFieldY.*NormalizedPupilFunctionArrayY.*((NormalizedPupilFunctionArray).^2);

OPD_131y = (2.*pi./lambda).*WFE_131y;

PhaseFunction_131y = exp(1i*OPD_131y);

end

switch WavesOfAstigmatismX

case 0

PhaseFunction_222x = 1;

otherwise

W_222x = WavesOfAstigmatismX*lambda;

WFE_222x = W_222x.*(NormalizedFieldX.^2).*(NormalizedPupilFunctionArrayX.^2);

OPD_222x = (2.*pi./lambda).*WFE_222x;

PhaseFunction_222x = exp(1i*OPD_222x);

end

switch WavesOfAstigmatismY

case 0

PhaseFunction_222y = 1;

otherwise

W_222y = WavesOfAstigmatismY*lambda;

WFE_222y = W_222y.*(NormalizedFieldY.^2).*(NormalizedPupilFunctionArrayY.^2);

OPD_222y = (2.*pi./lambda).*WFE_222y;

PhaseFunction_222y = exp(1i*OPD_222y);

end

switch WavesOfFieldCurvatureX

case 0

PhaseFunction_220x = 1;

otherwise

W_220x = WavesOfFieldCurvatureX*lambda;

WFE_220x = W_220x.*(NormalizedFieldX.^2).*(NormalizedPupilFunctionArray.^2);

OPD_220x = (2.*pi./lambda).*WFE_220x;

PhaseFunction_220x = exp(1i*OPD_220x);

end

switch WavesOfFieldCurvatureY

case 0

PhaseFunction_220y = 1;

otherwise

W_220y = WavesOfFieldCurvatureY*lambda;

WFE_220y = W_220y.*(NormalizedFieldY.^2).*(NormalizedPupilFunctionArray.^2);

OPD_220y = (2.*pi./lambda).*WFE_220y;

PhaseFunction_220y = exp(1i*OPD_220y);

end

switch WavesOfDistortionX

case 0

PhaseFunction_311x = 1;

otherwise

W_311x = WavesOfDistortionX*lambda;

WFE_311x = W_311x.*(NormalizedFieldX.^3).*(NormalizedPupilFunctionArrayX);

OPD_311x = (2.*pi./lambda).*WFE_311x;

PhaseFunction_311x = exp(1i*OPD_311x);

end

switch WavesOfDistortionY

case 0

PhaseFunction_311y = 1;

otherwise

W_311y = WavesOfDistortionY*lambda;

WFE_311y = W_311y.*(NormalizedFieldY.^3).*(NormalizedPupilFunctionArrayY);

OPD_311y = (2.*pi./lambda).*WFE_311y;

PhaseFunction_311y = exp(1i*OPD_311y);

end

disp('All Aberration Phase Functions Calculated')

% aberrated transmission function

%{AberratedPupilFunction = PupilFunctionArray.*PhaseFunction_111x.*PhaseFunction_111y...

.*PhaseFunction_020.*PhaseFunction_040...

.*PhaseFunction_131x.*PhaseFunction_131y...

.*PhaseFunction_222x.*PhaseFunction_222y...

.*PhaseFunction_220x.*PhaseFunction_220y...

.*PhaseFunction_311x.*PhaseFunction_311y;

% calculation of the aberrated amplitude PSF

AberratedAmplitudePSF = fftshift(ifft2(ifftshift(AberratedPupilFunction))).*n.*m;

% calculation of the intensity PSF

AberratedIntensityPSF = AberratedAmplitudePSF.*conj(AberratedAmplitudePSF);

disp('Aberrated Intensity PSF calculated')

%% Calculation of the aberrated OTF

AberratedOTF = ifftshift(fft2(fftshift(AberratedIntensityPSF)))./n./m;

disp('Aberrated OTF Calculated')

%% Calculation of the Pixel Response Function

%PixelPitch = 1; % 1 um

% PupilFunctionArray = zeros(PupilSampling+(2*ZeroPaddingSize),PupilSampling+(2*ZeroPaddingSize));

% PSFdimensionX = SpatialFrequenciesX.*FocalLength.*lambda;

% PSFdimensionY = SpatialFrequenciesY.*FocalLength.*lambda;

%%

timeElapsed = toc;

%% code for various plots

% % plots Pupil Function

% figure, imagesc(PupilDimensionsHorizontal,PupilDimensionsVertical,PupilFunctionArray),colormap(gray)

% title('Diffraction Limited Pupil Function')

% xlabel('Aperture in X-Dimension (mm)')

% ylabel('Aperture in Y-Dimension (mm)')

%

% % plot of the natural log of the Intensity PSF

% figure, imagesc(PSFdimensionX,PSFdimensionY,log(IntensityPSF./(max(max(IntensityPSF))))),colormap(gray)

% title('Diffraction Limited Log Intensity PSF')

% xlabel('PSF in X-Dimension (\mum)')

% ylabel('PSF in Y-Dimension (\mum)')

%

% figure, imagesc(PSFdimensionX,PSFdimensionY,real(IntensityPSF)./(max(max(real(IntensityPSF))))),colormap(gray)

% title('Diffraction Limited Intensity PSF')

% xlabel('PSF in X-Dimension (\mum)')

% ylabel('PSF in Y-Dimension (\mum)')

%

% % cross section across the center of the real part of the amplitude PSF

% figure, plot(PSFdimensionX,(real(AmplitudePSF(arrayHmiddle,:))./max(max(AmplitudePSF))));

% title('Diffraction Limited Amplitude PSF Cross-section')

% xlabel('PSF in X-Dimension (\mum)')

% ylabel('Normalized Amplitude')

%

% % plot of the diffraction limited OTF

% figure, imagesc(MTF_SpatialFrequenciesX,MTF_SpatialFrequenciesX,real(OTF)),colormap(gray)

% title('Diffraction Limited OTF')

% xlabel('X-Spatial Frequencies (1/mm)')

% ylabel('Y-Spatial Frequencies (1/mm)')

%

% % cross section across the center of the real part of the OTF

% % THIS IS THE MTF

% figure, plot(MTF_SpatialFrequenciesX,abs((OTF(arrayHmiddle,:)))./max(max(real(OTF))));

% title('Diffraction Limited MTF')

% xlabel('X-Spatial Frequencies (1/mm)')

% ylabel('Normalized MTF')

% plot of the aberrated Pupil Function

figure,imagesc(real(AberratedPupilFunction)),colormap(gray);

title('Aberrated Pupil Function')

xlabel('Aperture in X-Dimension (mm)')

ylabel('Aperture in Y-Dimension (mm)')

% % plot of the natural log of the aberrated intensity PSF

% figure, imagesc(PSFdimensionX,PSFdimensionY,log(real(AberratedIntensityPSF)./(max(max(real(AberratedIntensityPSF)))))),colormap(gray)

% title('Aberrated Log Intensity PSF')

% xlabel('PSF in X-Dimension (\mum)')

% ylabel('PSF in Y-Dimension (\mum)')

% plot of aberrated intensity PSF

figure, imagesc(PSFdimensionX,PSFdimensionY,real(AberratedIntensityPSF)./(max(max(real(AberratedIntensityPSF))))),colormap(gray)

title('Aberrated Intensity PSF')

xlabel('PSF in X-Dimension (\mum)')

ylabel('PSF in Y-Dimension (\mum)')

% % cross section across the center of the real part of the amplitude PSF

% figure, plot(PSFdimensionX,(real(AberratedAmplitudePSF(arrayHmiddle,:))./max(max(AmplitudePSF))));

% title('Aberrated Amplitude PSF Cross-section')

% xlabel('PSF in X-Dimension (\mum)')

% ylabel('Normalized Amplitude')

%

% % comparison between the aberrated PSF and the diffraction limited PSF.

% % Both normalized to the max of the diffraction limited PSF

figure, plot(PSFdimensionX,(real(IntensityPSF(arrayHmiddle,:))./max(max(IntensityPSF))),PSFdimensionX,(real(AberratedIntensityPSF(arrayHmiddle,:))./max(max(IntensityPSF))));

legend('Diff. Lim.','Aberrated');

title('PSF Comparison')

xlabel('PSF in X-Dimension (\mum)')

ylabel('Normalized Intensity')

axis([-10 10 -0.5 1.2]);

% % comparison between the aberrated PSF and the diffraction limited PSF.

% % Both normalized to the max of the diffraction limited PSF and x-axis

% % normalized to diffraction limited spot size

figure, plot(PSFdimensionX./DiffLimSpotSize,(real(IntensityPSF(arrayHmiddle,:))./max(max(IntensityPSF))),PSFdimensionX./DiffLimSpotSize,(real(AberratedIntensityPSF(arrayHmiddle,:))./max(max(IntensityPSF))));

legend('Diff. Lim.','Aberrated');

title('PSF Comparison (Normalized Spatial Dim)')

xlabel('PSF Width Norm to Diff Lim Spot Size')

ylabel('Normalized Intensity')

axis([-8 8 -0.5 1.2]);}%

%

% % plot of the aberrated OTF

% figure, imagesc(MTF_SpatialFrequenciesX,MTF_SpatialFrequenciesX,real(AberratedOTF)),colormap(gray)

% title('Aberrated OTF')

% xlabel('X-Spatial Frequencies (1/mm)')

% ylabel('Y-Spatial Frequencies (1/mm)')

%

% % cross section across the center of the real part of the aberrated OTF

% % THIS IS THE ABERRATED MTF

% figure, plot(MTF_SpatialFrequenciesX,abs((AberratedOTF(arrayHmiddle,:)))./max(max(real(AberratedOTF))));

% title('Aberrated MTF')

% xlabel('X-Spatial Frequencies (1/mm)')

% ylabel('Normalized MTF')

%

% % comparison plot of the diffraction limited MTF and aberrated MTF

% figure, plot(MTF_SpatialFrequenciesX,(abs((AberratedOTF(arrayHmiddle,:)))./max(max(real(OTF)))),MTF_SpatialFrequenciesX,(abs(real(OTF(arrayHmiddle,:)))./max(max(real(OTF)))));

% legend('Diff. Lim.','Aberrated');

% title('MTF Comparison')

% xlabel('X-Spatial Frequencies (1/mm)')

% ylabel('Normalized MTF')
