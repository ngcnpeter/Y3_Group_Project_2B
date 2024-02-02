%%
% MC_FLIM_FRET_analysis
% WeiYue Chen, Feb 2015
% Copyright (C) 2015, Laser Analytics Group, University of Cambridge
% License - GPL Version 3 or later, http://www.gnu.org/licenses/gpl.html

% Refer to WeiYue Chen et al. Biophysical Journal 2015.
% http://dx.doi.org/10.1016/j.bpj.2015.01.012

% This code is for Multi-Channel FLIM-FRET analysis for the FLIM data 
% obtained with Becker and Hickl SPCM software.

% The function for importing FLIM data from B&H SPCM (.sdt files),
% loadBHfileusingmeasDescBlock, is written by Sean Warren. 

% How to use this code:
% All the parameters that need to be inputted are listed in the cell "INPUT
% PARAMETERS", and commented with the word "INPUT". Please input all the
% parameters needed in this section according to the comments and run the code.

%% INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LASER AND SPECTRA PARAMETERS

f = 40E6; %INPUT: Laser repetition frequency/Hz
    w = 2*pi*f; % Laser repetition rate(angular speed)
t_chan = 25/(1E9*256); %INPUT: Time interval of each time bin
    t_serie = ((0:255).*t_chan)';
beta=5.03; %INPUT: beta value for spectra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLAG SWITCHES FOR PHASOR ANALYSIS

flagAcceptorSingle=1; %INPUT: choose 1 if assuming pure acceptor decay is
% single exponential; choose 0 if assuming pure acceptor decay is
% bi-exponential.
    tauA=1.39234E-9;%INPUT: average acceptor lifetime
flagInputDonorParameter=0; %INPUT: choose 1 if it is necessary to input donor 
% phasor trajectory parameters manually - useful when the active and passive 
% donor phasors/lifetimes are known, while the data to be analyzed does not 
% have enough distribution of bound donor fraction to recover the active and 
% passive donor phasors. Choose 0 if analyzing the active and passive donor
% phasors directly from the inputted data.
    if flagInputDonorParameter==1
       donorPoly(1)=-0.426436505126505;%INPUT
       donorPoly(2)=0.715451119728850;%INPUT
       % Notes: donor phasor trajectory: 
       % S=donorPoly(1).*G+donorPoly(2)
       [GDonor,SDonor,tauDonor]= tauComp(donorPoly(1),donorPoly(2),w); 
       % Recover the phasors and lifetimes of the active and passive donors.
       % e.g. GDonor(1):G coordinate of active donor, GDonor(2): G coordinate of passive donor
    end
flagOffsetComp_dd=1; %INPUT: choose 1 to calculate the intensity offset (baseline) 
% of donor channel FLIM data according to the first "num_offset" data points 
% of a decay curve, and subtract this number from all time bins in 
% the decay curve (the negative values in the results will be changed to 0).
flagOffsetComp_da=1;%INPUT: choose 1 to compensate intensity offset for FRET channel data
flagOffsetComp_aa=1;%INPUT: choose 1 to compensate intensity offset for acceptor channel data
flagOffsetCompRef=1;%INPUT: choose 1 to compensate intensity offset for standard referenece data
    num_offset=9;
flagIntenKdCal=1;%INPUT: choose from 1-3 to calculate dissociation constant
% Kd or donor and acceptor concentrations ([D] and [A]). Choose other numbers to skip
% this analysis.
% When choosing 1, know [A] according to the acceptor channel intensity, calculate
% Kd and [D].
% When choosing 2, donor is homogeneously distributed and [D] is known, calculate
% Kd and [A]
% When choosing 3, Kd is know, calculate [A] and [D].
    slopeAinten=0.16806; % INPUT: For [A] calculation with aa intensity (case1)
    interceptAinten=11.26327; % INPUT: For [A] calculation with aa intensity (case1)
    %Note: [A]/uM=slopeAinten*photon counts+interceptAinten. This
    %equation is for each one image pixel.
    inputConcenD=2; %INPUT: [D]/uM (case2)
    inputKd=27.1857; %INPUT: Kd/uM (case3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STANDARD SAMPLE INFORMATIONS (FOR IRF CALIBRATION)

%Lifetime of reference data (obtained by fitting with B&H or with other methods)
% NOTE: for our case it is easy to find mono-exponential decay standard
% sample for donor channel and FRET channel, but not for acceptor channel.
% Therefore here we assume that the standard sample for donor and FRET
% channel has mono-exponential decay, while the sample for acceptor channel 
% could be multi-exponential decay. It is easy to enable the
% multi-exponential possibility for donor and FRET channel standard sample.

ddRefTau=3.8361E-9; %INPUT: lifetime of standard sample measured in donor channel
daRefTau=3.8768E-9; %INPUT: lifetime of standard sample measured in donor channel

aaRefTau1=2.0112E-9;%INPUT: lifetime component 1 of standard sample measured in acceptor channel
aaRefTau2=1.1081E-9;%INPUT: lifetime component 2 of standard sample measured in acceptor channel
aaRefTau3=1.0621E-9;%INPUT: lifetime component 3 of standard sample measured in acceptor channel
Contri_refaa1=0.1533;%INPUT: intensity fraction of lifetime component 1
Contri_refaa2=0.3523;%INPUT: intensity fraction of lifetime component 2
Contri_refaa3=0.4944;%INPUT: intensity fraction of lifetime component 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTENSITY FILTERS AND BINNING SETUP

% Intensity filters and binning settings for donor channel data
pmin_dd = 60; %INPUT: threshold for minimum photon counts for the FLIM image
pmax_dd = 10000000000; %INPUT: threshold for maximum photon counts for the FLIM image
pminRef_dd = 0; %INPUT: threshold for minimum photon counts for the standard sample data
binIm_dd = 7; %INPUT: binning value for the sample data.Can only be ODD number or 256.
binRef_dd = 256; %INPUT: binning value for the standard sample data.Can only be ODD number or 256.
 
% Intensity filters and binning settings for FRET channel data
pmin_da = 65; %INPUT: threshold for minimum photon counts for the FLIM image
pmax_da = 50000000;%INPUT: threshold for maximum photon counts for the FLIM image
pminRef_da = 0; %INPUT: threshold for minimum photon counts for the standard sample data
binIm_da = binIm_dd; %INPUT: binning value for the sample.Can only be ODD number or 256.
binRef_da = 256; %INPUT: binning value for the standard reference sample.Can only be ODD number or 256.

% Intensity filters and binning settings for acceptor channel data.
pmin_aa = 60; %INPUT: threshold for minimum photon counts for the FLIM image
pmax_aa = 80000;%INPUT: threshold for maximum photon counts for the FLIM image
pminRef_aa = 0; %INPUT: threshold for minimum photon counts for the standard sample data
binIm_aa = 256; %INPUT: binning value for the sample.Can only be ODD number or 256.
binRef_aa = 256; %INPUT: binning value for the standard reference sample.Can only be ODD number or 256.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTOFLUORESCENCE PARAMETER SETUP

% Autofluorescence data need to be analysed separately. 
% (See Szmacinski H. et al. J. Biomed. Opt. 2014)
GAuto_dd=0.5957;%INPUT: Autofluorescence phasor G coordinate in donor channel
SAuto_dd=0.4274;%INPUT: Autofluorescence phasor S coordinate in donor channel
intenAuto_dd=3;%INPUT: average photon counts of autofluorescence in donor channel
% for each pixel (measured with the same parameter as the fluorescence sample)
GAuto_da=0.6425;%INPUT: Autofluorescence phasor G coordinate in FRET channel
SAuto_da=0.3531;%INPUT: Autofluorescence phasor S coordinate in FRET channel
intenAuto_da=9;%INPUT: average photon counts of autofluorescence in FRET channel
% in each pixel.

%NOTE: here there is not autofluorescence compensation for acceptor channel
%measurements, because in our case the autofluorescence is negligible in
%longer wavelengths. But this function can be easily added if needed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLAG AND PARAMETERS FOR PLOTTING IMAGES

flagDA=1; %INPUT: Choose 0 if only want to analyze donor channel phasor. Choose 1
% to analyse FRET channel phasors and do the complete MC-FLIM-FRET analysis.
flagPlotPhasor=1; %INPUT: Choose 1 to plot the phasor plots
flagPlotDonorFrac=1;%INPUT: Choose 1 to plot the donor bound fraction images
flagPlotAcceptorFrac=1;%INPUT: Choose 1 to plot the acceptor bound fraction images
flagPlotKd=0;%INPUT: Choose 1 to plot the dissociation constant distribution images
flagPlotConcenD=1;%INPUT: Choose 1 to plot the donor concentration variation images
flagCompetitor=1; %INPUT: Choose 1 to calculate and plot competitor concentration images 
    if flagDA~=1
        flagPlotAcceptorFrac=0;
        flagPlotKd=0;
        flagPlotConcenD=0;
        flagCompetitor=0;
        flagIntenKdCal=4;
    end
phasor_xlim=[0.48,0.9];%INPUT: x axis range of the phasor plot
phasor_ylim=[0.3,0.67];%INPUT: y axis range of the phasor plot
phasor_xtick=[0.5:0.1:0.9];%INPUT: x axis ticks of phasor plot
phasor_ytick=[0.3:0.1:0.6];%INPUT: y axis ticks of phasor plot

load('jet_new')
CLIM_d=[0.4,0.76]; %INPUT: Color bar range for donor bound fraction image
CLIM_a=[0,0.1]; %INPUT: Color bar range for acceptor bound fraction image
CLIM_Kd=[0,100]; %INPUT: Color bar range for Kd image
CLIM_Compet=[0,150]; %INPUT: Color bar range for competitor concentration image
CLIM_concenD=[0,10]; %INPUT: Color bar [D] according to [A]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FILE LOCATION

%INPUT: Donor channel sample data file location
ddSdt='F:\sample\sample.sdt';
%INPUT: FRET channel sample data file location
daSdt='F:\sample\sample.sdt';
%INPUT: Acceptor channel sample data file location
aaSdt='F:\sample\sample.sdt';

%INPUT: Standard sample donor channel file location
ddRefSdt='F:\sample\sample.sdt';
%INPUT: Standard sample FRET channel file location
daRefSdt='F:\sample\sample.sdt';
%INPUT: Standard sample acceptor channel file location
aaRefSdt='F:\sample\sample.sdt';

%% CALCULATIONS FOR IRF CALIBRATION PARAMETERS
mRef_dd = (1 + (w*ddRefTau).^2).^(-0.5); % Modulation value that the donor channel standard sample should have.
phiRef_dd = atan(w*ddRefTau); % Phase value that the donor channel standard sample should have.

mRef_da = (1 + (w*daRefTau).^2).^(-0.5); % Modulation value that the FRET channel standard sample should have.
phiRef_da = atan(w*daRefTau); % Phase value that the FRET channel standard sample should have.

GRef_aa1=1/(1 + (w*aaRefTau1).^2);
SRef_aa1=(w*aaRefTau1)/(1 + (w*aaRefTau1).^2);
GRef_aa2=1/(1 + (w*aaRefTau2).^2);
SRef_aa2=(w*aaRefTau2)/(1 + (w*aaRefTau2).^2);
GRef_aa3=1/(1 + (w*aaRefTau3).^2);
SRef_aa3=(w*aaRefTau3)/(1 + (w*aaRefTau3).^2);
GRef_aa_cali=(GRef_aa1*Contri_refaa1*aaRefTau1+GRef_aa2*Contri_refaa2*aaRefTau2+GRef_aa3*Contri_refaa3*aaRefTau3)/...
    (Contri_refaa1*aaRefTau1+Contri_refaa2*aaRefTau2+Contri_refaa3*aaRefTau3);
SRef_aa_cali=(SRef_aa1*Contri_refaa1*aaRefTau1+SRef_aa2*Contri_refaa2*aaRefTau2+SRef_aa3*Contri_refaa3*aaRefTau3)/...
    (Contri_refaa1*aaRefTau1+Contri_refaa2*aaRefTau2+Contri_refaa3*aaRefTau3);

mRef_aa = (GRef_aa_cali^2+SRef_aa_cali^2).^0.5; % Modulation value that the acceptor channel standard sample should have.
phiRef_aa = atan(SRef_aa_cali/GRef_aa_cali); % Phase value that the acceptor channel standard sample should have.

%% IRF CALIBRATION FOR DONOR EXCITATION DONOR EMISSION

% Import data. 
% Note that the function "loadBHfileusingmeasDescBlock" for
% importing Becker & Hickl FLIM data is written by Dr. Sean Warren.
ddref = loadBHfileusingmeasDescBlock(ddRefSdt);
ddref=permute(ddref,[3,2,1]); 
% dataMatrix(y dimension, x dimension, time dimension)
ddref=double(ddref);

% Intensity Matrix of the sample
inten_Ref_dd = sum(ddref,3);

% Create intensity mask
mask_Ref_dd = inten_Ref_dd>pminRef_dd;
mask_Ref_dd = repmat(mask_Ref_dd,[1,1,256]);
ref_mask_dd = ddref.*mask_Ref_dd;

% Binning the image
ref2_dd=convn(ref_mask_dd, ones(binRef_dd),'valid');

% Filter out intensity 0 points (including points filtered out by the 
% intensity mask) before doing the Fourier transform
filter_Ref_dd = find(sum(ref2_dd,3));
ref3_dd = reshape(ref2_dd,[],256);
ref3_dd = ref3_dd';
[sizeFilterRef_dd,test] = size(filter_Ref_dd);
if flagOffsetCompRef==1
    %Subtracting the offset of the decay curve
    ref3_dd=ref3_dd-repmat(mean(ref3_dd([1:num_offset],:)),256,1);
    %Assign negative values 0 value instead.
    ref3_dd=ref3_dd-ref3_dd.*(ref3_dd<0);
end
intenRef_dd = sum(ref3_dd);

% Calculate the phasor representation of the raw data of the standard sample
gRef_dd = (sum(ref3_dd(:,filter_Ref_dd).*repmat(cos(w.*t_serie),[1,sizeFilterRef_dd]))./intenRef_dd(filter_Ref_dd));
sRef_dd = (sum(ref3_dd(:,filter_Ref_dd).*repmat(sin(w.*t_serie),[1,sizeFilterRef_dd]))./intenRef_dd(filter_Ref_dd));

% Calculate the averaged phasor parameters of the standard sample raw data
GRef_dd = mean(gRef_dd);
SRef_dd = mean(sRef_dd);

% Calculate the modulation and phase correction from the standard sample
mCor_dd = mRef_dd./(sqrt(GRef_dd.^2+SRef_dd.^2));
phiCor_dd = phiRef_dd - atan2(SRef_dd,GRef_dd);

%%
% IRF CALIBRATION FOR DONOR EXCITATION ACCEPTOR EMISSION
if flagDA==1
    % Import data. 
    % Note that the function "loadBHfileusingmeasDescBlock" for
    % importing Becker & Hickl FLIM data is written by Dr. Sean Warren.
    daref = loadBHfileusingmeasDescBlock(daRefSdt);
    daref = permute(daref,[3,2,1]);
    % dataMatrix(y dimension, x dimension, time dimension)
    daref = double(daref);

    % Intensity Matrix of the sample
    inten_Ref_da = sum(daref,3);

    % Create intensity mask
    mask_Ref_da = inten_Ref_da>pminRef_da;
    mask_Ref_da = repmat(mask_Ref_da,[1,1,256]);
    ref_mask_da = daref.*mask_Ref_da;

    % Binning the image
     ref2_da=convn(ref_mask_da, ones(binRef_da),'valid');

    % Filter out intensity 0 points (including points filtered out by the 
    % intensity mask) before doing the Fourier transform
    filter_Ref_da = find(sum(ref2_da,3));
    ref3_da = reshape(ref2_da,[],256);
    ref3_da = ref3_da';
    [sizeFilterRef_da,test] = size(filter_Ref_da);
    if flagOffsetCompRef==1
        %Subtracting the offset of the decay curve
        ref3_da=ref3_da-repmat(mean(ref3_da([1:num_offset],:)),256,1);
        %Assign negative values 0 value instead.
        ref3_da=ref3_da-ref3_da.*(ref3_da<0);
    end
    intenRef_da = sum(ref3_da);

    % Calculate the phasor representation of the raw data of the standard sample
    gRef_da = (sum(ref3_da(:,filter_Ref_da).*repmat(cos(w.*t_serie),[1,sizeFilterRef_da]))./intenRef_da(filter_Ref_da));
    sRef_da = (sum(ref3_da(:,filter_Ref_da).*repmat(sin(w.*t_serie),[1,sizeFilterRef_da]))./intenRef_da(filter_Ref_da));

    % Calculate the averaged phasor parameters of the standard sample raw data
    GRef_da = mean(gRef_da);
    SRef_da = mean(sRef_da);
    % Calculate the modulation and phase correction from the standard sample
    mCor_da = mRef_da./(sqrt(GRef_da.^2+SRef_da.^2));
    phiCor_da = phiRef_da - atan2(SRef_da,GRef_da);
end

%%
% IRF CALIBRATION FOR ACCEPTOR ONLY SAMPLE
if flagDA==1
    % Import data. 
    % Note that the function "loadBHfileusingmeasDescBlock" for
    % importing Becker & Hickl FLIM data is written by Dr. Sean Warren.
    aaref = loadBHfileusingmeasDescBlock(aaRefSdt);
    aaref = permute(aaref,[3,2,1]);
    aaref = double(aaref);

    % Intensity Matrix of the sample
    inten_Ref_aa = sum(aaref,3);

    % Create intensity mask
    mask_Ref_aa = inten_Ref_aa>pminRef_aa;
    mask_Ref_aa = repmat(mask_Ref_aa,[1,1,256]);
    ref_mask_aa = aaref.*mask_Ref_aa;

    % Binning the image
    ref2_aa=convn(ref_mask_aa, ones(binRef_aa),'valid');


    % Filter out intensity 0 points (including points filtered out by the 
    % intensity mask) before doing the Fourier transform
    filter_Ref_aa = find(sum(ref2_aa,3));
    ref3_aa = reshape(ref2_aa,[],256);
    ref3_aa = ref3_aa';
    [sizeFilterRef_aa,test] = size(filter_Ref_aa);
    if flagOffsetCompRef==1
        %Subtracting the offset of the decay curve
        ref3_aa=ref3_aa-repmat(mean(ref3_aa([1:num_offset],:)),256,1);
        %Assign negative values 0 value instead
        ref3_aa=ref3_aa-ref3_aa.*(ref3_aa<0);
    end
    intenRef_aa = sum(ref3_aa);

    % Calculate the phasor representation of the raw data of the standard sample
    gRef_aa = (sum(ref3_aa(:,filter_Ref_aa).*repmat(cos(w.*t_serie),[1,sizeFilterRef_aa]))./intenRef_aa(filter_Ref_aa));
    sRef_aa = (sum(ref3_aa(:,filter_Ref_aa).*repmat(sin(w.*t_serie),[1,sizeFilterRef_aa]))./intenRef_aa(filter_Ref_aa));

    % Calculate the averaged phasor parameters of the standard sample raw data
    GRef_aa = mean(gRef_aa);
    SRef_aa = mean(sRef_aa);
    % Calculate the modulation and phase correction from the standard sample
    mCor_aa = mRef_aa./(sqrt(GRef_aa.^2+SRef_aa.^2));
    phiCor_aa = phiRef_aa - atan2(SRef_aa,GRef_aa);
end

%% DONOR (and FRET) CHANNEL PHASOR CALCULATION

% Import donor channel data
dataMatrix_dd= loadBHfileusingmeasDescBlock(ddSdt);
dataMatrix_dd = permute(dataMatrix_dd,[3,2,1]);
dataMatrix_dd = double(dataMatrix_dd);
intenMatrix_dd = sum(dataMatrix_dd,3);

% Applying the intensity mask
mask_P_dd = (intenMatrix_dd>pmin_dd)&(intenMatrix_dd<pmax_dd);
mask_P_dd = repmat(mask_P_dd,[1,1,256]);
dataMatrix2_dd = dataMatrix_dd.*mask_P_dd;

% Binning the image
dataMatrix2_dd = convn(dataMatrix2_dd, ones(binIm_dd),'valid');
dataMatrix3_dd = reshape(dataMatrix2_dd,[],256);
dataMatrix3_dd = dataMatrix3_dd';

%Subtract the intensity offset from the decay curve:
if flagOffsetComp_dd==1
    dataMatrix3_dd=dataMatrix3_dd-repmat(mean(dataMatrix3_dd([1:num_offset],:)),256,1);
    dataMatrix3_dd=dataMatrix3_dd-dataMatrix3_dd.*(dataMatrix3_dd<0);
end
intenData_dd = sum(dataMatrix3_dd);

if flagDA==1 % FRET channel phasor calculation
    % Import FRET channel data
    dataMatrix_da= loadBHfileusingmeasDescBlock(daSdt);
    dataMatrix_da = permute(dataMatrix_da,[3,2,1]);
    dataMatrix_da = double(dataMatrix_da);
    intenMatrix_da = sum(dataMatrix_da,3);

    % Applying the intensity mask
    mask_P_da = (intenMatrix_da>pmin_da)&(intenMatrix_da<pmax_da);
    mask_P_da = repmat(mask_P_da,[1,1,256]);
    dataMatrix2_da = dataMatrix_da.*mask_P_da;

    % Binning the image
    dataMatrix2_da = convn(dataMatrix2_da, ones(binIm_da),'valid');
    dataMatrix3_da = reshape(dataMatrix2_da,[],256);
    dataMatrix3_da = dataMatrix3_da';
    
    % Subtract the intensity offset from the decay curve:
    if flagOffsetComp_da==1
        dataMatrix3_da=dataMatrix3_da-repmat(mean(dataMatrix3_da([1:num_offset],:)),256,1);
        dataMatrix3_da=dataMatrix3_da-dataMatrix3_da.*(dataMatrix3_da<0);
    end
    intenData_da = sum(dataMatrix3_da);
    % Filter out intensity 0 points before doing the Fourier Transform
    filter = (find(intenData_da.*intenData_dd))';
    [sizeFilter,test] = size(filter);

    % Calculate the phasor representation of sample
    g_da = (sum(dataMatrix3_da(:,filter).*repmat(cos(w.*t_serie),[1,sizeFilter]))./intenData_da(filter));
    s_da = (sum(dataMatrix3_da(:,filter).*repmat(sin(w.*t_serie),[1,sizeFilter]))./intenData_da(filter));

    % IRF correction
    G_da = (g_da.*cos(phiCor_da) - s_da.*sin(phiCor_da)).*mCor_da;
    S_da = (g_da.*sin(phiCor_da) + s_da.*cos(phiCor_da)).*mCor_da;
    
    % Calculate the number of pixels not filtered out in each binned pixel
    % (For following auto-fluorescence compensation)
    dataMatrix2_mask_da = convn(mask_P_da(:,:,1), ones(binIm_da),'valid');
    dataMatrix3_mask_da = reshape(dataMatrix2_mask_da,[],1);
    dataMatrix3_mask_da=dataMatrix3_mask_da';
    num_pixel_auto_da=dataMatrix3_mask_da(filter);
    
    % Autofluorescence compensation
    G_da=(G_da.*intenData_da(filter')-intenAuto_da.*num_pixel_auto_da.*GAuto_da)./(intenData_da(filter')-num_pixel_auto_da.*intenAuto_da);
    S_da=(S_da.*intenData_da(filter')-intenAuto_da.*num_pixel_auto_da.*SAuto_da)./(intenData_da(filter')-num_pixel_auto_da.*intenAuto_da);
else
    filter = (find(intenData_dd))';
    [sizeFilter,test] = size(filter);
end

% Calculate the phasor representation of sample
g_dd = (sum(dataMatrix3_dd(:,filter).*repmat(cos(w.*t_serie),[1,sizeFilter]))./intenData_dd(filter));
s_dd = (sum(dataMatrix3_dd(:,filter).*repmat(sin(w.*t_serie),[1,sizeFilter]))./intenData_dd(filter));

% IRF correction
G_dd = (g_dd.*cos(phiCor_dd) - s_dd.*sin(phiCor_dd)).*mCor_dd;
S_dd = (g_dd.*sin(phiCor_dd) + s_dd.*cos(phiCor_dd)).*mCor_dd;

% Calculate the number of pixels not filtered out in each binned pixel
% (For following auto-fluorescence compensation)
dataMatrix2_mask_dd = convn(mask_P_dd(:,:,1), ones(binIm_dd),'valid');
dataMatrix3_mask_dd = reshape(dataMatrix2_mask_dd,[],1);
dataMatrix3_mask_dd=dataMatrix3_mask_dd';
num_pixel_auto_dd=dataMatrix3_mask_dd(filter);

% Autofluorescence compensation
G_dd=(G_dd.*intenData_dd(filter')-intenAuto_dd.*num_pixel_auto_dd.*GAuto_dd)./(intenData_dd(filter')-num_pixel_auto_dd.*intenAuto_dd);
S_dd=(S_dd.*intenData_dd(filter')-intenAuto_dd.*num_pixel_auto_dd.*SAuto_dd)./(intenData_dd(filter')-num_pixel_auto_dd.*intenAuto_dd);

%%
% CALCULATING PHASOR FOR ACCEPTOR-ONLY SAMPLE
if flagDA==1
    %Import B&H data
    dataMatrix_aa= loadBHfileusingmeasDescBlock(aaSdt);
    dataMatrix_aa= permute(dataMatrix_aa,[3,2,1]);
    dataMatrix_aa= double(dataMatrix_aa);

    % Intensity Matrix of the sample
    intenMatrix_aa = sum(dataMatrix_aa,3);

    % Create intensity mask
    mask_P_aa = (intenMatrix_aa>pmin_aa)&(intenMatrix_aa<pmax_aa);
    mask_P_aa = repmat(mask_P_aa,[1,1,256]);
    data_mask_aa = dataMatrix_aa.*mask_P_aa;

    % Binning the image
    dataMatrix2_aa = convn(data_mask_aa, ones(binIm_aa),'valid');
    inten_conc_aa = convn(intenMatrix_aa, ones(binIm_dd),'valid')./(binIm_dd^2);

    % Filter out intensity 0 points before doing the Fourier Transform
    filter_aa = find(sum(dataMatrix2_aa,3));
    dataMatrix3_aa = reshape(dataMatrix2_aa,[],256);
    dataMatrix3_aa = dataMatrix3_aa';
    [sizeFilter_aa,test] = size(filter_aa);

    % Subtract the intensity offset from the decay curve:
    if flagOffsetComp_aa==1
        dataMatrix3_aa=dataMatrix3_aa-repmat(mean(dataMatrix3_aa([1:num_offset],:)),256,1);
        dataMatrix3_aa=dataMatrix3_aa-dataMatrix3_aa.*(dataMatrix3_aa<0);
    end
    intenData_aa = sum(dataMatrix3_aa);

    % Calculate the phasor representation of sample
    g_aa = (sum(dataMatrix3_aa(:,filter_aa).*repmat(cos(w.*t_serie),[1,sizeFilter_aa]))./intenData_aa(filter_aa));
    s_aa = (sum(dataMatrix3_aa(:,filter_aa).*repmat(sin(w.*t_serie),[1,sizeFilter_aa]))./intenData_aa(filter_aa));

    % IRF calibration
    G_acceptor = (g_aa.*cos(phiCor_aa) - s_aa.*sin(phiCor_aa)).*mCor_aa;
    S_acceptor = (g_aa.*sin(phiCor_aa) + s_aa.*cos(phiCor_aa)).*mCor_aa;

    if binIm_aa~=1
        G_aa=mean(G_acceptor);
        S_aa=mean(S_acceptor);
    else
        G_aa=G_acceptor;
        S_aa=S_acceptor;
    end
end
%%
if flagInputDonorParameter~=1
    [donorPoly,donorStv] = polyfit(G_dd,S_dd,1);
    [GDonor,SDonor,tauDonor]= tauComp(donorPoly(1),donorPoly(2),w);
end
Din = [GDonor(1);SDonor(1)];% Phasor for active donor
Dnon = [GDonor(2);SDonor(2)];% Phasor for passive donor
tauIn_dd=min(tauDonor);% Lifetime of active donor
tauNon_dd=max(tauDonor);% Lifetime of passive donor
E=1-tauIn_dd/tauNon_dd;% FRET efficiency

%%
% CALCULATING THE ACTIVE ACCEPTOR PHASOR
if flagDA==1
    if flagAcceptorSingle~=1 %Bi-exponential decay acceptor
        [aA,bA,accpos1,accpos2,acctau1,acctau2]=decomBiexp(G_aa,S_aa,w);
        % S=aA*G+bA is the line connecting the two acceptor components.
        % accpos1 and accpos2 are the positions of the two lifetime components
        % acctau1 and acctau2 are the corresponding lifetime of the two
        % components
        [GAfret,SAfret]=aFret(Din,accpos1,accpos2,w);
        %GAfret and SAfret are the G and S coordinates of the active acceptor
    else % Single exponential decay acceptor
        [GAfret,SAfret]=aFretSingle(Din,tauA,w);
    end
    Ain = [GAfret;SAfret];% Active acceptor phasor
    Anon = [G_aa;S_aa];% Passive acceptor phasor
end

%% CALCULATING THE ACCEPTOR-CONTRIBUTION-ONLY PHASOR tauAM (MC-FLIM-FRET)
if flagDA==1
    [a,b,G_AM,S_AM] = acceptorFLIMFRET_batch(G_dd,S_dd,G_da,S_da,Ain,Anon);
end

%% DONOR ACTIVE AND INTERACTING FRACTION CALCULATION
% Active donor fraction is the same as bound donor fraction
Gfoot_dd= (G_dd+donorPoly(1).*S_dd-donorPoly(1)*donorPoly(2))./(donorPoly(1)^2+1);
Sfoot_dd= (donorPoly(1).*G_dd+donorPoly(1)^2.*S_dd+donorPoly(2))./(donorPoly(1)^2+1);

fracInten_dIn=(Gfoot_dd-Dnon(1))./(Din(1)-Dnon(1));
% Molecular fraction of bound donor:
fracMol_dIn=fracInten_dIn.*tauNon_dd./(fracInten_dIn.*tauNon_dd+(1-fracInten_dIn).*tauIn_dd);

%% ACCEPTOR ACTIVE AND INTERACTING FRACTION CALCULATION
if flagDA==1
    fracInten_aIn=(G_AM-G_aa)./(GAfret-G_aa);
    fracMol_aIn=fracInten_aIn./(fracInten_aIn.*(1-beta*E)+beta*E);% Molecular fraction of active acceptor
    fracMol_aIn_true=fracMol_aIn./(1-fracMol_aIn);%Molecular fraction of bound acceptor
end

%% ASSIGN VALUE BACK TO ORIGINAL IMAGE
[subx,suby]=ind2sub([257-binIm_dd,257-binIm_dd],filter);
fracMol_dIn_Mat=accumarray([subx,suby],fracMol_dIn,[257-binIm_dd,257-binIm_dd]);
sub_Mat=accumarray([subx,suby],1,[257-binIm_dd,257-binIm_dd]);

G_dd_Mat=accumarray([subx,suby],G_dd,[257-binIm_dd,257-binIm_dd]);
S_dd_Mat=accumarray([subx,suby],S_dd,[257-binIm_dd,257-binIm_dd]);

if flagDA==1
    fracMol_aIn_Mat=accumarray([subx,suby],fracMol_aIn_true,[257-binIm_dd,257-binIm_dd]);
    G_da_Mat=accumarray([subx,suby],G_da,[257-binIm_dd,257-binIm_dd]);
    S_da_Mat=accumarray([subx,suby],S_da,[257-binIm_dd,257-binIm_dd]);
end

%% MAKING PHASOR PLOTS
if flagPlotPhasor==1;
    figure(1)
    bcirc=linspace(0,1,1000);
    acirc=sqrt(0.25-(bcirc-0.5).^2);
    plot(bcirc,acirc,'k','LineWidth',2)
    set(gca,'Fontsize',15);
    xlabel('B');
    ylabel('A')
    hold on

    if flagDA==1
        load('cmap_phasor2')
        test=linspace(min(GDonor)-0.005,max(GDonor)+0.005,100000);
        plot(test,donorPoly(1).*test+donorPoly(2),'g','linewidth',2);
        a_slope=(Ain(2)-Anon(2))/(Ain(1)-Anon(1));
        a_intersection=(Ain(1)*Anon(2)-Anon(1)*Ain(2))/(Ain(1)-Anon(1));
        a_test=linspace(0.5,0.9,100000);
        plot(a_test,a_slope.*a_test+a_intersection,'color','r','linewidth',3);
        dscatter_multicolor(G_dd',S_dd',G_da',S_da',G_AM',S_AM',1,cmap_phasor2)
        scatter(GDonor(1),SDonor(1),110,[0,0.3,0],'filled','s');
        scatter(GDonor(2),SDonor(2),110,[0,0.3,0],'filled','^');
        scatter(G_aa,S_aa,110,[0.5,0,0],'filled','^');
        scatter(Ain(1),Ain(2),110,[0.5,0,0],'filled','s');
    else
        load('cmap_dd')
        test=linspace(min(GDonor)-0.005,max(GDonor)+0.005,100000);
        plot(test,donorPoly(1).*test+donorPoly(2),'g','linewidth',2);
        dscatter_dd(G_dd',S_dd',1,cmap_dd)
        scatter(GDonor(1),SDonor(1),110,[0,0.3,0],'filled','s');
        scatter(GDonor(2),SDonor(2),110,[0,0.3,0],'filled','^');
    end
    set(gca,'xlim',phasor_xlim,'ylim',phasor_ylim,'xtick',phasor_xtick,'ytick',phasor_ytick)
end
%% MAKING BOUND FRACTION FIGURES
if flagDA==1
    binary_image=mask_P_dd(:,:,1).*mask_P_da(:,:,1).*mask_P_aa(:,:,1);
else
    binary_image=mask_P_dd(:,:,1);
end
cutEdge=(binIm_dd-1)/2;
if flagPlotDonorFrac==1;
    inten_dd=sum(dataMatrix_dd,3);
    inten_dd=inten_dd./(max(inten_dd(:)));

    fracMol_dIn_color=fracMol_dIn_Mat;
    colorInterval_d=(CLIM_d(2)-CLIM_d(1))/(64-1);
    fracMol_dIn_color = (fracMol_dIn_color-CLIM_d(1))./colorInterval_d+2;
    fracMol_dIn_color = round(fracMol_dIn_color);

    mask_color_d_up=fracMol_dIn_color>65;
    mask_color_d_low=fracMol_dIn_color<2;
    mask_color_fd_low=fracMol_dIn_Mat<0;
    mask_color_fd_up=fracMol_dIn_Mat>1;

    test_fd_low=mask_color_fd_low+mask_color_d_low;
    mask_color_d_low=mask_color_d_low-(test_fd_low>1);

    test_fd_up=mask_color_fd_up+mask_color_d_up;
    mask_color_d_up=mask_color_d_up-(test_fd_up>1);

    fracMol_dIn_color=fracMol_dIn_color.*(1-mask_color_d_up-mask_color_d_low-mask_color_fd_low-mask_color_fd_up)...
        +mask_color_d_up.*65+mask_color_d_low.*2+mask_color_fd_low.*1+mask_color_fd_up.*66;

    fracMol_dIn_rgb = ind2rgb(fracMol_dIn_color,jet_new);
    fracMol_dIn_rgb=fracMol_dIn_rgb.*(repmat(sub_Mat,[1,1,3]));

    fracMol_dIn_rgb_inten=fracMol_dIn_rgb.*repmat(inten_dd(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]);
    % Bound fraction image of donor
    figure(2)
    imshow(fracMol_dIn_rgb.*repmat(binary_image(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]))
    title('Donor bound fraction')
    
    % Bound fraction image of donor (multiplied with the intensity image)
    figure(3)
    imshow(fracMol_dIn_rgb_inten.*repmat(binary_image(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]))
    title('Donor bound fraction - intensity merged')
end

if flagPlotAcceptorFrac==1; 
    inten_aa=sum(dataMatrix_aa,3);
    inten_aa=inten_aa./(max(inten_aa(:)));
    
    fracMol_aIn_color=fracMol_aIn_Mat;
    colorInterval_a=(CLIM_a(2)-CLIM_a(1))/(64-1);
    fracMol_aIn_color = (fracMol_aIn_color-CLIM_a(1))./colorInterval_a+2;
    fracMol_aIn_color = round(fracMol_aIn_color);

    mask_color_a_up=fracMol_aIn_color>65;
    mask_color_a_low=fracMol_aIn_color<2;
    mask_color_fa_low=fracMol_aIn_Mat<0;
    mask_color_fa_up=fracMol_aIn_Mat>1;

    test_fa_low=mask_color_fa_low+mask_color_a_low;
    mask_color_a_low=mask_color_a_low-(test_fa_low>1);

    test_fa_up=mask_color_fa_up+mask_color_a_up;
    mask_color_a_up=mask_color_a_up-(test_fa_up>1);

    fracMol_aIn_color=fracMol_aIn_color.*(1-mask_color_a_up-mask_color_a_low-mask_color_fa_low-mask_color_fa_up)...
        +mask_color_a_up.*65+mask_color_a_low.*2+mask_color_fa_low.*1+mask_color_fa_up.*66;

    fracMol_aIn_rgb = ind2rgb(fracMol_aIn_color,jet_new);
    fracMol_aIn_rgb=fracMol_aIn_rgb.*(repmat(sub_Mat,[1,1,3]));

    fracMol_aIn_rgb_inten=fracMol_aIn_rgb.*repmat(inten_aa(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]);
    % Bound fraction image of the acceptor
    figure(4)
    imshow(fracMol_aIn_rgb.*repmat(binary_image(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]))
    title('Acceptor bound fraction')
    % Bound fraction image of the acceptor (multiplied with the intensity image)
    figure(5)
    imshow(fracMol_aIn_rgb_inten.*repmat(binary_image(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]))
    title('Acceptor bound fraction - intensity merged')
end
%% CALCULATE AND PLOT KD
if flagIntenKdCal==1
    concenA=inten_conc_aa.*slopeAinten+interceptAinten;
    concenD=concenA.*fracMol_aIn_Mat./fracMol_dIn_Mat;
    Kd=(1-fracMol_aIn_Mat).*(1-fracMol_dIn_Mat)./fracMol_dIn_Mat.*concenA;
elseif flagIntenKdCal==2
    concenD=ones(size(fracMol_dIn_Mat)).*inputConcenD;
    concenA=concenD.*fracMol_dIn_Mat./fracMol_aIn_Mat;
    Kd=(1-fracMol_aIn_Mat).*(1-fracMol_dIn_Mat)./fracMol_aIn_Mat.*concenD;
elseif flagIntenKdCal==3
    Kd=ones(size(fracMol_dIn_Mat)).*inputKd;
    concenD=Kd.*fracMol_aIn_Mat./(1-fracMol_aIn_Mat)./(1-fracMol_dIn_Mat);
    concenA=concenD./fracMol_aIn_Mat.*fracMol_dIn_Mat;
end

if flagPlotKd==1;
    Kd_color=Kd;
    colorInterval_Kd=(CLIM_Kd(2)-CLIM_Kd(1))/(64-1);
    Kd_color = (Kd_color-CLIM_Kd(1))./colorInterval_Kd+2;
    Kd_color = round(Kd_color);

    mask_color_Kd_up=Kd_color>65;
    mask_color_Kd_low=Kd_color<2;
    mask_color_f=1-(fracMol_dIn_Mat>0).*(fracMol_dIn_Mat<1).*(fracMol_aIn_Mat>0).*(fracMol_aIn_Mat<1);

    test_fd_low_Kd=mask_color_f+mask_color_Kd_low;
    mask_color_Kd_low=mask_color_Kd_low-(test_fd_low_Kd>1);

    test_fd_up_Kd=mask_color_f+mask_color_Kd_up;
    mask_color_Kd_up=mask_color_Kd_up-(test_fd_up_Kd>1);

    Kd_color=Kd_color.*(1-mask_color_Kd_up-mask_color_Kd_low-mask_color_f)...
        +mask_color_Kd_up.*65+mask_color_Kd_low.*2+mask_color_f.*66;

    fracMol_Kd_rgb = ind2rgb(Kd_color,jet_new);
    fracMol_Kd_rgb=fracMol_Kd_rgb.*(repmat(sub_Mat,[1,1,3]));
    
    figure(6)
    imshow(fracMol_Kd_rgb.*repmat(binary_image(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]))
    title('Dissociation constant Kd/\muM')
end

%% CALCULATE COMPETETOR CONCENTRATION
if flagCompetitor==1
    KdTrue=ones(size(fracMol_dIn_Mat)).*inputKd;
    % Competitor concentration:
    concenCompet=(concenA.*(1-fracMol_aIn_Mat)-KdTrue.*fracMol_dIn_Mat./(1-fracMol_dIn_Mat))./fracMol_dIn_Mat;
      
    Compet_color=concenCompet;
    colorInterval_Kd=(CLIM_Compet(2)-CLIM_Compet(1))/(64-1);
    Compet_color = (Compet_color-CLIM_Compet(1))./colorInterval_Kd+2;
    Compet_color = round(Compet_color);

    mask_color_Compet_up=Compet_color>65;
    mask_color_Compet_low=Compet_color<2;
    mask_color_f=1-(fracMol_dIn_Mat>0).*(fracMol_dIn_Mat<1).*(fracMol_aIn_Mat>0).*(fracMol_aIn_Mat<1);

    test_fd_low_Compet=mask_color_f+mask_color_Compet_low;
    mask_color_Compet_low=mask_color_Compet_low-(test_fd_low_Compet>1);

    test_fd_up_Compet=mask_color_f+mask_color_Compet_up;
    mask_color_Compet_up=mask_color_Compet_up-(test_fd_up_Compet>1);

    Compet_color=Compet_color.*(1-mask_color_Compet_up-mask_color_Compet_low-mask_color_f)...
        +mask_color_Compet_up.*65+mask_color_Compet_low.*2+mask_color_f.*66;

    Compet_rgb = ind2rgb(Compet_color,jet_new);
    Compet_rgb=Compet_rgb.*(repmat(sub_Mat,[1,1,3]));
    
    figure(7)
    imshow(Compet_rgb.*repmat(binary_image(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]))
    title('Competitor concentration/\muM')
end
%%
if flagPlotConcenD==1
    concenD_color=concenD;
    colorInterval_concenD=(CLIM_concenD(2)-CLIM_concenD(1))/(64-1);
    concenD_color = (concenD_color-CLIM_concenD(1))./colorInterval_concenD+2;
    concenD_color = round(concenD_color);

    mask_color_concenD_up=concenD_color>65;
    mask_color_concenD_low=concenD_color<2;
    mask_color_f=1-(fracMol_dIn_Mat>0).*(fracMol_dIn_Mat<1).*(fracMol_aIn_Mat>0).*(fracMol_aIn_Mat<1);

    test_fd_low_concenD=mask_color_f+mask_color_concenD_low;
    mask_color_concenD_low=mask_color_concenD_low-(test_fd_low_concenD>1);

    test_fd_up_concenD=mask_color_f+mask_color_concenD_up;
    mask_color_concenD_up=mask_color_concenD_up-(test_fd_up_concenD>1);

    concenD_color=concenD_color.*(1-mask_color_concenD_up-mask_color_concenD_low-mask_color_f)...
        +mask_color_concenD_up.*65+mask_color_concenD_low.*2+mask_color_f.*66;

    concenD_rgb = ind2rgb(concenD_color,jet_new);
    concenD_rgb=concenD_rgb.*(repmat(sub_Mat,[1,1,3]));

    figure(8)
    imshow(concenD_rgb.*repmat(binary_image(cutEdge+1:256-cutEdge,cutEdge+1:256-cutEdge),[1,1,3]))
    title('Donor concentration/\muM')
end