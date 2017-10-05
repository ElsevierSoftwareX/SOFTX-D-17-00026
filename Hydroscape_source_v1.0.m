%%%HYDROSCAPE 1.0 Source Code
%%%AUTHOURS: Sean P. Funk, M.Sc.; Danny Hnatyshin, M.Sc.; Daniel S. Alessi,
%%%Ph.D.
%%%Department of Earth & Atmospheric Sciences, University of Alberta,
%%%Edmonton, Alberta, Canada
%%%Last updated: July 29, 2017

%%%This MATLAB code simulates the concentration distribution (plume) using
%%%the Karanovic et al. (2007) solution, allowing the user to input an
%%%arbitrary source region, source function and implement simple geology;

%%%HYDROSCAPE, THE SOURCE CODE, THE USER MANUAL AND ALL OTHER
%%%SUPPLEMENTARY MATERIALS/MEDIA (SUCH AS, BUT NOT LIMITED TO, YOUTUBE(TM)
%%%VIDEOS, POSTERS, POWERPOINTS PRESENTATIONS, ETC.) ARE MADE AVAILABLE
%%%UNDER AN AS-IS BASIS, WITHOUT GUARENTEE OR WARRENTY OF ANY KIND, EXPRESS
%%%OR IMPLIDE. NEITHER THE AUTHORS (SEAN P. FUNK, DANNY HNATYSHIN, DANIEL
%%%S. ALESSI) OF HYDROSCAPE, NOR THOSE WHO AIDED IN ITS DEVELOPMENT OR
%%%REVIEW, ACCEPT ANY LIABILITY RESULTING FROM THE USE OF HYDROSCAPE (OR
%%%ITS SOURCE CODE) OR ITS DOCUMENTATION (E.G. USER MANUAL OR OTHER MEDIA).
%%%IMPLEMTATION OF HYDROSCAPE (OR OTHER MATERIALS) AND INTERPRETATION OF
%%%THE PREDICTIONS OF THE MODEL/SIMULATION ARE THE SOLE RESPONSIBILITY OF
%%%THE USER
%%%THIS SOURCE CODE IS PROVIDED FREE OF CHARGE, AND IS NOT SUBJECT TO ANY
%%%TECH SUPPORT BY ITS AUTHORS. USERS ARE HIGHLY ENCOURAGED TO CHANGE ONLY
%%%THE PARAMETERS THAT ARE NEEDED TO BE CHANGED, AND NOTHING MORE. 

%Domain Inputs
tic
XPos = 100; %X-Length (m)
YPos = 25; %Y-Length (m)
ZPos = 10; %Z-Height (m)
t = 4; %Time (yr)
t = t*3.15569e7; %Converts the time in years to seconds

DisX = 20; %Discretization X
DisY = 20; %Discretization Y
DisZ = 5; %Discretization Z
DisTime = 10; %Discretization t
numofelements = (DisX+1)*(DisY+1)*(DisZ+1)*(DisTime+1); %Calculates the number of points needed to be calculated %DO NOT EDIT

BackgroundNoiseLevel = 0.01; %Background Noise Level (mg/L)

%Nodes (only available option is the default 64)
%Highly recommended that you do not change this
f = [-0.024350293	0.024350293	-0.072993122	0.072993122	-0.121462819	0.121462819	-0.16964442	0.16964442	-0.217423644	0.217423644	-0.264687162	0.264687162	-0.311322872	0.311322872	-0.357220158	0.357220158	-0.402270158	0.402270158	-0.446366017	0.446366017	-0.489403146	0.489403146	-0.531279464	0.531279464	-0.571895646	0.571895646	-0.611155355	0.611155355	-0.648965471	0.648965471	-0.685236313	0.685236313	-0.71988185	0.71988185	-0.752819907	0.752819907	-0.783972359	0.783972359	-0.813265315	0.813265315	-0.840629296	0.840629296	-0.865999398	0.865999398	-0.889315446	0.889315446	-0.910522137	0.910522137	-0.929569172	0.929569172	-0.946411375	0.946411375	-0.9610088	0.9610088	-0.973326828	0.973326828	-0.983336254	0.983336254	-0.991013371	0.991013371	-0.996340117	0.996340117	-0.999305042	0.999305042];
w = [0.048690957	0.048690957	0.048575467	0.048575467	0.048344762	0.048344762	0.047999389	0.047999389	0.047540166	0.047540166	0.046968183	0.046968183	0.046284797	0.046284797	0.045491628	0.045491628	0.044590558	0.044590558	0.043583725	0.043583725	0.042473515	0.042473515	0.041262563	0.041262563	0.039953741	0.039953741	0.038550153	0.038550153	0.037055129	0.037055129	0.035472213	0.035472213	0.033805162	0.033805162	0.032057928	0.032057928	0.030234657	0.030234657	0.028339673	0.028339673	0.02637747	0.02637747	0.024352703	0.024352703	0.022270174	0.022270174	0.020134823	0.020134823	0.017951716	0.017951716	0.01572603	0.01572603	0.013463048	0.013463048	0.011168139	0.011168139	0.00884676	0.00884676	0.006504458	0.006504458	0.004147033	0.004147033	0.001783281	0.001783281];

%Source Region Inputs
%Input data below in the following order Source Position Y (m); Source
%Position Z (m); Source Width (m); Source Thickness (m); Concentration
%(mg/L); For multiple sources, keep match sources in the same order for
%each parameter (e.g. Source1 is first position, Source2 second, etc. 
%Each column in data is its own patch source
data                        = [0 0; -1 -3; 10 10; 2 2; 1000 800]; %Data matrix for the source region geometry
   
YPosSource                  = (data(1,:)); %Source Position Y (m) %DO NOT EDIT
ZPosSource                  = (data(2,:)); %Source Position Z (m) %DO NOT EDIT
Width                       = (data(3,:)); %Source Width (m) %DO NOT EDIT
Thickness                   = (data(4,:)); %Source Thickness (m) %DO NOT EDIT
Concentration               = (data(5,:)); %Source Concentration (mg/L) %DO NOT EDIT

Conc3                       = zeros(numel(Concentration),numofelements); %DO NOT EDIT

%Input the geology data here. First input how many layers you want in your
%domain. Then, in the Geology matrix, input the top and bottom boundary
%data. Each column in the matrix should correspond with the top and bottom
%of a particular layer
%Make sure that tops and bottoms coincide!!
NumberOfLayers = 2; %Input how many geological layers you want in your domain
Geology = [0 3; 3 10]; %The first row is the boundary tops, the second row is the boundary bottoms

%Input the hydraulic and transport parameters for the simulation here. The
%first input in the row vector corresponds with the first layer in the
%simple geology input, the second value corresponds with the second layer,
%etc. 
RetardationInput = [1 1]; %Retardation Factors (-)
SoluteDecayInput = [0 0]; %Solute Decay Constant (1/yr)
VelocityInput = [4E-7 4E-6]; %Average Linear Groundwater Velocity (m/s)
DiffusionInput = [0 0]; %Effective Diffusion Coefficent (m^2/s)
AlphaLInput = [2 2]; %Longitudinal Dispersivity (m)
AlphaTHInput = [0.3 0.3]; %Transverse Horizontal Dispersivity (m)
AlphaTVInput = [0.01 0.01]; %Transverse Vertical Dispersivity (m)

for j = 1:length(Concentration)
    LayerTops               = zeros(1,NumberOfLayers); %DO NOT EDIT
    LayerBottoms            = zeros(1,NumberOfLayers); %DO NOT EDIT
       for k = 1:NumberOfLayers
        LayerBoundaryTop    = Geology(k); %DO NOT EDIT
        LayerBoundaryBottom = Geology(k+NumberOfLayers); %DO NOT EDIT
        
        LayerTops(k)        = LayerBoundaryTop; %DO NOT EDIT
        LayerBottoms(k)     = LayerBoundaryBottom; %DO NOT EDIT
       end

    SourceBottom            = -ZPosSource(j) + Thickness(j)/2; %DO NOT EDIT
    SourceTop               = -ZPosSource(j) - Thickness(j)/2; %DO NOT EDIT

    TopDif                  = SourceTop - LayerTops; %DO NOT EDIT
    TopDif(TopDif<0)        = NaN; %DO NOT EDIT
    FirstLayer              = find(TopDif == min(TopDif)); %DO NOT EDIT

    if isempty(FirstLayer);
        FirstLayer          = 1; %DO NOT EDIT
    end

    BottomDif               = LayerBottoms - SourceBottom; %DO NOT EDIT
    BottomDif(BottomDif < 0) = NaN; %DO NOT EDIT
    LastLayer               = find(BottomDif == min(BottomDif)); %DO NOT EDIT

    if isempty(LastLayer);
        LastLayer           = 1; %DO NOT EDIT
    end

    ValidLayerIndex         = FirstLayer:LastLayer;    %DO NOT EDIT
    Conc2                   = zeros(numel(k), numofelements); %DO NOT EDIT

    for k = ValidLayerIndex                                           % For each layer
        
        Retardation         = RetardationInput(k); %Retardation Factor DO NOT EDIT
        SoluteDecay         = SoluteDecayInput(k); %Solute Decay Constant DO NOT EDIT         
        Velocity            = VelocityInput(k); %Average Linear Groundwater Velocity DO NOT EDIT         
        DiffusionCoefficient = DiffusionInput(k); %Effective Diffusion Coefficient DO NOT EDIT       
        AlphaL              = AlphaLInput(k); %Longitudinal Dispersivity  DO NOT EDIT      
        AlphaTH             = AlphaTHInput(k); %Transverse Horizontal Dispersivity DO NOT EDIT       
        AlphaTV             = AlphaTVInput(k); %Transverse Vertical Dispersivity DO NOT EDIT
        
        Velocity2           = Velocity/Retardation; %Calculates the new average linear groundwater velocity factoring in R (m/s)
        DL                  = DiffusionCoefficient/Retardation + AlphaL*Velocity2; %Calculates the Longitudinal Dispersion Coefficient (m^2/s)
        DTH                 = DiffusionCoefficient/Retardation + AlphaTH*Velocity2; %Calculates the Transverse Horizontal Dispersion Coefficient (m^2/s)
        DTV                 = DiffusionCoefficient/Retardation + AlphaTV*Velocity2; %Calculates the Transverse Vertical Dispersion Coefficient (m^2/s)
        lambda              = SoluteDecay/3.15569E7; %Converts the Solute Decay Constant to (1/sec)
        
        ks = 0; %Source Decay Constant (1/yr); must be set to zero (0) for the arbitrary source function (see below) to work

        LeftWidth           = YPosSource(j) - Width(j)/2; %DO NOT EDIT
        RightWidth          = YPosSource(j) + Width(j)/2; %DO NOT EDIT
        
        Time                = 0.1:t/DisTime:t+0.1;                                      % Discritized Time Vector (s) %DO NOT EDIT
        X                   = 0:XPos/DisX:XPos;                                         % Discritized X vector %DO NOT EDIT
        Y                   = -YPos:YPos/DisY*2:YPos;                                   % Discritized Y vector %DO NOT EDIT
        Z                   = 0:ZPos/DisZ:ZPos;                                         % Discritized X vector %DO NOT EDIT
        
        [X1,Y1,Z1,T]        = ndgrid(X,Y,Z,Time);                                  % 4D Grid Construction %DO NOT EDIT
        
        t0                  = 0; %DO NOT EDIT
        X2                  = reshape(X1,numel(X1),1); %DO NOT EDIT
        Y2                  = reshape(Y1,numel(Y1),1); %DO NOT EDIT
        Z2                  = reshape(Z1,numel(Z1),1); %DO NOT EDIT
        T2                  = reshape(T,numel(T),1); %DO NOT EDIT
        
        
        if LayerBottoms(k) > SourceBottom
            Bottom          = SourceBottom; %DO NOT EDIT
        else
            Bottom          = LayerBottoms(k); %DO NOT EDIT
        end
        
        if LayerTops(k) < SourceTop
            Top             = SourceTop; %DO NOT EDIT
        else
            Top             = LayerTops(k); %DO NOT EDIT
        end
        
        
        INT1        = zeros(1,numel(X1)); %DO NOT EDIT
        C1          = zeros(1,numel(X1)); %DO NOT EDIT
        %Calcualates the plume based on the Karanovic et al. (2007) eqn
        for i = 1:numel(X1)
            tildew  = (T2(i)-t0)/2*w; %computs new weights because of the different bounds of integration
            tildef  = (T2(i)-t0)/2*f+(T2(i)+t0)/2; %computs new points because of the different bounds of integration
            INT1(i) = sum(tildew.*(((1./tildef.^(3/2))).*(exp(((ks-lambda).*tildef)-(((X2(i)-Velocity2.*tildef).^2)./(4.*DL.*tildef)))).*...
                ((erfc((Y2(i)-(RightWidth))./(2.*sqrt(DTH.*tildef)))-erfc((Y2(i)-(LeftWidth))./(2.*sqrt(DTH.*tildef))))).*...
                ((erfc((Z2(i)-(Bottom))./(2.*sqrt(DTV.*tildef)))-erfc((Z2(i)-(Top))./(2.*sqrt(DTV.*tildef)))-erfc((Z2(i)+(Bottom))./(2.*sqrt(DTV.*tildef)))+erfc((Z2(i)+(Top))./(2.*sqrt(DTV.*tildef))))))); %Calculates the numerical integration
            C1(i)   = ((Concentration(j).*X2(i))./(8*sqrt(pi*DL)).*exp(-ks.*T2(i))).*INT1(i); %Calculates the plume
        end
        
        C2          =  zeros(size(C1));        
        
        %Inputs for Arbitrary Source Function
            FunctionTime = [0; 2; 4]'; %Row vector for times (yr) when source function changes; include the beginning (0) and end (t)
            FunctionTime = FunctionTime*3.15569E7; %Converts to sec
            SourceFunctionConc = [1; 0.5; 0.5]'; %Row vector for the relative concentration (mg/L) change; start is always 1; include the last concentration at time t
            FunctionConc = SourceFunctionConc*Concentration(j);
            Steepness    = 20; %defines the curvature of the logistic function; DO NOT EDIT
            tX           = FunctionTime(2:end-1);
            ConcX        = [Concentration(j) FunctionConc(2:end-1)'];
            
            for m = 1:length(tX)
                INT2        = zeros(1,numel(X1));
                for i = 1:numel(X1)
                    tildew2 = (T2(i)-tX(m))/2*w; %computs new weights because of the different bounds of integration
                    tildef2 = (T2(i)-tX(m))/2*f+(T2(i)-tX(m))/2; %computs new points because of the different bounds of integration
                    INT2(i) = real(sum(tildew2.*(((1./tildef2.^(3/2))).*(exp(((ks-lambda).*tildef2)-(((X2(i)-Velocity2.*tildef2).^2)./(4.*DL.*tildef2)))).*...
                        ((erfc(real((Y2(i)-(RightWidth))./(2.*sqrt(DTH.*tildef2))))-erfc(real((Y2(i)-(LeftWidth))./(2.*sqrt(DTH.*tildef2)))))).*...
                        ((erfc(real((Z2(i)-(Bottom))./(2.*sqrt(DTV.*tildef2))))-erfc(real((Z2(i)-(Top))./(2.*sqrt(DTV.*tildef2))))-erfc(real((Z2(i)+(Bottom))./(2.*sqrt(DTV.*tildef2))))+erfc(real((Z2(i)+(Top))./(2.*sqrt(DTV.*tildef2)))))))));
                    C2(m,i) = real((1./(1+exp(-Steepness.*(T2(i)-tX(m))))).*(((ConcX(m+1)-ConcX(m)).*(X2(i)))./(8*sqrt(pi*DL)).*exp(-ks.*(T2(i)-tX(m)))).*INT2(i));
                end
            end
            C2(isnan(C2))   = 0; %DO NOT EDIT
            C2              = sum(C2,1); %DO NOT EDIT          
        
        C3                  = C1 + C2; %DO NOT EDIT
        Conc2(k,:)          = C3; %DO NOT EDIT
        size(Conc2);
        Conc2               = sum(Conc2, 1); %DO NOT EDIT
    end
    Conc3(j,:)              = Conc2; %DO NOT EDIT
    size(Conc3);
    Conc2                   = []; %DO NOT EDIT
end
toc
Conc4                       = sum(Conc3, 1); %DO NOT EDIT
C_surf                      = reshape(Conc4,size(X1)); %This is the 4-dimensional data volume for the plume
%%%From here, using C_surf, you can plot any slice of the data volume as
%%%you see fit or represent the plume however you wish.