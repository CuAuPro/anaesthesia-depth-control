                            % MAIN PROGRAM FOR BENCHMARK PATIENT SIMULATOR
                                    % CLARA IONESCU AND DANA COPOT
                                % GHENT UNIVERSITY, BELGIUM, JANUARY 2021

%when this program is used please cite the following: 

% An Open Source Patient Simulator for Design and Evaluation of Computer  Based Multiple Drug Dosing Control for Anesthetic and Hemodynamic Variables
% C. Ionescu, M. Neckebroek, Mi. Ghita and D. Copot
% IEEE Access; 2021; doi 10.1109/ACCESS.2021.3049880  
                                
                                

%% Reset
clear all; clear memory; close all; clc;

%% Parameters
% Sampling time
Ts = 1;

% Constant declaration anesthetic models
gammaP = 2; C50P = 4.16; %Hill function Propofol
gammaR = 2.4; C50R = 8.84; %Hill function Remifentanil
betaBIS = 2.13; gammaBIS = 4; sigmaBIS = 8.2;  E0 = 100; Emax = 100; %Surface interaction model Propofol and Remifentanil
BISdelay = 19.7; 
betaPSI = 2.13; gammaPSI = 4; sigmaPSI = 8.2;
PSIdelay = 19.7; 

% patient database
Patients=[1   74   164   88   1  32.7  60;
          2   67   161   69   1  26.6  53;
          3   75   176   101  1  32.6  69;
          4   69   173   97   1  32.4  67;
          5   45   171   64   1  21.9  52;
          6   57   182   80   1  24.2  62;
          7   74   155   55   1  22.9  44;
          8   71   172   78   1  26.4  60;
          9   65   176   77   1  24.9  60;
          10  72   192   73   1  19.8  62;
          11  69   168   84   1  29.8  60;
          12  60   190   92   1  25.5  71;
          13  61   177   81   1  25.9  62; 
          14  54   173   86   1  27.5  63;
          15  71   172   83   1  28.1  62;
          16  53   186   114  1  33    77;
          17  72   162   87   1  33.2  59;
          18  61   182   93   1  28.1  69;
          19  70   167   77   1  27.6  58;
          20  69   168   82   1  29.1  60;
          21  69   158   81   1  32.4  55;
          22  60   165   85   1  31.2  60;
          23  70   173   69   1  23.1  56;
          24  56   186   99   1  28.6  73];     
[m,n]=size(Patients); 
%%
%%%%%%%select patient%%%%%%%%%
k=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
age    = Patients(k,2); height = Patients(k,3); weight = Patients(k,4); Gender = Patients(k,5); BMI = Patients(k,6); lbm = Patients(k,7);

% Transfer function (PK) from Propofol input to effect site concentration   
V1p = 4.27; V2p = 18.9 - 0.391*(age-53); V3p = 238; %[l]
Cl1p = 1.89 + 0.0456*(weight - 77) - 0.0681*(lbm - 59)+ 0.0264*(height - 177); Cl2p = 1.29 - 0.024*(age - 53); Cl3p = 0.836; % [l.min^(-1)]
k10p = Cl1p/V1p; k12p = Cl2p/V1p; k13p = Cl3p/V1p; k21p = Cl2p/V2p; k31p = Cl3p/V3p; %[min^(-1)]
k1ep = 0.456; ke0p = 0.456; %[min^(-1)]

% Transfer function (PK) from Remifentanil input to effect site concentration
V1r = 5.1 - 0.0201*(age-40) + 0.072*(lbm-55); V2r = 9.82 - 0.0811*(age-40) + 0.108*(lbm-55); V3r = 5.42; %[l]
Cl1r = 2.6 + 0.0162*(age - 40) + 0.0191*(lbm-55); %[l*min^(-1)]
Cl2r = 2.05 - 0.0301*(age - 40); Cl3r = 0.076 - 0.00113*(age - 40); %[l*min^(-1)]
k10r = Cl1r/V1r; k12r = Cl2r/V1r; k13r = Cl3r/V1r; k21r = Cl2r/V2r; k31r = Cl3r/V3r; %[min^(-1)]
ke0r = 0.595 - 0.007*(age - 40); k1er = 0.456; %[min^(-1)]


% Transfer function from Remifentanil effect site concentration to level RASS sedation score
k1r = 0.81; k0r = 0.81;

% Transfer function from Atracurium input to effect site concentration
k1a = 1; k2a = 4; k3a = 10; alphaNMB = 0.0374; gammaN = 2.6677; C50N = 3.2425;

%%% change here if you want to use rocuronium instead of atracurium  %%%%%%% needs to be adapted the values are not yet known
% k1a = 1; k2a = 4; k3a = 10; alphaNMB = 0.0374; gammaN = 2.6677; C50N = 3.2425; 

% Interaction Prop to hemodynamic
%%%============================================
%Interaction model from Propofol to CO
k1Pco = 0.81; k0Pco = 0.81; E0Pco = 5; EmaxPco = 5; gainPco=10; gammaPco = 4.5; C50Pco = 8;

%Interaction model from Propofol to MAP
k1Pmap = 0.61; k0Pmap = 0.81; E0Pmap = 5; EmaxPmap = 5; gainPmap=15; gammaPmap = 4.5; C50Pmap = 6;


% Interaction Remi to hemodynamic
%Interaction model from Remi to CO
k1Rco = 0.51; k0Rco = 0.51; E0Rco = 15; EmaxRco = 5; gainRco=10; gammaRco = 4.5; C50Rco = 12;

%Interaction model from Remi to MAP
k1Rmap = 0.31; k0Rmap = 0.31; E0Rmap = 70; EmaxRmap = 70; gainRmap=10; gammaRmap = 4.5; C50Rmap = 17;

% Hemodynamic model for a nominal patient
K11 = 5; tau11 = 300; T11 = 60; K21 = 12; tau21 = 150; T21 = 50; %dopamine
K11 = 5; tau11 = 300; T11 = 60; K21 = 12; tau21 = 150; T21 = 50; % %dobutamine %% change here if you want to use dobutamine instead of dopamine %%%% needs to be adapted (the TF for this are not yet known)
K12 = 3; tau12 = 40; T12 = 60; K22 = -15; tau22 = 40; T22 = 50; %snp

% Set basis MAP and CO
Cobasis = 5; MAPbasis = 80;

% Constant declaration disturbance models
%Stimulus (disrtubance)
ao=[zeros(1,100) 100*ones(1,50) zeros(1,50) 100*ones(1,50) 100*ones(1,50) zeros(1,100) 100*ones(1,30) zeros(1,50) 100*ones(1,50) zeros(1,200) 100*ones(1,100) zeros(1,200)] ;

%Nociceptor pathway model
INIT = [2*0.3*150 2*0.22*149 150^2 149^2 2*0.08*165 2*0.15*163 165^2 163^2 2*0.1*155 2*0.1*155 155^2 155^2 0.2]; 
Z1 = INIT(1); Z2 = INIT(3); Z3 = INIT(5); Z4 = INIT(7); Z5 = INIT(9); Z6 = INIT(11); %zeros
P1 = INIT(2); P2 = INIT(4); P3 = INIT(6); P4 = INIT(8); P5 = INIT(10); P6 = INIT(12); %poles
K = INIT(13); %gain

%Anesthesiologist in the loop (bolus)
anestS = [zeros(1,80) 10*ones(1,20) zeros(1,60) 10*ones(1,20) zeros(1,40) [0:-1:-10] zeros(1,30) 10*ones(1,10) zeros(1,130) 10*ones(1,20) zeros(1,60) 10*ones(1,20) zeros(1,220) 10*ones(1,20) zeros(1,310)] ;

% Initialise models
anesthesia_model %Initialise anesthetic models
hemodyn_model %Initialise hemodynamical models
dist_model %Initialise disturbance models

% Padé approximation of BIS delay
dBIS = pade(exp(-BISdelay*s),4);

% Padé approximation of PSI delay
dPSI = pade(exp(-PSIdelay*s),4);

% Calculate the gain for RASS
dcgRASS = dcgain(tf(numRass,denRass)*tf(1,[k1r*15 k0r]));

% Approximation of first order plus dead time hemodynamic models
g11 = ((K11)/(1+tau11*s))*(1)/(1+T11*s+(T11*s)^2/2+(T11*s)^3/6 + (T11*s)^4/24);
g12 = ((K12)/(1+tau12*s))*(1)/(1+T12*s+(T12*s)^2/2+(T12*s)^3/6 + (T12*s)^4/24);
g21 = ((K21)/(1+tau21*s))*(1)/(1+T21*s+(T21*s)^2/2+(T21*s)^3/6 + (T21*s)^4/24);
g22 = ((K22)/(1+tau22*s))*(1)/(1+T22*s+(T22*s)^2/2+(T22*s)^3/6 + (T22*s)^4/24);

%% 
% here starts the input changes for open loop simulation/ interaction 
% change one input at a time
Tsim=1200;
Tinject = 10;

CER_REF = 0;
inputForRemi = CER_REF./dcgain(remiSS)*60;

% Input vectors for Propofol, Remifentanil, Atracurium,...
inputStepP=[zeros(1,Tinject) 2*ones(1,Tsim-Tinject)]; %% 
inputStepR=[zeros(1,Tinject) inputForRemi*ones(1,Tsim-Tinject)]; %% 
inputStepA=[zeros(1,Tinject) 0*ones(1,Tsim-Tinject)]; %% 
inputStepD=[zeros(1,Tinject) 0*ones(1,Tsim-Tinject)]; %% 
inputStepS=[zeros(1,Tinject) 0*ones(1,Tsim-Tinject)]; %% 

% define inputs for use in simulink
Propinput=timeseries(inputStepP); Reminput=timeseries(inputStepR); Atracinput=timeseries(inputStepA); Dopinput=timeseries(inputStepD); SNPinput=timeseries(inputStepS); 
anestS=timeseries(anestS); ao=timeseries(ao); 





%% References

BIS_REF = 50;
CEP_REF = 4.16;
%% Matlab-Simulink Exchange Variables

% Hill inverse, parameters for use in Simulink internal models
sigma = sigmaBIS;
gamma = gammaBIS;


% Propofol PK model for internal use in Simulink Model Reference Control
model_prop = {};
propSSdiscrete = c2d(propSS, Ts);
model_prop.A = propSSdiscrete.A;
model_prop.B = propSSdiscrete.B;
model_prop.C = propSSdiscrete.C;
model_prop.D = propSSdiscrete.D;

% model_prop_ss = ss(model_prop.A, ...
%     model_prop.B, ...
%     model_prop.C, ...
%     model_prop.D, ...
%     Ts);

% Remifentanil PK model for internal use in Simulink Model Reference Control
model_remi = {};
remiSSdiscrete = c2d(remiSS, Ts);
model_remi.A = remiSSdiscrete.A;
model_remi.B = remiSSdiscrete.B;
model_remi.C = remiSSdiscrete.C;
model_remi.D = remiSSdiscrete.D;

%% PSO

%Ref model
% propSS poles:
% pole(propSS)
% ans =
%    -0.4560
%    -0.8616
%    -0.0566
%    -0.0025
% -1./pole(propSS)
%    2.1930
%    1.1606
%   17.6637
%  407.2016

%Reference trajectory (first order)
Ar = [0.945]; % Ar = exp(-Ts/Tm) exp(-1/(1/pole))
Br = 1*(eye(size(Ar,1)) - Ar);
Cr = eye(size(Ar,1));
Dr = 0;

H = 10; % Prediction horizont
model_cerefSS = ss(Ar,Br,Cr,Dr, Ts);
dcgain(model_cerefSS)

% Reference model for internal use in Simulink Model Reference Control
model_ceref = {};
model_ceref.A = model_cerefSS.A;
model_ceref.B = model_cerefSS.B;
model_ceref.C = model_cerefSS.C;
model_ceref.D = model_cerefSS.D;

%%
% Parametri
PART_NUM = 50;
max_nr_iter = 20;
w = 1.0; %inertia coefficient
w_damping = 0.99; % inertia damping coefficient
c1 = 3; %social coefficient in own vector
c2 = c1; %social coefficient in cluster

X_max = 6;
X_min = 0;
V_max = (X_max-X_min)/5;
uBest = 0;
nr_inputs = length(uBest);

H = 10; %Horizont

%template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

particle_init = repmat(empty_particle, PART_NUM, 1);

GlobalBest_init.Cost = inf;
GlobalBest_init.Position = 0;

initGlobalBest = false;
for i=1:PART_NUM

    particle_init(i).Position = X_min + (X_max)*rand(1, nr_inputs);
    particle_init(i).Velocity = zeros(nr_inputs,1);

    particle_init(i).Cost = inf;
    particle_init(i).Best.Position = particle_init(i).Position;
    particle_init(i).Best.Cost = particle_init(i).Cost;

    if particle_init(i).Best.Cost < GlobalBest_init.Cost
        GlobalBest_init = particle_init(i).Best;
        initGlobalBest = true;
    end

end
if ~initGlobalBest
    GlobalBest_init = particle_init(i).Best;
    GlobalBest_init.Position = X_min; %Safety...
    initGlobalBest = true;
end


%%
out = sim('Patient_TCI_PSO');
%% Plotting


CEP_REF_plot = CEP_REF*ones(length(out.time), 1);
BIS_REF_plot = ContentrationToBIS(CEP_REF_plot./C50P, CER_REF./C50R, gammaBIS, sigmaBIS, Emax);

figure;

subplot(3,1,1);
plot(out.time./60, BIS_REF_plot, 'k');
title("TCI - PSO", 'Interpreter', 'latex');

hold on;
%plot(out.time./60, out.y_refModel, 'r');
plot(out.time./60, out.BIS_combined, 'b');
%legend('ref', 'ref\_m','y', 'Interpreter', 'latex');
legend('ref','y', 'Interpreter', 'latex');
xlabel("$t [min]$", 'Interpreter', 'latex');
ylabel("$BIS [\%]$", 'Interpreter', 'latex');
ylim([0,Emax*1.05]);

subplot(3,1,2);
plot(out.time./60, CEP_REF_plot, 'k');
hold on;
plot(out.time./60, out.y_refModel, 'r');
plot(out.time./60, out.ce_est, 'b');
legend('$ref$', '$ref_{model}$', '$x_{est}$', 'Interpreter', 'latex');
xlabel("$t [min]$", 'Interpreter', 'latex');
ylabel("$x_e [mg/ml]$", 'Interpreter', 'latex');
%ylim([0, max(xe_ref)*1.05]);
ylim([0, max(out.ce_est)*1.05]);

subplot(3,1,3);
plot(out.time./60, out.u_reg_sat, 'b');
xlabel("$t [min]$", 'Interpreter', 'latex');
ylabel("$u_{Prop} [mg/ml/min]$", 'Interpreter', 'latex');
ylim([-0.5, 6.5]);


