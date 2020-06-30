%% SET PARAMETERS.
% Parameters for Cell Density:
diam = 7;        % Units: um. Average cell diameter.
cpmL = 2*10^6;   % Units: cells/mL. Several parameters below will depend on this cell density. 
pack_eff = cpmL*(((diam/2)^3)*(4/3)*pi)/(10^12);    % Unitless. Converts cell density (cells/mL) to the packing efficiency of cells in space. Pack_eff for randomly packed spheres is 0.64, which back-calculates to 1.22*10^9 cells/mL if cell diameter is 10um, 3.56*10^9 cells/mL if cell diameter is 7um (Fang et al 2013).
cvpc = (1-pack_eff)/pack_eff;       % Unitless. Converts packing efficiency to the number of cell-volumes-worth of extracellular space per cell.
% Parameters for TF Equations:
b = 0.2;    % Units: molecules/hr. Basal rate of TF transcription due to TCR/general stimulation.
p1 = 10;     % Units: molecules/hr. Maximum rate of TF1 transcription due to self-promotion.
p2 = 10;     % Units: molecules/hr. Maximum rate of TF2 transcription due to self-promotion.
P1 = 10;    % Units: molecules. Half-maximal amount of TF1 for self-promotion.
P2 = 10;    % Units: molecules. Half-maximal amount of TF2 for self-promotion.
hp = 2;     % Unitless. Hill exponent for self-promotion.
X1 = 30;     % Units: molecules. Half-maximal amount of TF1 for cross-inhibition.
X2 = 30;     % Units: molecules. Half-maximal amount of TF2 for cross-inhibition.
hx = 2;     % Unitless. Hill exponent for cross-inhibition.
s1 = 20;    % Units: molecules/hr. Maximum rate of TF1 transcription due to CY1 signaling.
s2 = 20;    % Units: molecules/hr. Maximum rate of TF2 transcription due to CY2 signaling.
S1 = 200;    % Units: molecules. Half-maximal amount of CY1 for signaling.
S2 = 200;    % Units: molecules. Half-maximal amount of CY2 for signaling.
hs = 1;     % Unitless. Hill exponent for signaling.
Z1 = 2000;  % Units: molecules. Half-maximal amount of CY1 to inhibit CY2-driven TF2 production. 2400
Z2 = 2000;  % Units: molecules. Half-maximal amount of CY2 to inhibit CY1-driven TF1 production. 1600
hz = 1;     % Unitless. Hill exponent for signaling interference.
dTF1 = 0.15; % Units: /hr. Decay rate of TF1.
dTF2 = 0.15; % Units: /hr. Decay rate of TF2.
nTF = 0.18;    % Unitless. Standard deviation in TF expression as a fraction of the expression rates.
nTF1 = nTF;
nTF2 = nTF;
% Parameters for CY Equations:
a1 = 20*60*60/cvpc;   % Units: molecules/hr. Maximum rate of CY1 production due to TF1 activation, adjusting for cell density.
a2 = 20*60*60/cvpc;   % Units: molecules/hr. Maximum rate of CY2 production due to TF2 activation, adjusting for cell density.
A1 = 80;              % Units: molecules. Half-maximal amount of TF1 for activation.
A2 = 80;              % Units: molecules. Half-maximal amount of TF2 for activation.
ha = 1;               % Unitless. Hill exponent for activation.
R1 = 80;              % Units: molecules. Half-maximal amount of TF1 for repression. 112
R2 = 80;              % Units: molecules. Half-maximal amount of TF2 for repression. 48
hr = 1;               % Unitless. Hill exponent for repression.
U1 = 2000;              % Units: molecules: Half-maximal amount of CY1 for undercutting.
U2 = 2000;              % Units: molecules: Half-maximal amount of CY2 for undercutting.
hu = 1;               % Unitless. Hill exponent for undercutting.
dCY1 = (0.8/cvpc + 0.015);% Units: /hr. Uptake rate of CY1, adjusting for cell density, plus decay rate for CY1. SHOULD BE 0.8/CVPC
dCY2 = (0.8/cvpc + 0.015);% Units: /hr. Uptake rate of CY2, adjusting for cell density, plus decay rate for CY2. SHOULD BE 0.8/CVPC
nCY = 0.18;   % Unitless. Standard deviation in CY expression as a fraction of the expression rates.
nCY1 = nCY;
nCY2 = nCY;
% Extra Sources of Stimulation:
DC1 = 0*(nthroot(49*S1^hs,hs))*0.17;   % Units: molecules. Stimulation through contact with a Type 1 DC. 
                % Fraction of Th1 DCs * value required to make 98% of T Cells express T-bet by probability of DC contact * Frequency of DCs.
DC2 = 0*(nthroot(49*S2^hs,hs))*0.17;   % Units: molecules. Stimulation through contact with a Type 2 DC. 
                % Fraction of Th2 DCs * value required to make 98% of T Cells express GATA3 by probability of DC contact * Frequency of DCs.
EX1 = 0*((1.66054*10^3)*(((diam/2)^3)*(4/3)*pi))/((75*10^3)*(1-pack_eff)); % Units: molecules/cell-volume extracellular space.
    % Converts ng/mL IL-12 (mol weight 75 kDa) to # of IL-12 molecules per cell-volume of extracellular space.
EX2 = 0*((1.66054*10^3)*(((diam/2)^3)*(4/3)*pi))/((13.5*10^3)*(1-pack_eff)); % Units: molecules/cell-volume extracellular space.
    % Converts ng/mL IL-4 (mol weight 13.5 kDa) to # of IL-4 molecules per cell-volume of extracellular space.    
% Non-Dimensional Parameters for Time-Separated TF Subsystem
Asym = (X1*p2)/(X2*p1);
B1 = b/p1;
B2 = b/p2;
H1 = P1/X1;
H2 = P2/X2;
D1 = dTF1*X1/p1;
D2 = dTF2*X2/p2;
F1 = s1/p1;
F2 = s2/p2;
E1 = S1/U1;
E2 = S2/U2;
W1 = Z1/U1;
W2 = Z2/U2;
fM1 = 0; % Fixed M value you're interested in.
fM2 = 0; % Fixed M value you're interested in.
K1 = (F1*(fM1^hs)/(E1^hs + fM1^hs))*(W2^hz/(W2^hz + fM2^hz));
K2 = (F2*(fM2^hs)/(E2^hs + fM2^hs))*(W1^hz/(W1^hz + fM1^hz));
% Non-Dimensional Parameters for Time-Separated CY Subsystem
asym = U1*a2/(U2*a1);
L1 = dCY1*U1/a1;
L2 = dCY2*U2/a2;
V1 = A1/X1;
V2 = A2/X2;
Q1 = R1/X1;
Q2 = R2/X2;
fN1 = 0; % Fixed N value you're interested in.
fN2 = 0; % Fixed N value you're interested in.
C1 = ((fN1^ha)/(V1^ha + fN1^ha))*(Q2^hr/(Q2^hr + fN2^hr));
C2 = ((fN2^ha)/(V2^ha + fN2^ha))*(Q1^hr/(Q1^hr + fN1^hr));
% Adjustments to Non-Dimensional Parameters if Timescales Are NOT Separated.
G1 = X1*a1/(U1*p1);
G2 = X2*a2/(U2*p2);
L1 = dCY1*X1/p1;
L2 = dCY2*X2/p2;
% Parameters for Plotting
color1 = [0.8500, 0.3250, 0.0980];
color2 = [0, 0.4470, 0.7410];
color3 = [0.9722 0.6076 0.3869];
color4 = [0.2778 0.1736 0.1106];
color5 = [0.8715 0.9028 0.9028];
color6 = [0.2917 0.3333 0.4167];

%% FIGURE 2A,B: Distributions of TF Expression in Presence vs. Absence of Cytokines
% PARAMETERS TO ADJUST: a1,2 = 20 or 0, nTF = nCY = 0.18, cpmL = 2*10^6
% CODE TO ADJUST: FIGURE TITLE, HISTOGRAM BINNING, AND PANEL LABEL
F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(DC1+Y(3))^hs)/(S1^hs+(DC1+Y(3))^hs))*(Z2^hz)/(Z2^hz+Y(4)^hz) - dTF1*Y(1);
            (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(DC2+Y(4))^hs)/(S2^hs+(DC2+Y(4))^hs))*(Z1^hz)/(Z1^hz+Y(3)^hz) - dTF2*Y(2);
            (a1*Y(1)^ha)/(A1^ha+Y(1)^ha)*(R2^hr)/(R2^hr+Y(2)^hr)*(U2^hu)/(U2^hu+Y(4)^hu) - dCY1*Y(3);
            (a2*Y(2)^ha)/(A2^ha+Y(2)^ha)*(R1^hr)/(R1^hr+Y(1)^hr)*(U1^hu)/(U1^hu+Y(3)^hu) - dCY2*Y(4)];
G = @(t,Y) [nTF1*Y(1);
            nTF2*Y(2);
            nCY1*Y(3);
            nCY2*Y(4)];
TSPAN = 0:0.1:24;                          
Y0 = [0 0 EX1 EX2];                         
options = sdeset('SDEType','Ito');         
numSims = 1000;
histogramData = zeros(numSims, 1); 
i = 1;
while i <= numSims
    fprintf('Beginning Loop: %d\n', i);
    Y = sde_euler(F,G,TSPAN,Y0,options);
    if (Y(241,1) >= 0) && (Y(241,2) >= 0) && (Y(241,3) >= 0) && (Y(241,4) >= 0)
        histogramData(i) = atan(Y(241,2)/Y(241,1));
        i = i+1
    else
        fprintf('Error!');
    end
end
histogram(histogramData, 10); % 10 for Fig 2a, 20 for Fig 2b
xticks([0 pi/4 pi/2])
xticklabels({'Th1','Mixed','Th2'})
ylim([0 300])
title('Cytokine Signaling Permitted, 1000 SDE Sample Paths','FontSize',14) % Permitted for Fig 2a, Blocked for 2b
xlabel('Balance of Transcription Factor Expression','FontSize',14)
ylabel('Relative Number of Sample Paths','FontSize',14)
yticks([])
text(-0.2,300,'a','FontSize',24); % 'a' for Fig 2a, 'b' for Fig 2b

%% FIGURE 2C (Data): Bifurcation of low density system as cytokines are removed.
numSamples = 20;
avals = [linspace(0,3,2) linspace(3.1,5.9,numSamples-6) linspace(6,20,4)].*(60*60/cvpc);
Results = NaN([numSamples 4 2]);
sol_thresh = 1e-2;
sameness_thresh = 0.5;
eig_thresh = -1e-4;
searchGrid = linspace(0,10,11);
options = optimset('Display','off','MaxFunEvals',100000,'TolFun', 1e-9,'MaxIter',1e4); 
hval = [hp hx hs hz ha hr hu];
for k=1:numSamples
    solutions = zeros(length(searchGrid)^2,5);
    a1 = avals(k);
    a2 = avals(k);
    G1 = X1*a1/(U1*p1);
    G2 = X2*a2/(U2*p2);
    pars = [B1 B2 H1 H2 F1 F2 E1 E2 W1 W2 D1 D2 G1 G2 V1 V2 Q1 Q2 L1 L2 Asym];
    for i = 1:length(searchGrid)
        for j = 1:length(searchGrid)
            startPoint = [searchGrid(i) searchGrid(j) searchGrid(i)/10 searchGrid(j)/10];
            solutions((i-1)*length(searchGrid)+j,1:4) = (fsolve(@(x)root4dSUPP(x,pars,hval),startPoint, options)).^2;
        end
    end
    solutions = unique(round(solutions,4), 'rows');   % Delete duplicate fixpoints.
    TEST = @(x) [(pars(1) + (x(1)^hval(1))/(pars(3)^hval(1)+x(1)^hval(1)))/(1+x(2)^hval(2)) + ((pars(5)*x(3)^hval(3))/(pars(7)^hval(3)+x(3)^hval(3)))*(pars(10)^hval(4))/(pars(10)^hval(4)+x(4)^hval(4)) - pars(11)*x(1);
        ((pars(2) + (x(2)^hval(1))/(pars(4)^hval(1)+x(2)^hval(1)))/(1+x(1)^hval(2)) + ((pars(6)*x(4)^hval(3))/(pars(8)^hval(3)+x(4)^hval(3)))*(pars(9)^hval(4))/(pars(9)^hval(4)+x(3)^hval(4)) - pars(12)*x(2))*pars(21);
        ((pars(13)*x(1)^hval(5))/(pars(15)^hval(5)+x(1)^hval(5)))*((pars(18)^hval(6))/(pars(18)^hval(6)+x(2)^hval(6)))/(1+x(4)^hval(7)) - pars(19)*x(3);
        (((pars(14)*x(2)^hval(5))/(pars(16)^hval(5)+x(2)^hval(5)))*((pars(17)^hval(6))/(pars(17)^hval(6)+x(1)^hval(6)))/(1+x(3)^hval(7)) - pars(20)*x(4))*pars(21)];
    for i = 1:size(solutions,1)   % Test whether candidate fixpoints are truly fixpoints.
        solutions(i,5) = sum(abs(TEST(solutions(i,1:4)))); % Fill the 5th column with a metric of how close to 0 all 4 rates of change really are.
    end
    solutions(solutions(:,5)>=sol_thresh,:) = []; % Remove false candidate fixpoints (those for which all 4 rates of change are not close to 0).
    [C, ia, ic] = unique(round(solutions,1), 'rows'); % Even at this point, there may be some duplicate fixpoints that differ slightly due to numerical approximation.
    for i = 1:length(ia)                                                   % Among sets of as-yet undetected duplicate fixpoints
        solutions(solutions(:,5)~=min(solutions(ic==i,5)) & ic==i, 5) = 1; % mark those that are not closest to total rate-of-change = 0.
    end
    solutions(solutions(:,5)==1,:) = [];
    for i = 1:length(solutions(:,1))        % Now, hopefully, only true and distinct fixpoints remain.
        if (~isreal(nondim4dEIGS(solutions(i,1:4), pars, hval))) || (any(nondim4dEIGS(solutions(i,1:4), pars, hval) >= eig_thresh))
            solutions(i,5) = 1;             % Mark and subsequently remove unstable fixpoints.
        end
    end
    solutions(solutions(:,5)==1,:) = [];
    while sum(pdist(solutions(:,1:4)) < sameness_thresh) > 0           % One last check for duplicate fixpoints
        [ind1, ind2] = find(squareform(pdist(solutions(:,1:4)))<sameness_thresh & squareform(pdist(solutions(:,1:4)))>0, 1);
        cands1 = (solutions(:,5)~=min(solutions([ind1 ind2],5)));
        cands2 = false(length(solutions(:,1)),1);
        cands2([ind1 ind2]) = true;
        if sum(cands1&cands2)==0    % In the rare case that two duplicate fixpoints have identical TEST values, just eliminate the first of the two.
            cands1(find(cands1==0 & cands2==1, 1)) = true;
        end
        solutions(cands1&cands2, :) = [];
    end
    Results(k,:,1) = solutions(1,1:4);
    Results(k,:,2) = solutions(length(solutions(:,1)),1:4);
    readout = ['Progress Made: ', num2str(100*k/numSamples), '%'];
    disp(readout);
end
Fig2cResults = [atan(Results(:,2,1)./Results(:,1,1)) atan(Results(:,2,2)./Results(:,1,2))];
%save('Fig2cResults_FOR_MS2.mat','Fig2cResults');
%load('Fig2cResults_FOR_MS2.mat');

%% FIGURE 2C (Plot)
crv1 = plot(avals./(60*60/cvpc).*5, (pi/4).*ones(1,numSamples), 'LineWidth',8,'Color','k','LineStyle',':');
hold on
crv2 = plot(avals./(60*60/cvpc).*5, Fig2cResults(:,1), 'LineWidth',8,'Color','k');
crv3 = plot(avals./(60*60/cvpc).*5, Fig2cResults(:,2), 'LineWidth',8,'Color','k'); 
title('Equilibrium State at Low Cell Density','FontSize',18)
xlabel('% of Cytokine Secretion Permitted','FontSize',14)
ylabel('Balance of Transcription Factor Expression','FontSize',14)
ylim([0 pi/2])
yticks([0, pi/4, pi/2])
yticklabels({'Th1', 'Mixed', 'Th2'})
text(-16.5,1.7,'c','FontSize',24);
annotation('arrow',[0.15 0.15], [0.55 0.88], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.15 0.15], [0.48 0.15], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.28 0.28], [0.55 0.88], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.28 0.28], [0.48 0.15], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.35 0.35], [0.88 0.55], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.35 0.35], [0.15 0.48], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.517 0.517], [0.88 0.55], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.517 0.517], [0.15 0.48], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.683 0.683], [0.88 0.55], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.683 0.683], [0.15 0.48], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.85 0.85], [0.8 0.55], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.85 0.85], [0.15 0.48], 'Color', [0.75 0.75 0.75]);
legend([crv1 crv3],{'Unstable','Stable'},'FontSize',10)
hold off

%% FIGURE 3A,B (Data): Model vs. Data of TF Expression Time Course.
% PARAMETERS TO ADJUST: Z1=2400, Z2=1600, R1=112, R2=48, cpmL = 2*10^6
numSamples = 1 + 1000;
F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(DC1+Y(3))^hs)/(S1^hs+(DC1+Y(3))^hs))*(Z2^hz)/(Z2^hz+(DC2+Y(4))^hz) - dTF1*Y(1);
            (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(DC2+Y(4))^hs)/(S2^hs+(DC2+Y(4))^hs))*(Z1^hz)/(Z1^hz+(DC1+Y(3))^hz) - dTF2*Y(2);
            (a1*Y(1)^ha)/(A1^ha+Y(1)^ha)*(R2^hr)/(R2^hr+Y(2)^hr)*(U2^hu)/(U2^hu+(DC2+Y(4))^hu) - dCY1*Y(3);
            (a2*Y(2)^ha)/(A2^ha+Y(2)^ha)*(R1^hr)/(R1^hr+Y(1)^hr)*(U1^hu)/(U1^hu+(DC1+Y(3))^hu) - dCY2*Y(4)];
G = @(t,Y) [nTF1*Y(1);
            nTF2*Y(2);
            nCY1*Y(3);
            nCY2*Y(4)];
TSPAN = 0:0.1:96;                        
Y0 = [0 0 0 0];                     
options = sdeset('SDEType','Ito');
Fig3Results = NaN(length(TSPAN), 4, numSamples);
Y = sde_euler(F,G,TSPAN,Y0,options);
Fig3Results(:,:,1) = Y;
nTF1 = 0.18;
nTF2 = 0.18;
nCY1 = 0.18;
nCY2 = 0.18;
G = @(t,Y) [nTF1*Y(1);
            nTF2*Y(2);
            nCY1*Y(3);
            nCY2*Y(4)];
for i = 2:numSamples
    Y = sde_euler(F,G,TSPAN,Y0,options);
    Fig3Results(:,:,i) = Y;
end
nTF1 = 0;
nTF2 = 0;
nCY1 = 0;
nCY2 = 0;
%save('Fig3Results_FOR_MS2.mat','Fig3Results');

%% FIGURE 3A (Plot)
hold on
shp1 = patch('XData', [TSPAN fliplr(TSPAN)], 'YData', [(prctile(Fig3Results(:,1,2:numSamples),25, 3))' fliplr((prctile(Fig3Results(:,1,2:numSamples),75, 3))')], 'FaceColor', color1.^0.5, 'EdgeColor', 'none')
shp2 = patch('XData', [TSPAN fliplr(TSPAN)], 'YData', [(prctile(Fig3Results(:,2,2:numSamples),25, 3))' fliplr((prctile(Fig3Results(:,2,2:numSamples),75, 3))')], 'FaceColor', color2.^0.5, 'EdgeColor', 'none')
alpha(0.5)
crv1 = plot(TSPAN,Fig3Results(:,1,1),'Color',color1,'LineWidth',8)
crv2 = plot(TSPAN,Fig3Results(:,2,1),'Color',color2,'Linewidth',8)
crv3 = plot(TSPAN,mean(Fig3Results(:,1,2:numSamples), 3),'Color', color1.^2,'LineWidth',3)
crv4 = plot(TSPAN,mean(Fig3Results(:,2,2:numSamples), 3),'Color', color2.^2,'LineWidth',3)
TFtimepoints = [0, 6, 12, 16, 24, 36, 48, 72];
TF1data = [0, 1.2, 8, 19, 34, 49, 41, 34];
TF2data = [0, 1, 9, 18.5, 36, 60, 80, 90];
pts1 = plot(TFtimepoints, TF1data, 'o', 'MarkerEdgeColor', color1.^0.5, 'MarkerSize', 8, 'LineWidth', 3)
pts2 = plot(TFtimepoints, TF2data, 'o', 'MarkerEdgeColor', color2.^0.5, 'MarkerSize', 8, 'LineWidth', 3)
title('')
title('Transcription Factor Co-Expression in Vitro','FontSize', 18)
xlabel('Time (h)','FontSize', 14)
ylabel('Copies Per Th Cell', 'FontSize', 14)
legend([crv1 crv2 crv3 crv4 shp1 shp2 pts1 pts2],{'T-bet ODE', 'GATA3 ODE', 'T-bet SDE Mean', 'GATA3 SDE Mean', 'T-bet SDE IQ Range', 'GATA3 SDE IQ Range', 'T-bet Data', 'GATA3 Data'},'FontSize',10,'Location','northwest')
text(-14,140,'a','FontSize',24);
hold off

%% FIGURE 3B (Plot)
% Don't include this figure if you cannot get data from Fang et al authors
% As it is displayed in their paper, the data is too arbtrarily measured.
CYFracs = NaN(length(TSPAN), 2, numSamples);
for i = 1:numSamples
    CYFracs(:,1,i) = ((Fig3Results(:,1,i).^ha)./(A1^ha+Fig3Results(:,1,i).^ha)).*((R2^hr)./(R2^hr+Fig3Results(:,2,i).^hr)).*((U2^hu)./(U2^hu+Fig3Results(:,4,i).^hu));
    CYFracs(:,2,i) = ((Fig3Results(:,2,i).^ha)./(A2^ha+Fig3Results(:,2,i).^ha)).*((R1^hr)./(R1^hr+Fig3Results(:,1,i).^hr)).*((U1^hu)./(U1^hu+Fig3Results(:,3,i).^hu));
end
hold on
shp1 = patch('XData', [TSPAN fliplr(TSPAN)], 'YData', [(prctile(CYFracs(:,1,2:numSamples),25, 3))' fliplr((prctile(CYFracs(:,1,2:numSamples),75, 3))')], 'FaceColor', color1.^0.5, 'EdgeColor', 'none')
shp2 = patch('XData', [TSPAN fliplr(TSPAN)], 'YData', [(prctile(CYFracs(:,2,2:numSamples),25, 3))' fliplr((prctile(CYFracs(:,2,2:numSamples),75, 3))')], 'FaceColor', color2.^0.5, 'EdgeColor', 'none')
alpha(0.5)
crv1 = plot(TSPAN,CYFracs(:,1,1),'Color',color1,'LineWidth',8)
crv2 = plot(TSPAN,CYFracs(:,2,1),'Color',color2,'Linewidth',8)
crv3 = plot(TSPAN,mean(CYFracs(:,1,2:numSamples), 3),'Color', color1.^2,'LineWidth',3)
crv4 = plot(TSPAN,mean(CYFracs(:,2,2:numSamples), 3),'Color', color2.^2,'LineWidth',3)
CYtimepoints = [0, 12, 16, 24, 36, 48, 96];
CY1data = [0, 0.05, 0.075, 0.05, 0.048, 0.046, 0.04];
CY2data = [0, 0.01, 0.02, 0.04, 0.21, 0.23, 0.29];
pts1 = plot(CYtimepoints, CY1data, 'o', 'MarkerEdgeColor', color1.^0.5, 'MarkerSize', 8, 'LineWidth', 3)
pts2 = plot(CYtimepoints, CY2data, 'o', 'MarkerEdgeColor', color2.^0.5, 'MarkerSize', 8, 'LineWidth', 3)
title('Cytokine Co-Expression in Vitro','FontSize', 18)
xlabel('Time (h)','FontSize', 14)
ylabel('Fraction of Th Cells', 'FontSize', 14)
legend([crv1 crv2 crv3 crv4 shp1 shp2 pts1 pts2],{'IFNg ODE', 'IL-4 ODE', 'IFNg SDE Mean', 'IL-4 SDE Mean', 'IFNg SDE IQ Range', 'IL-4 SDE IQ Range', 'IFNg Data', 'IL-4 Data'},'FontSize',10,'Location','northwest')
text(-14,0.5,'b','FontSize',24);
hold off

%% FIGURE 4A (Data): Bifurcation of the Full System with Cell Density. 
numSamples = 100;
Fig4aResults = NaN([numSamples 4 2]);
cellDensities = logspace(6,log10(3*10^9),numSamples);
sol_thresh = 1e-2;
sameness_thresh = 0.5;
searchGrid = linspace(0,10,11);
options = optimset('Display','off'); 
hval = [hp hx hs hz ha hr hu];
for k=1:numSamples
    solutions = zeros(length(searchGrid)^2,5);
    cpmL = cellDensities(k);
    pack_eff = cpmL*(((diam/2)^3)*(4/3)*pi)/(10^12);
    cvpc = (1-pack_eff)/pack_eff;
    a1 = 20*60*60/cvpc;
    a2 = 20*60*60/cvpc;
    dCY1 = 0.8/cvpc + 0.015;
    dCY2 = 0.8/cvpc + 0.015;
    G1 = X1*a1/(U1*p1);
    G2 = X2*a2/(U2*p2);
    L1 = dCY1*X1/p1;
    L2 = dCY2*X2/p2;
    pars = [B1 B2 H1 H2 F1 F2 E1 E2 W1 W2 D1 D2 G1 G2 V1 V2 Q1 Q2 L1 L2 Asym];
    for i = 1:length(searchGrid)
        for j = 1:length(searchGrid)
            startPoint = [searchGrid(i) searchGrid(j)*3 searchGrid(i) searchGrid(j)*3];
            solutions((i-1)*length(searchGrid)+j,1:4) = (fsolve(@(x)root4dSUPP(x,pars,hval),startPoint, options)).^2;
        end
    end
    solutions = unique(round(solutions,4), 'rows');   % Delete duplicate fixpoints.
    TEST = @(x) [(pars(1) + ((x(1)^2)^hval(1))/(pars(3)^hval(1)+(x(1)^2)^hval(1)))/(1+(x(2)^2)^hval(2)) + ((pars(5)*(x(3)^2)^hval(3))/(pars(7)^hval(3)+(x(3)^2)^hval(3)))*(pars(10)^hval(4))/(pars(10)^hval(4)+(x(4)^2)^hval(4)) - pars(11)*(x(1)^2);
        ((pars(2) + ((x(2)^2)^hval(1))/(pars(4)^hval(1)+(x(2)^2)^hval(1)))/(1+(x(1)^2)^hval(2)) + ((pars(6)*(x(4)^2)^hval(3))/(pars(8)^hval(3)+(x(4)^2)^hval(3)))*(pars(9)^hval(4))/(pars(9)^hval(4)+(x(3)^2)^hval(4)) - pars(12)*(x(2)^2))*pars(21);
        ((pars(13)*(x(1)^2)^hval(5))/(pars(15)^hval(5)+(x(1)^2)^hval(5)))*((pars(18)^hval(6))/(pars(18)^hval(6)+(x(2)^2)^hval(6)))/(1+(x(4)^2)^hval(7)) - pars(19)*(x(3)^2);
        (((pars(14)*(x(2)^2)^hval(5))/(pars(16)^hval(5)+(x(2)^2)^hval(5)))*((pars(17)^hval(6))/(pars(17)^hval(6)+(x(1)^2)^hval(6)))/(1+(x(3)^2)^hval(7)) - pars(20)*(x(4)^2))*pars(21)];
    for i = 1:size(solutions,1)   % Test whether candidate fixpoints are truly fixpoints.
        solutions(i,5) = sum(abs(TEST(sqrt(solutions(i,1:4))))); % Fill the 5th column with a metric of how close to 0 all 4 rates of change really are.
    end
    solutions(solutions(:,5)>=sol_thresh,:) = []; % Remove false candidate fixpoints (those for which all 4 rates of change are not close to 0).
    [C, ia, ic] = unique(round(solutions,1), 'rows'); % Even at this point, there may be some duplicate fixpoints that differ slightly due to numerical approximation.
    for i = 1:length(ia)                                                   % Among sets of as-yet undetected duplicate fixpoints
        solutions(solutions(:,5)~=min(solutions(ic==i,5)) & ic==i, 5) = 1; % mark those that are not closest to total rate-of-change = 0.
    end
    solutions(solutions(:,5)==1,:) = [];
    for i = 1:length(solutions(:,1))        % Now, hopefully, only true and distinct fixpoints remain.
        if (~isreal(nondim4dEIGS(solutions(i,1:4), pars, hval))) || (any(nondim4dEIGS(solutions(i,1:4), pars, hval) >= 0))
            solutions(i,5) = 1;             % Mark and subsequently remove unstable fixpoints.
        end
    end
    solutions(solutions(:,5)==1,:) = [];
    while sum(pdist(solutions(:,1:4)) < sameness_thresh) > 0           % One last check for duplicate fixpoints
        [ind1, ind2] = find(squareform(pdist(solutions(:,1:4)))<sameness_thresh & squareform(pdist(solutions(:,1:4)))>0, 1);
        cands1 = (solutions(:,5)~=min(solutions([ind1 ind2],5)));
        cands2 = false(length(solutions(:,1)),1);
        cands2([ind1 ind2]) = true;
        if sum(cands1&cands2)==0    % In the rare case that two duplicate fixpoints have identical TEST values, just eliminate the first of the two.
            cands1(find(cands1==0 & cands2==1, 1)) = true;
        end
        solutions(cands1&cands2, :) = [];
    end
    Fig4aResults(k,:,1) = solutions(1,1:4);
    Fig4aResults(k,:,2) = solutions(length(solutions(:,1)),1:4);
    readout = ['Progress Made: ', num2str(100*k/numSamples), '%'];
    disp(readout);
end
%save('Fig4aResults_FOR_MS2', 'Fig4aResults');
%load('Fig4aResults_FOR_MS2');

%% FIGURE 4B (Data): Cell Density Affects Cytokine Production:Removal Ratio
numSamples = 100;
cellDensities = logspace(6,log10(3*10^9),numSamples);
cvpcs = (1-(cellDensities.*(((diam/2)^3).*(4/3).*pi)./(10^12)))./(cellDensities.*(((diam/2)^3).*(4/3).*pi)./(10^12));
avals = ((20*60*60)./cvpcs);
dvals = (0.8./cvpcs + 0.015);
ratios = avals./dvals;

%% Figure 4C (Data): Cytokine Production:Removal Ratio Impacts CY:TF Ratio at the Mixed Equilibrium
Fig4cResults = NaN(4,length(cellDensities));
for i = 1:length(cellDensities)
    F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(DC1+Y(3))^hs)/(S1^hs+(DC1+Y(3))^hs))*(Z2^hz)/(Z2^hz+(DC2+Y(4))^hz) - dTF1*Y(1);
        (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(DC2+Y(4))^hs)/(S2^hs+(DC2+Y(4))^hs))*(Z1^hz)/(Z1^hz+(DC1+Y(3))^hz) - dTF2*Y(2);
        (avals(i)*Y(1)^ha)/(A1^ha+Y(1)^ha)*(R2^hr)/(R2^hr+Y(2)^hr)*(U2^hu)/(U2^hu+(DC2+Y(4))^hu) - dvals(i)*Y(3);
        (avals(i)*Y(2)^ha)/(A2^ha+Y(2)^ha)*(R1^hr)/(R1^hr+Y(1)^hr)*(U1^hu)/(U1^hu+(DC1+Y(3))^hu) - dvals(i)*Y(4)];
    G = @(t,Y) [nTF1*Y(1);
        nTF2*Y(2);
        nCY1*Y(3);
        nCY2*Y(4)];
    TSPAN = 0:0.1:960;                          
    Y0 = [0 0 0 0];                       
    options = sdeset('SDEType','Ito');          
    Y = sde_euler(F,G,TSPAN,Y0,options);        
    Fig4cResults(:,i) = (Y(9601,:))';
end
ratios2 = Fig4cResults(3,:)./Fig4cResults(1,:);

%% FIGURE 4D (Data): Strength of Stabilizing vs. Destabilizing Interactions with [CY]:[TF] Ratio
Gvals = (X1/(U1*p1)).*avals;
Lvals = (X1/p1).*dvals;
Nstar = Fig4cResults(1,:)./X1;
Mstar = Fig4cResults(3,:)./U1;
jacobian = NaN(4, 4, length(cellDensities));
for i = 1:length(cellDensities)
    Pars = [B1 B2 H1 H2 F1 F2 E1 E2 W1 W2 D1 D2 Gvals(i) Gvals(i) V1 V2 Q1 Q2 Lvals(i) Lvals(i) Asym];
    Hval = [hp hx hs hz ha hr hu];
    
    jacobian(:,:,i) = [((1/(1+Nstar(i)^Hval(2))))*((Pars(3)^Hval(1))*Hval(1)*(Nstar(i)^(Hval(1)-1)))/((Pars(3)^Hval(1)+Nstar(i)^Hval(1))^2)-Pars(11),...
    -(Pars(1)+(Nstar(i)^Hval(1))/(Pars(3)^Hval(1)+Nstar(i)^Hval(1)))*((Hval(2)*(Nstar(i))^(Hval(2)-1))/((1+Nstar(i)^Hval(2))^2)),...
    ((Pars(5)*Pars(10)^Hval(4))/(Pars(10)^Hval(4)+Mstar(i)^Hval(4)))*(((Pars(7)^Hval(3))*Hval(3)*Mstar(i)^(Hval(3)-1))/((Pars(7)^Hval(3)+Mstar(i)^Hval(3))^2)),...
    -((Pars(5)*Mstar(i)^Hval(3))/(Pars(7)^Hval(3)+Mstar(i)^Hval(3)))*(((Pars(10)^Hval(4))*Hval(4)*Mstar(i)^(Hval(4)-1))/((Pars(10)^Hval(4)+Mstar(i)^Hval(4))^2));...
        
    (-(Pars(2)+(Nstar(i)^Hval(1))/(Pars(4)^Hval(1)+Nstar(i)^Hval(1)))*((Hval(2)*(Nstar(i))^(Hval(2)-1))/((1+Nstar(i)^Hval(2))^2)))*Pars(21),...
    (((1/(1+Nstar(i)^Hval(2))))*((Pars(4)^Hval(1))*Hval(1)*(Nstar(i)^(Hval(1)-1)))/((Pars(4)^Hval(1)+Nstar(i)^Hval(1))^2)-Pars(12))*Pars(21),...
    (-((Pars(6)*Mstar(i)^Hval(3))/(Pars(8)^Hval(3)+Mstar(i)^Hval(3)))*(((Pars(9)^Hval(4))*Hval(4)*Mstar(i)^(Hval(4)-1))/((Pars(9)^Hval(4)+Mstar(i)^Hval(4))^2)))*Pars(21),...
    (((Pars(6)*Pars(9)^Hval(4))/(Pars(9)^Hval(4)+Mstar(i)^Hval(4)))*(((Pars(8)^Hval(3))*Hval(3)*Mstar(i)^(Hval(3)-1))/((Pars(8)^Hval(3)+Mstar(i)^Hval(3))^2)))*Pars(21);...
    
    ((Pars(13)*Pars(18)^Hval(6))/(Pars(18)^Hval(6)+Nstar(i)^Hval(6)))*(1/(1+Mstar(i)^Hval(7)))*((Hval(5)*(Nstar(i)^(Hval(5)-1))*Pars(15)^Hval(5))/((Pars(15)^Hval(5)+Nstar(i)^Hval(5))^2)),...
    -((Nstar(i)^Hval(5))/(Pars(15)^Hval(5)+Nstar(i)^Hval(5)))*(1/(1+Mstar(i)^Hval(7)))*((Pars(13)*(Pars(18)^Hval(6))*Hval(6)*Nstar(i)^(Hval(6)-1))/((Pars(18)^Hval(6)+Nstar(i)^Hval(6))^2)),...
    -Pars(19),...
    -((Pars(13)*Nstar(i)^Hval(5))/(Pars(15)^Hval(5)+Nstar(i)^Hval(5)))*((Pars(18)^Hval(6))/(Pars(18)^Hval(6)+Nstar(i)^Hval(6)))*((Hval(7)*Mstar(i)^(Hval(7)-1))/((1+Mstar(i)^Hval(7))^2));...
    
    (-((Nstar(i)^Hval(5))/(Pars(16)^Hval(5)+Nstar(i)^Hval(5)))*(1/(1+Mstar(i)^Hval(7)))*((Pars(14)*(Pars(17)^Hval(6))*Hval(6)*Nstar(i)^(Hval(6)-1))/((Pars(17)^Hval(6)+Nstar(i)^Hval(6))^2)))*Pars(21),...
    (((Pars(14)*Pars(17)^Hval(6))/(Pars(17)^Hval(6)+Nstar(i)^Hval(6)))*(1/(1+Mstar(i)^Hval(7)))*((Hval(5)*(Nstar(i)^(Hval(5)-1))*Pars(16)^Hval(5))/((Pars(16)^Hval(5)+Nstar(i)^Hval(5))^2)))*Pars(21),...
    (-((Pars(14)*Nstar(i)^Hval(5))/(Pars(16)^Hval(5)+Nstar(i)^Hval(5)))*((Pars(17)^Hval(6))/(Pars(17)^Hval(6)+Nstar(i)^Hval(6)))*((Hval(7)*Mstar(i)^(Hval(7)-1))/((1+Mstar(i)^Hval(7))^2)))*Pars(21),...
    (-Pars(20))*Pars(21)];
end
bminusa = (abs(squeeze(jacobian(1,2,:))) - abs(squeeze(jacobian(1,1,:))));
hminusg = (abs(squeeze(jacobian(3,4,:))) - abs(squeeze(jacobian(3,3,:))));
cplusd = (abs(squeeze(jacobian(1,3,:))) + abs(squeeze(jacobian(1,4,:))));
eplusf = (abs(squeeze(jacobian(3,1,:))) + abs(squeeze(jacobian(3,2,:))));

%% FIGURE 4A-D (Plot)
plotData1 = atan(Fig4aResults(:,1,1)./Fig4aResults(:,2,1));
plotData2 = atan(Fig4aResults(:,1,2)./Fig4aResults(:,2,2));
plotData3 = (repelem(pi/4, numSamples))';

tiles = tiledlayout(2,2);
tiles.Padding = 'none';
tiles.TileSpacing = 'none';
nexttile
crv1 = semilogx(cellDensities, plotData1, 'LineWidth',6,'Color','k');
hold on
crv2 = semilogx(cellDensities, plotData2, 'LineWidth',6,'Color','k');
crv3 = semilogx(cellDensities, plotData3, ':', 'LineWidth',6,'Color','k');
ylim([0 pi/2])
yticks([0 pi/4 pi/2])
yticklabels({'Th1','Mixed','Th2'})
title('Equilibrium States Across Cell Densities', 'FontSize', 10)
%xlabel('Cell Density (cells/mL)','FontSize',10)
xlim([10^6 3*10^9])
ylabel('Balance of Transcription Factor Expression','FontSize', 10)
legend([crv2 crv3],{'Stable','Unstable'}, 'FontSize', 8, 'Location', 'northeast');
hold off

nexttile
crv1 = semilogx(ratios2, bminusa.*hminusg, 'LineWidth', 6, 'Color', 'k');
hold on
crv2 = semilogx(ratios2, cplusd.*eplusf, 'LineWidth', 6, 'Color', 'k', 'LineStyle',':');
set(gca,'YTick', [])
title('Within-Scale vs. Cross-Scale Interactions','FontSize',10);
%xlabel('Ratio of [Cytokine]:[Transcription Factor] at Mixed Equilibrium','FontSize',10);
ylabel('Cumulative Strength of Interactions','FontSize',10);
xlim([3 100]);
ylim([0 0.04]);
legend([crv1 crv2 ],{'Within-Scale Net Stabilization','Cross-Scale Total Destabilization'}, 'FontSize', 8, 'Location', 'northwest');
hold off

nexttile
loglog(cellDensities, ratios, 'LineWidth', 6, 'Color', 'k');
ylabel('Ratio of Cytokine Production:Removal Rates', 'FontSize', 10);
xlabel('Cell Density (cells/mL)', 'FontSize', 10);
xlim([10^6 3*10^9]);
ylim([9*10^2 10^5]);

nexttile
loglog(ratios2, ratios, 'LineWidth', 6, 'Color', 'k');
%ylabel('Ratio of Cytokine Production:Removal Rates', 'FontSize', 10);
xlabel('Ratio of [Cytokine]:[Transcription Factor] at Mixed Equilibrium','FontSize',10);
xlim([3 100]);
ylim([9*10^2 10^5]);

annotation('ellipse',[0.19 0.737 0.03 0.04], 'Color', 'r', 'LineWidth',5)
annotation('ellipse',[0.757 0.639 0.03 0.04], 'Color', 'r', 'LineWidth',5)
annotation('arrow',[0.10 0.10], [0.955 0.77], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.17 0.17], [0.955 0.77], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.24 0.24], [0.77 0.955], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.31 0.31], [0.77 0.955], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.38 0.38], [0.77 0.955], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.10 0.10], [0.56 0.745], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.17 0.17], [0.56 0.745], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.24 0.24], [0.745 0.56], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.31 0.31], [0.745 0.56], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.38 0.38], [0.745 0.56], 'Color', [0.75 0.75 0.75]);
annotation('arrow',[0.205 0.205], [0.737 0.278], 'Color', 'r', 'LineWidth',5, 'HeadWidth', 15);
annotation('arrow',[0.21 0.767], [0.278 0.278], 'Color', 'r', 'LineWidth',5, 'HeadWidth', 15);
annotation('arrow',[0.772 0.772], [0.283 0.637], 'Color', 'r', 'LineWidth',5, 'HeadWidth', 15);
annotation('textbox',[0.01 0.96 0.05 0.05],'String','a','FitBoxToText','on','FontSize',12,'EdgeColor','none');
annotation('textbox',[0.01 0.485 0.05 0.05],'String','b','FitBoxToText','on','FontSize',12,'EdgeColor','none');
annotation('textbox',[0.52 0.96 0.05 0.05],'String','d','FitBoxToText','on','FontSize',12,'EdgeColor','none');
annotation('textbox',[0.52 0.49 0.05 0.05],'String','c','FitBoxToText','on','FontSize',12,'EdgeColor','none');

%% FIGURE 5A (Data): QS in vivo with APCs Are Present
% Set Parameters: Cell Density = 10^9 cells/mL.
TSPAN = linspace(0,72,721);
early = 0.003;
late = 0.17;
DC1data = ones(6,length(TSPAN)).*(early*(nthroot(49*S1,hs))).*0.6;
DC2data = ones(6,length(TSPAN)).*(early*(nthroot(49*S2,hs))).*0.4;
for i = 1:5
    DC1data(i,i*60+1:length(TSPAN)) = zeros(1,length(TSPAN)-(i*60));
    DC2data(i,i*60+1:length(TSPAN)) = ones(1,length(TSPAN)-(i*60)).*(late*(nthroot(49*S2,hs)));
end
Fig5Results = NaN(6,length(TSPAN),2);
for i = 1:6
    dc1 = @(t) interp1(TSPAN,DC1data(i,:),t);
    dc2 = @(t) interp1(TSPAN,DC2data(i,:),t);
    F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(Y(3)+dc1(t))^hs)/(S1^hs+(Y(3)+dc1(t))^hs))*((Z2^hz)/(Z2^hz+(Y(4)+dc2(t))^hz)) - dTF1*Y(1);
        (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(Y(4)+dc2(t))^hs)/(S2^hs+(Y(4)+dc2(t))^hs))*((Z1^hz)/(Z1^hz+(Y(3)+dc1(t))^hz)) - dTF2*Y(2);
        ((a1*Y(1)^ha)/(A1^ha+Y(1)^ha))*((R2^hr)/(R2^hr+Y(2)^hr))*((U2^hu)/(U2^hu+(Y(4)+dc2(t))^hu)) - dCY1*Y(3);
        ((a2*Y(2)^ha)/(A2^ha+Y(2)^ha))*((R1^hr)/(R1^hr+Y(1)^hr))*((U1^hu)/(U1^hu+(Y(3)+dc1(t))^hu)) - dCY2*Y(4)];
    G = @(t,Y) [nTF1*Y(1);
        nTF2*Y(2);
        nCY1*Y(3);
        nCY2*Y(4)];
    Y0 = [0 0 0 0];
    options = sdeset('SDEType','Ito');
    Y = sde_euler(F,G,TSPAN,Y0,options);
%     Fig5Results(i,:,1) = Y(:,1);
%     Fig5Results(i,:,2) = Y(:,2);
    Fig5Results(i,:,1) = ((Y(:,1).^ha)./(A1^ha+Y(:,1).^ha)).*((R2^hr)./(R2^hr+Y(:,2).^hr)).*((U2^hu)./(U2^hu+Y(:,4).^hu));
    Fig5Results(i,:,2) = ((Y(:,2).^ha)./(A2^ha+Y(:,2).^ha)).*((R1^hr)./(R1^hr+Y(:,1).^hr)).*((U1^hu)./(U1^hu+Y(:,3).^hu));
end
%save('Fig5Results_FOR_MS2', 'Fig5Results');

%% FIGURE 5A (Plot)
figure()
hold on
ln1 = line([6 6], [0 pi/2], 'Color',[0 0 0], 'LineWidth', 2, 'LineStyle','--')
ln2 = line([12 12], [0 pi/2], 'Color',[0.2, 0.2, 0.2],'LineWidth', 2, 'LineStyle','--')
ln3 = line([18 18], [0 pi/2], 'Color',[0.4, 0.4, 0.4],'LineWidth', 2, 'LineStyle','--')
ln4 = line([24 24], [0 pi/2], 'Color',[0.6, 0.6, 0.6],'LineWidth', 2, 'LineStyle','--')
ln5 = line([30 30], [0 pi/2], 'Color',[0.8, 0.8, 0.8],'LineWidth', 2, 'LineStyle','--')
data1 = atan(Fig5Results(1,:,2)./Fig5Results(1,:,1));
data2 = atan(Fig5Results(2,:,2)./Fig5Results(2,:,1));
data3 = atan(Fig5Results(3,:,2)./Fig5Results(3,:,1));
data4 = atan(Fig5Results(4,:,2)./Fig5Results(4,:,1));
data5 = atan(Fig5Results(5,:,2)./Fig5Results(5,:,1));
data6 = atan(Fig5Results(6,:,2)./Fig5Results(6,:,1));
data1(isnan(data1)) = pi/4;
data2(isnan(data2)) = pi/4;
data3(isnan(data3)) = pi/4;
data4(isnan(data4)) = pi/4;
data5(isnan(data5)) = pi/4;
data6(isnan(data6)) = pi/4;
crv1 = plot(TSPAN,data1,'Color',[0 0 0],'LineWidth',7) 
crv2 = plot(TSPAN,data2,'Color',[0.2, 0.2, 0.2], 'LineWidth',6)
crv3 = plot(TSPAN,data3,'Color',[0.4, 0.4, 0.4], 'LineWidth',5)
crv4 = plot(TSPAN,data4,'Color',[0.6, 0.6, 0.6], 'LineWidth',4)
crv5 = plot(TSPAN,data5,'Color',[0.8, 0.8, 0.8], 'LineWidth',3)
crv6 = plot(TSPAN,data6,'Color', [0,0,0], 'LineWidth',2)
title('APC Instruction Transitions to Th Quorum Sensing','FontSize',18)
xlabel('Time (h)', 'FontSize', 14)
xlim([0 72])
ylabel('Balance of Cytokine Expression','FontSize',14)
yticks([0, pi/4, pi/2])
yticklabels({'Th1', 'Mixed', 'Th2'})
legend([crv1 crv2 crv3 crv4 crv5 crv6], {'Switch at 6 hr','Switch at 12 hr', 'Switch at 18 hr', 'Switch at 24 hr', 'Switch at 30 hr', 'No Switch'})
annotation('textbox', [0.04 0.5 0.5 0.5], 'String', 'a', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
hold off

%% FIGURE 5B (Data): Speed of Quorum Formation Depends on Frequency and Polarization of DCs
TSPAN = linspace(0,100,1001);
flipTime = linspace(0.5,100,200);
DCfreq = logspace(log10(0.003), log10(0.17), 10);
DCbias = logspace(log10(0.6), log10(1), 10);
Fig5bResults = NaN(length(DCfreq),length(DCbias));
for i = 1:length(DCfreq)
    for j = 1:length(DCbias)
        DC1data = zeros(1,length(TSPAN)).*(nthroot(49*S1,hs));
        DC2data = ones(1,length(TSPAN)).*(nthroot(49*S2,hs)).*0.17;
        k = 1;
        while k <= length(flipTime)
            DC1data(1:(k*5+1)) = ones(1,(k*5+1)).*DCfreq(i).*DCbias(j).*(nthroot(49*S1,hs));
            DC2data(1:(k*5+1)) = ones(1,(k*5+1)).*DCfreq(i).*(1-DCbias(j)).*(nthroot(49*S2,hs));
            dc1 = @(t) interp1(TSPAN,DC1data,t);
            dc2 = @(t) interp1(TSPAN,DC2data,t);
            F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(Y(3)+dc1(t))^hs)/(S1^hs+(Y(3)+dc1(t))^hs))*((Z2^hz)/(Z2^hz+(Y(4)+dc2(t))^hz)) - dTF1*Y(1);
                (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(Y(4)+dc2(t))^hs)/(S2^hs+(Y(4)+dc2(t))^hs))*((Z1^hz)/(Z1^hz+(Y(3)+dc1(t))^hz)) - dTF2*Y(2);
                ((a1*Y(1)^ha)/(A1^ha+Y(1)^ha))*((R2^hr)/(R2^hr+Y(2)^hr))*((U2^hu)/(U2^hu+(Y(4)+dc2(t))^hu)) - dCY1*Y(3);
                ((a2*Y(2)^ha)/(A2^ha+Y(2)^ha))*((R1^hr)/(R1^hr+Y(1)^hr))*((U1^hu)/(U1^hu+(Y(3)+dc1(t))^hu)) - dCY2*Y(4)];
            G = @(t,Y) [nTF1*Y(1);
                nTF2*Y(2);
                nCY1*Y(3);
                nCY2*Y(4)];
            Y0 = [0 0 0 0];
            options = sdeset('SDEType','Ito');
            Y = sde_euler(F,G,TSPAN,Y0,options);
            CY1Frac = ((Y(:,1).^ha)./(A1^ha+Y(:,1).^ha)).*((R2^hr)./(R2^hr+Y(:,2).^hr)).*((U2^hu)./(U2^hu+Y(:,4).^hu));
            CY2Frac = ((Y(:,2).^ha)./(A2^ha+Y(:,2).^ha)).*((R1^hr)./(R1^hr+Y(:,1).^hr)).*((U1^hu)./(U1^hu+Y(:,3).^hu));
            CYBalance = (atan(CY2Frac./CY1Frac))./(pi/2);
            CYBalance(1) = 0.5;
            if (CYBalance(length(TSPAN)) > 0.5) || ((CYBalance(length(TSPAN))-CYBalance(length(TSPAN)-10)) > 0.001)
                k = k+1;
            else
                Fig5bResults(i,j) = k/2;
                k = length(flipTime) + 1;
            end
        end
        fprintf('Finished Combo %d \n', (i-1)*length(DCfreq)+j);
    end
end
%save('Fig5bResults_FOR_MS2', 'Fig5bResults');
%load('Fig5bResults_FOR_MS2');

%% FIGURE 5B (Plot)
hmData = Fig5bResults - (min(min(Fig5bResults)));
hmData = hmData./(max(max(hmData)));
hmData = flipud(hmData);
h = heatmap(100.*round(DCbias,2), round(logspace(3,1,10), 0), hmData, 'Colormap', jet, 'ColorbarVisible','off','CellLabelColor','none','ColorScaling','log');
h.XLabel = 'Effector Bias of APCs (%)';
h.YLabel = 'Number of APCs';
annotation('textbox', [0.025 0.5 0.5 0.5], 'String', 'b', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);

%% FIGURE 5B (Colorbar)
axis off
colormap(jet);
cb = colorbar('Ticks',linspace(0,1,5), 'TickLabels',round(logspace(log10(1),log10(21), 5),0));
ylabel(cb,'Time Until Quorum Formation (h)','FontSize',16)

%% FIGURE 6A-D (Data): Representative APC Time Course and SDE Sample Paths with Different Noise Levels
% Generate APC time-courses, and save an illustrative one.
% Set parameters appropriately: cell density = 10^9 cells/mL.
DCfreq = 0.17; % Probability of T cell encountering a DC. Constant.
threshold = [20; 40; 60; 80; 100; 120; 140; 160; 180; 200; 220; 240]; % Number of timesteps of 100%-biased DCs to opposite type for optimum to switch.
timecourselength = 7*24*10+1; % Length of time course in time steps.
DC1bias = nan(1,timecourselength); % Initialize the DC1 bias time course.
DC1bias(1) = binornd(1,0.5); % First time step is either 0 or 1.
count = 1; % This counts which time step we're on.
streak = poissrnd(120); % How many time steps the initial DC bias will last. 
DC1bias(count+1:count+streak) = DC1bias(1,1); % Set those first X time steps.
Th1response = repmat(DC1bias, [length(threshold),1]); % Initialize the "optimal" Th1 cell response. So far, it is the same as the DC1 bias time course.
flag = (ones(length(threshold),1)).*(count+streak);
curropt = Th1response(:,1); % This keeps track of the current optimal Th1 cell response, for each threshold value.
integral = zeros(length(threshold),1); % This counts the total integral that might cause a switch in the optimal response, for each threshold value.
count = count + streak; % Update the most recent time step of the DC1 bias time course to have been filled in.
while (count < timecourselength)
    new = round((betarnd(DC1bias(count)+0.05, 1.05-DC1bias(count))), 1);
    streak = min(poissrnd(100),timecourselength-count);
    DC1bias(count+1:count+streak) = new;
    Chal = (abs(new-curropt) > 0.5);
    nChalUpdate = repmat(curropt,[1,timecourselength]); 
    rowixn = find(~Chal); % For-loop matrix update code from Sean de Wolski at Matlab.
    for ii = 1:numel(rowixn) % For-loop matrix update code from Sean de Wolski at Matlab.
        rowiin = rowixn(ii); % For-loop matrix update code from Sean de Wolski at Matlab.
        colixn = flag(rowiin):(count+streak); % For-loop matrix update code from Sean de Wolski at Matlab.
        Th1response(rowiin, colixn) = nChalUpdate(rowiin, colixn); % For-loop matrix update code from Sean de Wolski at Matlab.
    end % For-loop matrix update code from Sean de Wolski at Matlab.
    flag(~Chal) = count + streak;
    integral(~Chal) = 0;
    yChalInt = (abs(new-curropt)).*streak;
    integral(Chal) = integral(Chal) + yChalInt(Chal);
    Cros = (Chal & (integral>=threshold));
    curropt(Chal & Cros) = 1 - curropt(Chal & Cros);
    yChalUpdate = repmat(curropt,[1,timecourselength]); % Start edits after this line
    rowixy = find(Chal & Cros); % For-loop matrix update code from Sean de Wolski at Matlab.
    for ii = 1:numel(rowixy) % For-loop matrix update code from Sean de Wolski at Matlab.
        rowiiy = rowixy(ii); % For-loop matrix update code from Sean de Wolski at Matlab.
        colixy = flag(rowiiy):(count+streak); % For-loop matrix update code from Sean de Wolski at Matlab.
        Th1response(rowiiy, colixy) = yChalUpdate(rowiiy, colixy); % For-loop matrix update code from Sean de Wolski at Matlab.
    end % For-loop matrix update code from Sean de Wolski at Matlab.
    integral(Chal & Cros) = 0;
    flag(Chal & Cros) = count + streak;
    if (count+streak >= timecourselength)
        Th1response((Chal & ~Cros), flag(Chal & ~Cros):(count+streak)) = yChalUpdate((Chal & ~Cros), flag(Chal & ~Cros):(count+streak));
    end
    count = count + streak;
end
DC2bias = 1-DC1bias;
hold on
crv1 = plot(1:timecourselength, DC1bias,'Color','k','LineWidth',8);
crv2 = plot(1:timecourselength, Th1response(1,:),'Color','r','LineWidth',6);
crv3 = plot(1:timecourselength, Th1response(6,:),'Color','c','LineWidth',4);
crv4 = plot(1:timecourselength, Th1response(12,:),'Color','g','LineWidth',2);
ylim([0 1]);
hold off
Fig6adResults = [DC2bias; Th1response];
%save('Fig6adResults_FOR_MS2','Fig6adResults');
%load('Fig6adResults_FOR_MS2');

%% FIGURE 6A (Data): Representative APC Time Course and ODE Trajectory
TSPAN = 0:0.1:168;  
DCfreq = 0.17;
DC1tc = @(t) interp1(TSPAN, ((1-Fig6adResults(1,:)).*DCfreq.*(nthroot(49*S1,hs))), t);   
DC2tc = @(t) interp1(TSPAN, (Fig6adResults(1,:).*DCfreq.*(nthroot(49*S1,hs))), t);
F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(Y(3)+DC1tc(t))^hs)/(S1^hs+(Y(3)+DC1tc(t))^hs))*((Z2^hz)/(Z2^hz+(Y(4)+DC2tc(t))^hz)) - dTF1*Y(1);
    (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(Y(4)+DC2tc(t))^hs)/(S2^hs+(Y(4)+DC2tc(t))^hs))*((Z1^hz)/(Z1^hz+(Y(3)+DC1tc(t))^hz)) - dTF2*Y(2);
    ((a1*Y(1)^ha)/(A1^ha+Y(1)^ha))*((R2^hr)/(R2^hr+Y(2)^hr))*((U2^hu)/(U2^hu+(Y(4)+DC2tc(t))^hu)) - dCY1*Y(3);
    ((a2*Y(2)^ha)/(A2^ha+Y(2)^ha))*((R1^hr)/(R1^hr+Y(1)^hr))*((U1^hu)/(U1^hu+(Y(3)+DC1tc(t))^hu)) - dCY2*Y(4)];
G = @(t,Y) [nTF1*Y(1);
    nTF2*Y(2);
    nCY1*Y(3);
    nCY2*Y(4)];                         
Y0 = [0 0 0 0];                        
options = sdeset('SDEType','Ito');          
Y = sde_euler(F,G,TSPAN,Y0,options);       
CY1Fracs = ((Y(:,1).^ha)./(A1^ha+Y(:,1).^ha)).*((R2^hr)./(R2^hr+Y(:,2).^hr)).*((U2^hu)./(U2^hu+Y(:,4).^hu));
CY2Fracs = ((Y(:,2).^ha)./(A2^ha+Y(:,2).^ha)).*((R1^hr)./(R1^hr+Y(:,1).^hr)).*((U1^hu)./(U1^hu+Y(:,3).^hu));
Fig6aResults = (atan(CY2Fracs./CY1Fracs))./(pi/2);
Fig6aResults(1) = 0.5;
%save('Fig6aResults_FOR_MS2','Fig6aResults');

%% FIGURE 6B-D (Data): Representative APC Time Course and SDE Sample Paths with Different Noise Levels
% Now compare SDE sample paths to the optimal Th responses. Vary noise.
% Set parameters appropriately: nTF = CY = 0.2, 0.6, 1.
% Also, change b vs c vs d based on noise level (2nd, 2nd-to-last lines)
numSims = 20;
Fig6dResults = nan(numSims,length(TSPAN));
TSPAN = 0:0.1:168; 
DC1tc = @(t) interp1(TSPAN, ((1-Fig6adResults(1,:)).*DCfreq.*(nthroot(49*S1,hs))), t);  % 
DC2tc = @(t) interp1(TSPAN, (Fig6adResults(1,:).*DCfreq.*(nthroot(49*S1,hs))), t); % 
for i = 1:numSims
    F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(Y(3)+DC1tc(t))^hs)/(S1^hs+(Y(3)+DC1tc(t))^hs))*((Z2^hz)/(Z2^hz+(Y(4)+DC2tc(t))^hz)) - dTF1*Y(1);
        (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(Y(4)+DC2tc(t))^hs)/(S2^hs+(Y(4)+DC2tc(t))^hs))*((Z1^hz)/(Z1^hz+(Y(3)+DC1tc(t))^hz)) - dTF2*Y(2);
        ((a1*Y(1)^ha)/(A1^ha+Y(1)^ha))*((R2^hr)/(R2^hr+Y(2)^hr))*((U2^hu)/(U2^hu+(Y(4)+DC2tc(t))^hu)) - dCY1*Y(3);
        ((a2*Y(2)^ha)/(A2^ha+Y(2)^ha))*((R1^hr)/(R1^hr+Y(1)^hr))*((U1^hu)/(U1^hu+(Y(3)+DC1tc(t))^hu)) - dCY2*Y(4)];
    G = @(t,Y) [nTF1*Y(1);
        nTF2*Y(2);
        nCY1*Y(3);
        nCY2*Y(4)];
    TSPAN = 0:0.1:168;                           
    Y0 = [0 0 0 0];                        
    options = sdeset('SDEType','Ito');          
    Y = sde_euler(F,G,TSPAN,Y0,options);        
    CY1Fracs = ((Y(:,1).^ha)./(A1^ha+Y(:,1).^ha)).*((R2^hr)./(R2^hr+Y(:,2).^hr)).*((U2^hu)./(U2^hu+Y(:,4).^hu));
    CY2Fracs = ((Y(:,2).^ha)./(A2^ha+Y(:,2).^ha)).*((R1^hr)./(R1^hr+Y(:,1).^hr)).*((U1^hu)./(U1^hu+Y(:,3).^hu));
    CYBalance = (atan(CY2Fracs./CY1Fracs))./(pi/2);
    CYBalance(1) = 0.5; 
    Fig6dResults(i,:) = CYBalance;
    fprintf('Finished Simulation %d \n', i);
end
%save('Fig6dResults_FOR_MS2','Fig6dResults');
%load('Fig6dResults_FOR_MS2');

%% FIGURE 6A-D (Plot)
tiledlayout(4,1)

nexttile
hold on
crv4 = plot(TSPAN, Fig6adResults(1,:), 'Color', [0 0.75 0 0.1],'LineWidth',1);
crv5 = plot(TSPAN, Fig6adResults(1,:), 'Color', [0 0.75 0], 'LineWidth', 3);
crv1 = plot(TSPAN, Fig6adResults(1,:), 'Color',[0.65 0.65 0.65], 'LineWidth', 12);
crv2 = plot(TSPAN, 1-Fig6adResults(6,:), 'Color','k', 'LineWidth', 6);
crv3 = plot(TSPAN, Fig6aResults, 'Color', [0 0.5 0],'LineWidth', 3);
%title('Deterministic Th Response to Sample APC Instruction');
%xlabel('Time (h)', 'FontSize', 14)
xlim([0 168])
%ylabel('Effector Bias', 'FontSize',14)
yticks([0 0.5 1])
yticklabels({'Th1', 'Mixed', 'Th2'})
legend([crv1 crv2 crv3 crv5 crv4],{'APC Instruction', '"Optimal" Th Effector Balance', 'Mean Realized Effector Balance', 'Median Realized Effector Balance' 'Sample Realized Effector Balances'},'FontSize',8)
text(-8, 1,'a','FontSize',14);
hold off

nexttile
hold on
for i = 1:numSims
    crv4 = plot(TSPAN, Fig6bResults(i,:),'Color',[0 0.75 0 0.1],'Linewidth',2);
end
crv1 = plot(TSPAN, Fig6adResults(1,:),'Color',[0.6 0.6 0.6],'LineWidth',12);
crv2 = plot(TSPAN, 1-Fig6adResults(7,:),'Color','k','LineWidth',6);
crv3 = plot(TSPAN, mean(Fig6bResults,1), 'Color',[0 0.5 0],'LineWidth',3);
crv5 = plot(TSPAN, median(Fig6bResults,1), 'Color',[0 0.75 0],'LineWidth',3);
%title('Stochastic Th Response to Sample APC Instruction');
%xlabel('Time (h)', 'FontSize', 14)
xlim([0 168])
ylabel('Effector Bias                               ', 'FontSize',14)
ylim([0 1])
yticks([0 0.5 1])
yticklabels({'Th1', 'Mixed', 'Th2'})
%legend([crv1 crv2 crv4 crv3],{'APC Instruction', '"Optimal" Th Cytokine Balance', 'Simulated Realized Th Cytokine Balances', 'Mean Realized Th Cytokine Balance'},'FontSize',10,'Location','southwest')
text(-8, 1,'b','FontSize',14);
hold off

nexttile
hold on
for i = 1:numSims
    crv4 = plot(TSPAN, Fig6cResults(i,:),'Color',[0 0.75 0 0.1],'Linewidth',2);
end
crv1 = plot(TSPAN, Fig6adResults(1,:),'Color',[0.6 0.6 0.6],'LineWidth',12);
crv2 = plot(TSPAN, 1-Fig6adResults(7,:),'Color','k','LineWidth',6);
crv3 = plot(TSPAN, mean(Fig6cResults,1), 'Color',[0 0.5 0],'LineWidth',3);
crv5 = plot(TSPAN, median(Fig6cResults,1), 'Color',[0 0.75 0],'LineWidth',3);
%title('Stochastic Th Response to Sample APC Instruction');
%xlabel('Time (h)', 'FontSize', 14)
xlim([0 168])
%ylabel('Effector Bias', 'FontSize',14)
ylim([0 1])
yticks([0 0.5 1])
yticklabels({'Th1', 'Mixed', 'Th2'})
%legend([crv1 crv2 crv4 crv3],{'APC Instruction', '"Optimal" Th Cytokine Balance', 'Simulated Realized Th Cytokine Balances', 'Mean Realized Th Cytokine Balance'},'FontSize',10,'Location','southwest')
text(-8, 1,'c','FontSize',14);
hold off

nexttile
hold on
for i = 1:numSims
    crv4 = plot(TSPAN, Fig6dResults(i,:),'Color',[0 0.75 0 0.1],'Linewidth',2);
end
crv1 = plot(TSPAN, Fig6adResults(1,:),'Color',[0.6 0.6 0.6],'LineWidth',12);
crv2 = plot(TSPAN, 1-Fig6adResults(7,:),'Color','k','LineWidth',6);
Fig6dResults(Fig6dResults<0) = NaN;
crv3 = plot(TSPAN, mean(Fig6dResults,1, 'omitnan'), 'Color',[0 0.5 0],'LineWidth',3);
crv5 = plot(TSPAN, median(Fig6dResults,1, 'omitnan'), 'Color',[0 0.75 0],'LineWidth',3);
%title('Stochastic Th Response to Sample APC Instruction');
xlabel('Time (h)', 'FontSize', 14)
xlim([0 168])
%ylabel('Effector Bias', 'FontSize',14)
ylim([0 1])
yticks([0 0.5 1])
yticklabels({'Th1', 'Mixed', 'Th2'})
%legend([crv1 crv2 crv4 crv3],{'APC Instruction', '"Optimal" Th Cytokine Balance', 'Simulated Realized Th Cytokine Balances', 'Mean Realized Th Cytokine Balance'},'FontSize',10,'Location','southwest')
text(-8, 1,'d','FontSize',14);
hold off

annotation('textarrow',[0.325 0.275], [0.85, 0.875],'String',['Sabotaged APCs,' newline 'Best to Ignore']);
annotation('textarrow',[0.555 0.505], [0.84, 0.8],'String',['Trustworthy APCs,' newline 'Best to Obey']);

%% FIGURE 6E-F (Data): For 100 APC time courses, measure SDE performance across noise levels (across sensitivity values).
numSims = 100;
DCfreq = 0.17; 
threshold = [20; 40; 60; 80; 100; 120; 140; 160; 180; 200; 220; 240]; 
timecourselength = 7*24*10+1; % Length of time course in time steps.
TSPAN = 0:0.1:168;
nTFvals = linspace(0,1,11);
nCYvals = linspace(0,1,11);
Fig6efResults = nan(length(nTFvals), length(nCYvals), length(threshold), numSims);
for j = 1:length(nTFvals)
    for k = 1:length(nCYvals)
        for i = 1:numSims
            % Make a DC bias time course, with "optimum" Th responses for all sensitivity ("threshold") values.
            DC1bias = nan(1,timecourselength); % Initialize the DC1 bias time course.
            DC1bias(1) = binornd(1,0.5); % First time step is either 0 or 1.
            count = 1; % This counts which time step we're on.
            streak = poissrnd(120); % How many time steps the initial DC bias will last.
            DC1bias(count+1:count+streak) = DC1bias(1,1); % Set those first X time steps.
            Th1response = repmat(DC1bias, [length(threshold),1]); % Initialize the "optimal" Th1 cell response. So far, it is the same as the DC1 bias time course.
            flag = (ones(length(threshold),1)).*(count+streak);
            curropt = Th1response(:,1); % This keeps track of the current optimal Th1 cell response, for each threshold value.
            integral = zeros(length(threshold),1); % This counts the total integral that might cause a switch in the optimal response, for each threshold value.
            count = count + streak; % Update the most recent time step of the DC1 bias time course to have been filled in.
            while (count < timecourselength)
                new = round((betarnd(DC1bias(count)+0.1, 1.1-DC1bias(count))), 1);
                streak = min(poissrnd(20+(400*((new-0.5)^2))),timecourselength-count);
                DC1bias(count+1:count+streak) = new;
                Chal = (abs(new-curropt) > 0.5);
                nChalUpdate = repmat(curropt,[1,timecourselength]);
                rowixn = find(~Chal); % For-loop matrix update code from Sean de Wolski at Matlab.
                for ii = 1:numel(rowixn) % For-loop matrix update code from Sean de Wolski at Matlab.
                    rowiin = rowixn(ii); % For-loop matrix update code from Sean de Wolski at Matlab.
                    colixn = flag(rowiin):(count+streak); % For-loop matrix update code from Sean de Wolski at Matlab.
                    Th1response(rowiin, colixn) = nChalUpdate(rowiin, colixn); % For-loop matrix update code from Sean de Wolski at Matlab.
                end % For-loop matrix update code from Sean de Wolski at Matlab.
                flag(~Chal) = count + streak;
                integral(~Chal) = 0;
                yChalInt = (abs(new-curropt)).*streak;
                integral(Chal) = integral(Chal) + yChalInt(Chal);
                Cros = (Chal & (integral>=threshold));
                curropt(Chal & Cros) = 1 - curropt(Chal & Cros);
                yChalUpdate = repmat(curropt,[1,timecourselength]); % Start edits after this line
                rowixy = find(Chal & Cros); % For-loop matrix update code from Sean de Wolski at Matlab.
                for ii = 1:numel(rowixy) % For-loop matrix update code from Sean de Wolski at Matlab.
                    rowiiy = rowixy(ii); % For-loop matrix update code from Sean de Wolski at Matlab.
                    colixy = flag(rowiiy):(count+streak); % For-loop matrix update code from Sean de Wolski at Matlab.
                    Th1response(rowiiy, colixy) = yChalUpdate(rowiiy, colixy); % For-loop matrix update code from Sean de Wolski at Matlab.
                end % For-loop matrix update code from Sean de Wolski at Matlab.
                integral(Chal & Cros) = 0;
                flag(Chal & Cros) = count + streak;
                if (count+streak >= timecourselength)
                    Th1response((Chal & ~Cros), flag(Chal & ~Cros):(count+streak)) = yChalUpdate((Chal & ~Cros), flag(Chal & ~Cros):(count+streak));
                end
                count = count + streak;
            end
            DC2bias = 1-DC1bias;
            DC1tc = @(t) interp1(TSPAN, (DC1bias.*DCfreq.*(nthroot(49*S1,hs))), t);
            DC2tc = @(t) interp1(TSPAN, (DC2bias.*DCfreq.*(nthroot(49*S1,hs))), t);
            % Run an SDE sample path.
            F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(Y(3)+DC1tc(t))^hs)/(S1^hs+(Y(3)+DC1tc(t))^hs))*((Z2^hz)/(Z2^hz+(Y(4)+DC2tc(t))^hz)) - dTF1*Y(1);
                (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(Y(4)+DC2tc(t))^hs)/(S2^hs+(Y(4)+DC2tc(t))^hs))*((Z1^hz)/(Z1^hz+(Y(3)+DC1tc(t))^hz)) - dTF2*Y(2);
                ((a1*Y(1)^ha)/(A1^ha+Y(1)^ha))*((R2^hr)/(R2^hr+Y(2)^hr))*((U2^hu)/(U2^hu+(Y(4)+DC2tc(t))^hu)) - dCY1*Y(3);
                ((a2*Y(2)^ha)/(A2^ha+Y(2)^ha))*((R1^hr)/(R1^hr+Y(1)^hr))*((U1^hu)/(U1^hu+(Y(3)+DC1tc(t))^hu)) - dCY2*Y(4)];
            G = @(t,Y) [nTFvals(j)*Y(1);
                nTFvals(j)*Y(2);
                nCYvals(k)*Y(3);
                nCYvals(k)*Y(4)];
            TSPAN = 0:0.1:168;                           
            Y0 = [0 0 0 0];                        
            options = sdeset('SDEType','Ito');          
            Y = sde_euler(F,G,TSPAN,Y0,options);        
            CY1Fracs = ((Y(:,1).^ha)./(A1^ha+Y(:,1).^ha)).*((R2^hr)./(R2^hr+Y(:,2).^hr)).*((U2^hu)./(U2^hu+Y(:,4).^hu));
            CY2Fracs = ((Y(:,2).^ha)./(A2^ha+Y(:,2).^ha)).*((R1^hr)./(R1^hr+Y(:,1).^hr)).*((U1^hu)./(U1^hu+Y(:,3).^hu));
            CYBalance = (atan(CY1Fracs./CY2Fracs))./(pi/2);
            CYBalance(1) = 0.5;
            % Measure MI or well-matchedness for every "optimal" Th response sensitivity.
            for p = 1:length(threshold)
                if ~isempty(TCDM_well_matchedness((Th1response(p,:))', CYBalance, timecourselength)) % Need some way to deal with when atan becomes negative or MI can't be calculated.
                    Fig6efResults(j,k,p,i) = TCDM_well_matchedness((Th1response(p,:))', CYBalance, timecourselength);
                end
            end
            fprintf('Finished a Simulation: %d \n', i);
        end
    fprintf('Finished a Noise Combo: %d \n', ((j-1)*10)+k);
    end
end
%save('Fig6efResults_FOR_MS2', 'Fig6efResults');
%load('Fig6efResults_FOR_MS2');

%% FIGURE 6E (Plot)
linePlotMeanHiSens = nan(1,length(Fig6efResults(:,1,1,1)));
linePlotMeanMeSens = nan(1,length(Fig6efResults(:,1,1,1)));
linePlotMeanLoSens = nan(1,length(Fig6efResults(:,1,1,1)));
linePlotStdDvHiSens = nan(1,length(Fig6efResults(:,1,1,1)));
linePlotStdDvMeSens = nan(1,length(Fig6efResults(:,1,1,1)));
linePlotStdDvLoSens = nan(1,length(Fig6efResults(:,1,1,1)));
for i = 1:length(linePlotMeanHiSens)
    for j = 1:3
        if j==1
            linePlotMeanLoSens(i) = nanmean(Fig6efResults(i,i,12,:));
            linePlotStdDvLoSens(i) = nanstd(Fig6efResults(i,i,12,:))./(sum(~isnan(Fig6efResults(i,i,12,:))));
        elseif j == 2
            linePlotMeanMeSens(i) = nanmean(Fig6efResults(i,i,6,:));
            linePlotStdDvMeSens(i) = nanstd(Fig6efResults(i,i,6,:))./(sum(~isnan(Fig6efResults(i,i,6,:))));
        else
            linePlotMeanHiSens(i) = nanmean(Fig6efResults(i,i,1,:));
            linePlotStdDvHiSens(i) = nanstd(Fig6efResults(i,i,1,:))./(sum(~isnan(Fig6efResults(i,i,1,:))));
        end
    end
end
%Subtract minimum value from every measurement.
[M,I] = min(linePlotMeanHiSens);
normConst = M - linePlotStdDvHiSens(I);
linePlotMeanHiSens = linePlotMeanHiSens - normConst;
linePlotMeanMeSens = linePlotMeanMeSens - normConst;
linePlotMeanLoSens = linePlotMeanLoSens - normConst;
%Divide every measurement by the maximum value.
[M,I] = max(linePlotMeanLoSens);
normConst = M + linePlotStdDvLoSens(I);
linePlotMeanHiSens = linePlotMeanHiSens./normConst;
linePlotMeanMeSens = linePlotMeanMeSens./normConst;
linePlotMeanLoSens = linePlotMeanLoSens./normConst;
linePlotStdDvHiSens = linePlotStdDvHiSens./normConst;
linePlotStdDvMeSens = linePlotStdDvMeSens./normConst;
linePlotStdDvLoSens = linePlotStdDvLoSens./normConst;
%Plot.
hold on
crv1 = errorbar(linspace(0,1,length(Fig6efResults(:,1,1,1))).*100, linePlotMeanLoSens, linePlotStdDvLoSens, 'Color', [0 0 0], 'LineWidth', 3);
crv2 = errorbar(linspace(0,1,length(Fig6efResults(:,1,1,1))).*100, linePlotMeanMeSens, linePlotStdDvMeSens, 'Color', [0.3 0.3 0.3], 'LineWidth', 3);
crv3 = errorbar(linspace(0,1,length(Fig6efResults(:,1,1,1))).*100, linePlotMeanHiSens, linePlotStdDvHiSens, 'Color', [0.6 0.6 0.6], 'LineWidth', 3);
title('Cell-Cell Variability Helps Track APC Instruction','FontSize',18)
xlabel('Percentage Volatility in Molecular Expression', 'FontSize', 14)
ylabel('Relative Performance', 'FontSize',14)
ylim([0 1])
legend({'Low Sensitivity', 'Medium Sensitivity', 'High Sensitivity'}, 'Location','southeast')
text(-10,1,'e','FontSize',24);
hold off

%% FIGURE 6F (Plot)
hmData = nan(length(Fig6efResults(:,1,1,1)),length(Fig6efResults(:,1,1,1)));
for i = 1:length(Fig6efResults(:,1,1,1))
    for j = 1:length(Fig6efResults(1,:,1,1))
        hmData(i,j) = nanmean(Fig6efResults(i,j,6,:));
    end
end
hmData = hmData - (min(min(hmData)));
hmData = hmData./(max(max(hmData)));
hmData = flipud(hmData);
xlabels = {'0','10','20','30','40','50','60','70','80','90','100'};
ylabels = {'100','90','80','70','60','50','40','30','20','10','0'};
h = heatmap(xlabels, ylabels, hmData, 'Colormap', jet, 'ColorbarVisible','off');
h.XLabel = 'Cytokine Percentage Volatility';
h.YLabel = 'Transcription Factor Percentage Volatility';
annotation('textbox', [0.025 0.5 0.5 0.5], 'String', 'f', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);

%% FIGURE 6F (Colorbar)
axis off
colormap(jet);
cb = colorbar('Ticks',[0 1], 'TickLabels',{'Worst','Best'});
ylabel(cb,'Relative Performance','FontSize',16)

%% SUPP FIGURE 1: Effect of Cell Density on Cytokine Production and Removal
cellDensities = logspace(6, 9, 1001);
cvpcs = (1-(cellDensities.*(((diam/2)^3).*(4/3).*pi)./(10^12)))./(cellDensities.*(((diam/2)^3).*(4/3).*pi)./(10^12));
aVals = ((20*60*60)./cvpcs); %./(((1.66054*10^3)*(((diam/2)^3)*(4/3)*pi))/((15.3*10^3)*(1-pack_eff)));
dVals = (0.8./cvpcs + 0.015); %./(((1.66054*10^3)*(((diam/2)^3)*(4/3)*pi))/((15.3*10^3)*(1-pack_eff)));
crv1 = loglog(cellDensities, aVals, 'LineWidth', 6, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
hold on
crv2 = loglog(cellDensities, dVals, 'LineWidth', 6, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
crv3 = loglog(cellDensities, aVals./dVals, 'LineWidth', 8, 'Color', 'k');
% title('Cytokine Production and Removal Scale Differently with Cell Density');
ylabel('Rate of Change of [Cytokine]', 'FontSize', 14);
xlabel('Cell Density (cells/mL)', 'FontSize', 14);
legend([crv1 crv2 crv3],{'Cytokine Production', 'Cytokine Removal', 'Production:Removal Ratio'}, 'FontSize',10)
hold off

%% SUPP FIGURE 2a-b (Data): Convergence to Mixed Effector Type Takes a Long Time In Vitro
% Set the appropriate parameters: cpmL = 2*10^6, R and Z asymmetric!
Cy1In = [0 0 2.5 15 90 540];
Cy2In = [0 0 2.5 15 90 540];
Tf1Out = zeros(length(Cy1In));
Tf2Out = zeros(length(Cy2In));
Cy1Out = zeros(length(Cy1In));
Cy2Out = zeros(length(Cy2In));
numSims = 1;
a1old = a1;
a2old = a2;
for i = 1:length(Cy1In)
    for j = 1:length(Cy2In)
        for k = 1:numSims
            if i==1
                a1 = 0.05*a1old; % Antibody against IFNg eliminates 95% of IFNg from the media.
            end
            if j==1
                a2 = 0.05*a2old; % Antibody against IL-4 eliminates 95% of IL-4 from the media.
            end
            F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(DC1+Y(3))^hs)/(S1^hs+(DC1+Y(3))^hs))*(Z2^hz)/(Z2^hz+Y(4)^hz) - dTF1*Y(1);
                (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(DC2+Y(4))^hs)/(S2^hs+(DC2+Y(4))^hs))*(Z1^hz)/(Z1^hz+Y(3)^hz) - dTF2*Y(2);
                (a1*Y(1)^ha)/(A1^ha+Y(1)^ha)*(R2^hr)/(R2^hr+Y(2)^hr)*(U2^hu)/(U2^hu+Y(4)^hu) - dCY1*Y(3);
                (a2*Y(2)^ha)/(A2^ha+Y(2)^ha)*(R1^hr)/(R1^hr+Y(1)^hr)*(U1^hu)/(U1^hu+Y(3)^hu) - dCY2*Y(4)];
            G = @(t,Y) [nTF1*Y(1);
                nTF2*Y(2);
                nCY1*Y(3);
                nCY2*Y(4)];
            TSPAN = 0:0.1:96;
            EX1 = Cy1In(i)*((1.66054*10^3)*(((diam/2)^3)*(4/3)*pi))/((17.1*10^3)*(1-pack_eff));  % Calculation for true concentration of IL-12
            EX2 = Cy1In(j)*((1.66054*10^3)*(((diam/2)^3)*(4/3)*pi))/((13.5*10^3)*(1-pack_eff));% Calculation for true concentration of IL-4
            Y0 = [0 0 EX1 EX2];
            options = sdeset('SDEType','Ito');
            Y = sde_euler(F,G,TSPAN,Y0,options);
            TSPANa = 0:0.1:72;
            bold = b;
            b = 0;
            Y0a = [Y(length(TSPAN),1) Y(length(TSPAN),2) 0 0];
            F = @(t,Ya) [(b + (p1*Ya(1)^hp)/(P1^hp+Ya(1)^hp))*(X2^hx)/(X2^hx+Ya(2)^hx) + ((s1*(DC1+Ya(3))^hs)/(S1^hs+(DC1+Ya(3))^hs))*(Z2^hz)/(Z2^hz+Ya(4)^hz) - dTF1*Ya(1);
                (b + (p2*Ya(2)^hp)/(P2^hp+Ya(2)^hp))*(X1^hx)/(X1^hx+Ya(1)^hx) + ((s2*(DC2+Ya(4))^hs)/(S2^hs+(DC2+Ya(4))^hs))*(Z1^hz)/(Z1^hz+Ya(3)^hz) - dTF2*Ya(2);
                (a1*Ya(1)^ha)/(A1^ha+Ya(1)^ha)*(R2^hr)/(R2^hr+Ya(2)^hr)*(U2^hu)/(U2^hu+Ya(4)^hu) - dCY1*Ya(3);
                (a2*Ya(2)^ha)/(A2^ha+Ya(2)^ha)*(R1^hr)/(R1^hr+Ya(1)^hr)*(U1^hu)/(U1^hu+Ya(3)^hu) - dCY2*Ya(4)];
            G = @(t,Y) [nTF1*Y(1);
                nTF2*Y(2);
                nCY1*Y(3);
                nCY2*Y(4)];
            Ya = sde_euler(F,G,TSPANa,Y0a,options);
            TSPANb = 0:0.1:6;
            b = bold;
            Y0b = [Ya(length(TSPANa),1) Ya(length(TSPANa),2) Ya(length(TSPANa),3) Ya(length(TSPANa),4)];
            F = @(t,Yb) [(b + (p1*Yb(1)^hp)/(P1^hp+Yb(1)^hp))*(X2^hx)/(X2^hx+Yb(2)^hx) + ((s1*(DC1+Yb(3))^hs)/(S1^hs+(DC1+Yb(3))^hs))*(Z2^hz)/(Z2^hz+Yb(4)^hz) - dTF1*Yb(1);
                (b + (p2*Yb(2)^hp)/(P2^hp+Yb(2)^hp))*(X1^hx)/(X1^hx+Yb(1)^hx) + ((s2*(DC2+Yb(4))^hs)/(S2^hs+(DC2+Yb(4))^hs))*(Z1^hz)/(Z1^hz+Yb(3)^hz) - dTF2*Yb(2);
                (a1*Yb(1)^ha)/(A1^ha+Yb(1)^ha)*(R2^hr)/(R2^hr+Yb(2)^hr)*(U2^hu)/(U2^hu+Yb(4)^hu) - dCY1*Yb(3);
                (a2*Yb(2)^ha)/(A2^ha+Yb(2)^ha)*(R1^hr)/(R1^hr+Yb(1)^hr)*(U1^hu)/(U1^hu+Yb(3)^hu) - dCY2*Yb(4)];
            G = @(t,Y) [nTF1*Y(1);
                nTF2*Y(2);
                nCY1*Y(3);
                nCY2*Y(4)];
            Yb = sde_euler(F,G,TSPANb,Y0b,options);
            Tf1Out(i,j) = Tf1Out(i,j) + Yb(length(TSPANb),1);
            Tf2Out(i,j) = Tf2Out(i,j) + Yb(length(TSPANb),2);
            Cy1Out(i,j) = Cy1Out(i,j) + Yb(length(TSPANb),3);
            Cy2Out(i,j) = Cy2Out(i,j) + Yb(length(TSPANb),4);
            if i==1
                a1 = a1old;
            end
            if j==1
                a2 = a2old;
            end
        end
        fprintf('Percent Completed: %d\n', round(100*((i-1)/length(Cy1In)+j/(length(Cy1In)^2))));
    end
end
TF1Out = (Tf1Out - min(Tf1Out(:)))./(max(Tf1Out(:)) - min(Tf1Out(:)));
TF2Out = (Tf2Out - min(Tf2Out(:)))./(max(Tf2Out(:)) - min(Tf2Out(:)));
CY1Out = (Cy1Out - min(Cy1Out(:)))./(max(Cy1Out(:)) - min(Cy1Out(:)));
CY2Out = (Cy2Out - min(Cy2Out(:)))./(max(Cy2Out(:)) - min(Cy2Out(:)));

%% SUPP FIGURE 2a (Plot)
xlabels = {'Ab','0','2.5','15','90','540'};
ylabels = {'540','90','15','2.5','0','Ab'};
h = heatmap(xlabels, ylabels, flipud(TF1Out), 'Colormap', jet, 'ColorbarVisible','off', 'CellLabelColor','none');
h.Title = 'Mean Cellular T-bet Expression';
h.XLabel = 'External IL-4 (ng/mL)';
h.YLabel = 'External IL-12 (ng/mL)';
annotation('textbox', [0.025 0.5 0.5 0.5], 'String', 'a', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
annotation('textbox', [0.262 0.658 0.125 0.13], 'String', 'c', 'FitBoxToText', 'off','EdgeColor','k', 'LineWidth', 6, 'FontSize', 30, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', 'w');
annotation('textbox', [0.390 0.525 0.125 0.13], 'String', 'd', 'FitBoxToText', 'off','EdgeColor','k', 'LineWidth', 6, 'FontSize', 30, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', 'k');
annotation('textbox', [0.52 0.388 0.125 0.13], 'String', 'e', 'FitBoxToText', 'off','EdgeColor','k', 'LineWidth', 6, 'FontSize', 30, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', 'k');
annotation('textbox', [0.65 0.255 0.125 0.13], 'String', 'f', 'FitBoxToText', 'off','EdgeColor','k', 'LineWidth', 6, 'FontSize', 30, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', 'w');

%% SUPP FIGURE 2b (Plot)
xlabels = {'Ab','0','2.5','15','90','540'};
ylabels = {'540','90','15','2.5','0','Ab'};
h = heatmap(xlabels, ylabels, flipud(TF2Out), 'Colormap', jet, 'ColorbarVisible','off', 'CellLabelColor','none');
h.Title = 'Mean Cellular GATA3 Expression';
h.XLabel = 'External IL-4 (ng/mL)';
h.YLabel = 'External IL-12 (ng/mL)';
annotation('textbox', [0.025 0.5 0.5 0.5], 'String', 'b', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
annotation('textbox', [0.262 0.658 0.125 0.13], 'String', 'c', 'FitBoxToText', 'off','EdgeColor','k', 'LineWidth', 6, 'FontSize', 30, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', 'w');
annotation('textbox', [0.390 0.525 0.125 0.13], 'String', 'd', 'FitBoxToText', 'off','EdgeColor','k', 'LineWidth', 6, 'FontSize', 30, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', 'k');
annotation('textbox', [0.52 0.388 0.125 0.13], 'String', 'e', 'FitBoxToText', 'off','EdgeColor','k', 'LineWidth', 6, 'FontSize', 30, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', 'k');
annotation('textbox', [0.65 0.255 0.125 0.13], 'String', 'f', 'FitBoxToText', 'off','EdgeColor','k', 'LineWidth', 6, 'FontSize', 30, 'FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', 'w');

%% SUPP FIGURE 2a-b (Colorbar)
axis off
colormap(jet);
cb = colorbar('Ticks',[0 1], 'TickLabels',{'Min','Max'});
ylabel(cb,'Relative Expression','FontSize',16)

%% SUPP FIGURE 2c-f (Data): Single Runs of Antebi's Experiments Across 4 Starting conditions.
% Switch back to symmetric paramters here, for illustration purposes.
% But in reality, better predictions come from using asymmetric parameters.
Ymaster = NaN(16800,4,4);
Cy1In = [90 15 2.5 0];
Cy2In = [0 2.5 15 90];
for i = 1:4 % Which of the 4 cytokine pairs you'll use: 1 = 90,0   2 = 15,2.5   3 = 2.5,15   4 = 0,90.
    F = @(t,Y) [(b + (p1*Y(1)^hp)/(P1^hp+Y(1)^hp))*(X2^hx)/(X2^hx+Y(2)^hx) + ((s1*(DC1+Y(3))^hs)/(S1^hs+(DC1+Y(3))^hs))*(Z2^hz)/(Z2^hz+Y(4)^hz) - dTF1*Y(1);
        (b + (p2*Y(2)^hp)/(P2^hp+Y(2)^hp))*(X1^hx)/(X1^hx+Y(1)^hx) + ((s2*(DC2+Y(4))^hs)/(S2^hs+(DC2+Y(4))^hs))*(Z1^hz)/(Z1^hz+Y(3)^hz) - dTF2*Y(2);
        (a1*Y(1)^ha)/(A1^ha+Y(1)^ha)*(R2^hr)/(R2^hr+Y(2)^hr)*(U2^hu)/(U2^hu+Y(4)^hu) - dCY1*Y(3);
        (a2*Y(2)^ha)/(A2^ha+Y(2)^ha)*(R1^hr)/(R1^hr+Y(1)^hr)*(U1^hu)/(U1^hu+Y(3)^hu) - dCY2*Y(4)];
    G = @(t,Y) [nTF1*Y(1);
        nTF2*Y(2);
        nCY1*Y(3);
        nCY2*Y(4)];
    TSPAN = 0:0.1:96;
    EX1 = Cy1In(i)*((1.66054*10^3)*(((diam/2)^3)*(4/3)*pi))/((17.1*10^3)*(1-pack_eff));% Calculation for true concentration of IL-12
    EX2 = Cy2In(i)*((1.66054*10^3)*(((diam/2)^3)*(4/3)*pi))/((13.5*10^3)*(1-pack_eff));% Calculation for true concentration of IL-4
    Y0 = [0 0 EX1 EX2];
    options = sdeset('SDEType','Ito');
    Y = sde_euler(F,G,TSPAN,Y0,options);
    bold = b;
    b = 0;
    TSPANa = 0:0.1:72;
    Y0a = [Y(length(TSPAN),1) Y(length(TSPAN),2) 0 0];
    F = @(t,Ya) [(b + (p1*Ya(1)^hp)/(P1^hp+Ya(1)^hp))*(X2^hx)/(X2^hx+Ya(2)^hx) + ((s1*(DC1+Ya(3))^hs)/(S1^hs+(DC1+Ya(3))^hs))*(Z2^hz)/(Z2^hz+Ya(4)^hz) - dTF1*Ya(1);
        (b + (p2*Ya(2)^hp)/(P2^hp+Ya(2)^hp))*(X1^hx)/(X1^hx+Ya(1)^hx) + ((s2*(DC2+Ya(4))^hs)/(S2^hs+(DC2+Ya(4))^hs))*(Z1^hz)/(Z1^hz+Ya(3)^hz) - dTF2*Ya(2);
        (a1*Ya(1)^ha)/(A1^ha+Ya(1)^ha)*(R2^hr)/(R2^hr+Ya(2)^hr)*(U2^hu)/(U2^hu+Ya(4)^hu) - dCY1*Ya(3);
        (a2*Ya(2)^ha)/(A2^ha+Ya(2)^ha)*(R1^hr)/(R1^hr+Ya(1)^hr)*(U1^hu)/(U1^hu+Ya(3)^hu) - dCY2*Ya(4)];
    G = @(t,Y) [nTF1*Y(1);
        nTF2*Y(2);
        nCY1*Y(3);
        nCY2*Y(4)];
    Ya = sde_euler(F,G,TSPANa,Y0a,options);
    b = bold;
    TSPANb = 0:0.1:1511.7;
    Y0b = [Ya(length(TSPANa),1) Ya(length(TSPANa),2) Ya(length(TSPANa),3) Ya(length(TSPANa),4)];
    F = @(t,Yb) [(b + (p1*Yb(1)^hp)/(P1^hp+Yb(1)^hp))*(X2^hx)/(X2^hx+Yb(2)^hx) + ((s1*(DC1+Yb(3))^hs)/(S1^hs+(DC1+Yb(3))^hs))*(Z2^hz)/(Z2^hz+Yb(4)^hz) - dTF1*Yb(1);
        (b + (p2*Yb(2)^hp)/(P2^hp+Yb(2)^hp))*(X1^hx)/(X1^hx+Yb(1)^hx) + ((s2*(DC2+Yb(4))^hs)/(S2^hs+(DC2+Yb(4))^hs))*(Z1^hz)/(Z1^hz+Yb(3)^hz) - dTF2*Yb(2);
        (a1*Yb(1)^ha)/(A1^ha+Yb(1)^ha)*(R2^hr)/(R2^hr+Yb(2)^hr)*(U2^hu)/(U2^hu+Yb(4)^hu) - dCY1*Yb(3);
        (a2*Yb(2)^ha)/(A2^ha+Yb(2)^ha)*(R1^hr)/(R1^hr+Yb(1)^hr)*(U1^hu)/(U1^hu+Yb(3)^hu) - dCY2*Yb(4)];
    G = @(t,Y) [nTF1*Y(1);
        nTF2*Y(2);
        nCY1*Y(3);
        nCY2*Y(4)];
    Yb = sde_euler(F,G,TSPANb,Y0b,options);
    TSPANmaster = 0:0.1:(length(TSPAN) + length(TSPANa) + length(TSPANb) - 1)/10;
    Ymaster(:,:,i) = vertcat(Y,Ya,Yb);
end

%% SUPP FIGURE 2c-f (Plot)
tiledlayout(4,1);

nexttile
hold on
crv1 = plot(TSPANmaster,Ymaster(:,1,1),'Color',color1,'LineWidth',5) 
crv2 = plot(TSPANmaster,Ymaster(:,2,1),'Color',color2,'LineWidth',5) 
ln1 = line([170 170], [0 200], 'Color', 'k', 'LineWidth', 3, 'LineStyle',':')
title('Initial Stimulation: [IFNg] = 90ng/mL; [IL-4] = 0ng/mL', 'FontSize',12)
xlim([0 1680])
xticks([0 96 168 336 504 672 840 1008 1176 1344 1512 1680])
xticklabels({'Stim+','Stim-','1','2','3','4','5','6','7','8','9','10'})
ylim([0 200])
yticks([0 50 100 150 200])
legend([crv1 crv2],{'T-bet', 'GATA3'},'FontSize',8)
annotation('textbox', [0.12 0.50 0.5 0.5], 'String', 'c', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
hold off

nexttile
hold on
crv1 = plot(TSPANmaster,Ymaster(:,1,2),'Color',color1,'LineWidth',5) 
crv2 = plot(TSPANmaster,Ymaster(:,2,2),'Color',color2,'LineWidth',5) 
ln1 = line([170 170], [0 200], 'Color', 'k', 'LineWidth', 3, 'LineStyle',':')
title('Initial Stimulation: [IFNg] = 15ng/mL; [IL-4] = 2.5ng/mL', 'FontSize',12)
xlim([0 1680])
xticks([0 96 168 336 504 672 840 1008 1176 1344 1512 1680])
xticklabels({'Stim+','Stim-','1','2','3','4','5','6','7','8','9','10'})
ylim([0 200])
yticks([0 50 100 150 200])
ylabel('Molecules Per Th Cell                            ','FontSize',12)
annotation('textbox', [0.12 0.272 0.5 0.5], 'String', 'd', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
hold off

nexttile
hold on
crv1 = plot(TSPANmaster,Ymaster(:,1,3),'Color',color1,'LineWidth',5) % Plot the numerical solution.
crv2 = plot(TSPANmaster,Ymaster(:,2,3),'Color',color2,'LineWidth',5) 
ln1 = line([170 170], [0 200], 'Color', 'k', 'LineWidth', 3, 'LineStyle',':')
title('Initial Stimulation: [IFNg] = 2.5ng/mL; [IL-4] = 15ng/mL', 'FontSize',12)
xlim([0 1680])
xticks([0 96 168 336 504 672 840 1008 1176 1344 1512 1680])
xticklabels({'Stim+','Stim-','1','2','3','4','5','6','7','8','9','10'})
ylim([0 200])
yticks([0 50 100 150 200])
annotation('textbox', [0.12 0.045 0.5 0.5], 'String', 'e', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
hold off

nexttile
hold on
crv1 = plot(TSPANmaster,Ymaster(:,1,4),'Color',color1,'LineWidth',5) % Plot the numerical solution.
crv2 = plot(TSPANmaster,Ymaster(:,2,4),'Color',color2,'LineWidth',5) 
ln1 = line([170 170], [0 200], 'Color', 'k', 'LineWidth', 3, 'LineStyle',':')
title('Initial Stimulation: [IFNg] = 0ng/mL; [IL-4] = 90ng/mL', 'FontSize',12)
xlabel('Time (wk)','FontSize',12)
xlim([0 1680])
xticks([0 96 168 336 504 672 840 1008 1176 1344 1512 1680])
xticklabels({'Stim+','Stim-','1','2','3','4','5','6','7','8','9','10'})
ylim([0 200])
yticks([0 50 100 150 200])
annotation('textbox', [0.12 0 0.5 0.5], 'String', 'f', 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
hold off

%% SUPP FIGURE 3a,b and 4 (Data): Sensitivity analyses near parameter points of interest
% Set parameters to reflect point of interest: symmetric parameters, cpmL =
% 10^6 or 10^9.
% Note that for wiggle >= 0.8, "unique" rounds to 1 decimal place instead of 0
Symmetry = true;    % Do you want to force parameters to be symmetric or not?
wiggle = 0.1;       % What fractional deviation from the parameter point-of-interest do you want to allow?
NumSamples = 1000;
sol_thresh = 1e-3;
sim_thresh = 1e-3;
searchGrid = linspace(0,10,11);
options = optimset('Display','off');
SAD = NaN(NumSamples, length(pars), 2);
SAD(:,:,1) = (rand(NumSamples, length(pars)).*(wiggle*2) + (1-wiggle)).*repmat(pars, NumSamples, 1);
if Symmetry
    SAD(:, 2:2:end, 1) = SAD(:, 1:2:end-1, 1);
end
for k = 1:NumSamples
    solutions = zeros(length(searchGrid)^2,5);
    for i = 1:length(searchGrid)    % For every starting point in the defined grid, locate the nearest fixpoint candidate.
        for j = 1:length(searchGrid)
            startPoint = [searchGrid(i) searchGrid(j)*3 searchGrid(i) searchGrid(j)*3];
            [A, B] = fsolve(@(x)root4dSUPP(x,SAD(k,:,1),hval),startPoint, options);
            solutions((i-1)*length(searchGrid)+j,1:4) = A.^2;       % Record the candidate fixpoint's 4-d coordinates...
            solutions((i-1)*length(searchGrid)+j,5) = sum(abs(B));  % ... and how close all 4 rates of change are to 0 ("B score").
        end
    end
    solutions(solutions(:,5)>=sol_thresh,:) = []; % Remove candidate fixpoints for which all 4 rates of change are not close enough to 0.
    [~, ~, ic] = unique(round(solutions(:,1:4),1), 'rows'); % Group which rows would be identical if rounded to the nearest whole #.
    for i = 1:max(ic)                                       % From each group, find the lowest B score. Flag the rest to be removed.
        solutions(solutions(:,5)>min(solutions(ic==i,5)) & ic==i, 5) = 1;
    end
    solutions(solutions(:,5)==1,:) = [];
    for i = 1:length(solutions(:,1))        % Now, analyze the Jacobian's eigenvalues for the candidate fixpoints.
        if any(real(nondim4dEIGS(solutions(i,1:4), SAD(k,:,1), hval)) >= 0) 
            solutions(i,5) = 1;             % Mark and subsequently remove unstable fixpoints.
        end
    end
    solutions(solutions(:,5)==1,:) = [];
    num_eq = NaN(length(solutions(:,1)), 4);
    for i = 1:length(solutions(:,1))      % In some cases, duplicates may persist (e.g. coordinates 6.49 & 6.51 appear different rounded.
        X0 = solutions(i,1:4);            % We can remove duplicates based on seeing them converge to identical equilibria numerically.
        F = @(t,x) [(SAD(k,1,1) + (x(1)^hval(1))/(SAD(k,3,1)^hval(1)+x(1)^hval(1)))/(1+x(2)^hval(2)) + ((SAD(k,5,1)*x(3)^hval(3))/(SAD(k,7,1)^hval(3)+x(3)^hval(3)))*(SAD(k,10,1)^hval(4))/(SAD(k,10,1)^hval(4)+x(4)^hval(4)) - SAD(k,11,1)*x(1);
            ((SAD(k,2,1) + (x(2)^hval(1))/(SAD(k,4,1)^hval(1)+x(2)^hval(1)))/(1+x(1)^hval(2)) + ((SAD(k,6,1)*x(4)^hval(3))/(SAD(k,8,1)^hval(3)+x(4)^hval(3)))*(SAD(k,9,1)^hval(4))/(SAD(k,9,1)^hval(4)+x(3)^hval(4)) - SAD(k,12,1)*x(2))*SAD(k,21,1);
            ((SAD(k,13,1)*x(1)^hval(5))/(SAD(k,15,1)^hval(5)+x(1)^hval(5)))*((SAD(k,18,1)^hval(6))/(SAD(k,18,1)^hval(6)+x(2)^hval(6)))/(1+x(4)^hval(7)) - SAD(k,19,1)*x(3);
            (((SAD(k,14,1)*x(2)^hval(5))/(SAD(k,16,1)^hval(5)+x(2)^hval(5)))*((SAD(k,17,1)^hval(6))/(SAD(k,17,1)^hval(6)+x(1)^hval(6)))/(1+x(3)^hval(7)) - SAD(k,20,1)*x(4))*SAD(k,21,1)];
        G = @(t,x) [0; 0; 0; 0];
        TSPAN = 0:0.1:2000;
        options2 = sdeset('SDEType','Ito');
        x = sde_euler(F,G,TSPAN,X0,options2);
        num_eq(i,:) = x(length(TSPAN),:);
        if sum(abs(solutions(i,1:4) - x(length(TSPAN), :))) > sim_thresh % If the simulation's end point doesn't equal the start point,
            solutions(i,5) = 1;                                          % flag that false solution to be removed.
        end
    end
    [~, ~, ic] = unique(round(num_eq,1), 'rows'); % Group which numerically found equilibria are equivalent if rounded to the nearest whole #.
    for i = 1:max(ic)                     % From each group, find the lowest B score. Flag the rest to be removed.
        solutions(solutions(:,5)>min(solutions(ic==i,5)) & ic==i, 5) = 1;
    end
    solutions(solutions(:,5)==1,:) = [];
    SAD(k, length(SAD(k,:,2)), 2) = length(solutions(:,1));
    for i = 1:length(solutions(:,1))
        SAD(k,i*4-3:i*4,2) = solutions(i,1:4);
    end
    if k/NumSamples == round(k/NumSamples,2)
        readout = ['Progress Made: ', num2str(100*k/NumSamples), '%'];
        disp(readout);
    end
end
%save('SensitivityAnalysis_LoDens_Wig1', 'SAD');
%From here, analyses and figures are easily executed in R.

%% SUPP FIGURE 3c,d (Data):
% Set parameters to reflect point of interest: symmetric parameters, cpmL =
% 10^6 or 10^9 for supp fig 3c or d, respectively.
% Be sure to change the name of the results array, below, to c or d.
Symmetry = true;    % Do you want to force parameters to be symmetric or not?
numSamples = 101;
SuppFig3dResults = NaN(numSamples, numSamples, 4);
hval = [hp hx hs hz ha hr hu];
Dvals = linspace(D1*0.95, D1*1.05, numSamples);
Fvals = linspace(F1*0.95, F1*1.05, numSamples);
for i = 1:numSamples
    for j = 1:numSamples
        pars = [B1 B2 H1 H2 Fvals(j) Fvals(j) E1 E2 W1 W2 Dvals(i) Dvals(i) G1 G2 V1 V2 Q1 Q2 L1 L2 Asym];
        F = @(t,Y) [(pars(1) + ((Y(1))^hval(1))/(pars(3)^hval(1)+(Y(1))^hval(1)))/(1+(Y(2))^hval(2)) + ((pars(5)*(Y(3))^hval(3))/(pars(7)^hval(3)+(Y(3))^hval(3)))*(pars(10)^hval(4))/(pars(10)^hval(4)+(Y(4))^hval(4)) - pars(11)*(Y(1));
            ((pars(2) + ((Y(2))^hval(1))/(pars(4)^hval(1)+(Y(2))^hval(1)))/(1+(Y(1)^2)^hval(2)) + ((pars(6)*(Y(4))^hval(3))/(pars(8)^hval(3)+(Y(4))^hval(3)))*(pars(9)^hval(4))/(pars(9)^hval(4)+(Y(3))^hval(4)) - pars(12)*(Y(2)))*pars(21);
            ((pars(13)*(Y(1))^hval(5))/(pars(15)^hval(5)+(Y(1))^hval(5)))*((pars(18)^hval(6))/(pars(18)^hval(6)+(Y(2))^hval(6)))/(1+(Y(4))^hval(7)) - pars(19)*(Y(3));
            (((pars(14)*(Y(2))^hval(5))/(pars(16)^hval(5)+(Y(2))^hval(5)))*((pars(17)^hval(6))/(pars(17)^hval(6)+(Y(1))^hval(6)))/(1+(Y(3))^hval(7)) - pars(20)*(Y(4)))*pars(21)];
        G = @(t,Y) [nTF1*Y(1);
            nTF2*Y(2);
            nCY1*Y(3);
            nCY2*Y(4)];
        TSPAN = 0:0.1:1000;
        Y0 = [0 0 0 0];
        options = sdeset('SDEType','Ito');
        Y = sde_euler(F,G,TSPAN,Y0,options);
        SuppFig3dResults(i,j,:) = Y(length(TSPAN), :);
        if j == numSamples
            fprintf('Percent Completed: %d\n', round(100*i/numSamples));
        end
    end
end
%save('SuppFig3dResults_FOR_MS2','SuppFig3dResults');
%load('SuppFig3dResults_FOR_MS2');

%% SUPP FIGURE 3c (Plot)
surf(Dvals./D1, Fvals./F1, SuppFig3cResults(:,:,1)./SuppFig3cResults(51,51,1));
xlim([0.95 1.05]);
xticks([0.95 1 1.05]);
xticklabels({'95%', '100%', '105%'});
xlabel('Relative Value of F');
ylim([0.95 1.05]);
ylabel('Relative Value of D');
yticks([0.95 1 1.05]);
yticklabels({'95%', '100%', '105%'});
zlim([0.95 1.05]);
zticks([0.95 1 1.05]);
zticklabels({'95%', '100%', '105%'});
zlabel('Relative Equilibrium Position');
colormap('jet');

%% SUPP FIGURE 3d (Plot)
surf(Dvals./D1, Fvals./F1, SuppFig3dResults(:,:,1)./SuppFig3dResults(51,51,1));
xlim([0.95 1.05]);
xticks([0.95 1 1.05]);
xticklabels({'95%', '100%', '105%'});
xlabel('Relative Value of F');
ylim([0.95 1.05]);
ylabel('Relative Value of D');
yticks([0.95 1 1.05]);
yticklabels({'95%', '100%', '105%'});
zlim([0.9 1.1]);
zticks([0.90 0.95 1 1.05 1.1]);
zticklabels({'90%', '95%', '100%', '105%', '110%'});
zlabel('Relative Equilibrium Position');
colormap('jet');
