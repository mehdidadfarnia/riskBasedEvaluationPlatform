%% Manufacturing Simulator Input Script, V6
%  Varying Machine Fail Times
%
%  Instructions
%  ------------
%
%  Testing different machine fail times, dependent on underlying 
%  fail time probability distributions
%
%  Parts 1-3 go through both Manufacturing Simulators & Visualizations
        %  Part 2 has 10 Steps
%  Parts 4-7 go through risk evaluation of each algorithm
%
%  This script calls on the following functions:
%       FunMfgSimulator_v02.m (which in turn calls on FunMeanFilt_v02.m)
%       FunMakePlots_v02.m
%       randfixedsum.m

%% ============== Part 1: Initialization & Machining Paths ================
% Step 1: Initialization
%clearvars -except Experiment;
clear all; clc; close all

% Step 2: Defining the machining paths
%  We start by defining the setup of the manufacturing system 


%Machines that perform Class A: Drilling operations
ClassAOpMach = [1 2 3]; 
%Machines that perform Class B: Boring operations
ClassBOpMach = [4 5 6 7 8]; 
%Machines that perform Class C: Laser Engraving operations
ClassCOpMach = [9 10];

permA=perms(ClassAOpMach); permA=permA(:,1); permA=unique(permA,'rows');
permB=perms(ClassBOpMach); permB=permB(:,1:2); permB=unique(permB,'rows');
permC=perms(ClassCOpMach); permC=permC(:,1); permC=unique(permC,'rows');

indMatLen=length(permA)*length(permB)*length(permC);
indMat=ones(indMatLen,3).*repmat((1:indMatLen)',1,3);

indMat(:,1)=(rem((indMat(:,1)-1),length(permA)))+1;
indMat(:,2)=(rem((indMat(:,2)-1-indMat(:,1)+1)/length(permA),length(permB)))+1;
mMat=(indMat(:,3)-1-indMat(:,1)+1)/length(permA);
indMat(:,3)=(rem((mMat-indMat(:,2)+1)/length(permB),length(permC)))+1;

permAmat=permA(indMat(:,1));
permBmat=permB(indMat(:,2),:);
permCmat=permC(indMat(:,3));

% Due to the order of machines in a path:
AllPossiblePaths=[permAmat permBmat(:,1) permCmat permBmat(:,2)];

%Select the paths considered by simulation:
NumRows=120;
SelRows=randperm(size(AllPossiblePaths,1),NumRows);
SelUniquePaths=AllPossiblePaths(SelRows,:);

% Now we can create a Unique Paths matrix with:
UniquePaths=mat2cell(SelUniquePaths,ones(1,size(SelUniquePaths,1)),size(SelUniquePaths,2));



%% =============== Part 2: Include Fail Time Probability Distributions  ================
% Example to incorporate probability distributions of machine failure times

% Look for failure times for our 10 machines, with different sets of 
% underlying fail-time probability distributions


% =============================
% Step #1: Initializations

% Total number of available machines
NumNodes = max(horzcat(UniquePaths{:})); 

% Simulation/operation duration is 45 minutes (2700 seconds)
NumParts3=2700;  

% Added value by each 'good' machine
ProdVal3 = 2;

% Subtracted value by each 'bad' machine (post-degradation)
BadLevel3 = -3;

% Acceptable limit of product value
GoodLim=0;


% ============================
% Step #2: Each class of machines (Class A, Class B, Class C) has 
% its own life distribution, by which where the failure times come from

% Class A-operation Machine Distributions:
theta=4800;
y=floor(sort(exprnd(theta*ones(1,length(ClassAOpMach)))));
fy=(1/theta)*exp(-y./theta);
ExpShapes=20;
yShape=floor(sort(exprnd(theta*ones(1,ExpShapes*length(ClassAOpMach))))); % ceil() or floor() is optional
fyShape=(1/theta)*exp(-yShape./theta);


% Class B-operation Machine Distributions:
% 1st Weibull distribution set:
    %Shape parameter, beta, unitless:
zwblB=1.3; 
    %Scale parameter, alpha, in seconds: 
    % close def'n: at what point of time does the curve has the smallest 
    % derivative determines "spread" of the distribution
zwblA=14400; 
z=floor(sort(wblrnd(zwblA,zwblB,[1,length(ClassBOpMach)]))); 
fz=(zwblB/(zwblA^zwblB)).*(z.^(zwblB-1)).*exp(-((z./zwblA).^(zwblB)));
% For a more clear shape of the distribution:
WeibullShapes=20;
zShape=floor(sort(wblrnd(zwblA,zwblB,[1,WeibullShapes*length(ClassBOpMach)]))); 
fzShape=(zwblB/(zwblA^zwblB)).*(zShape.^(zwblB-1)).*exp(-((zShape./zwblA).^(zwblB)));

% Class C-operation Machine Distributions:
% 2nd Weibull distribution set:
%     %Shape parameter, beta, unitless:
xwblB=4; 
    %Scale parameter, alpha, in seconds: 
    % close def'n: at what point of time does the curve has the smallest 
    % derivative determines "spread" of the distribution
xwblA=3800; 
x=floor(sort(wblrnd(xwblA,xwblB,[1,length(ClassCOpMach)]))); 
fx=(xwblB/(xwblA^xwblB)).*(x.^(xwblB-1)).*exp(-((x./xwblA).^(xwblB)));
% For a more clear shape of the distribution:
% WeibullShapes=20;
xShape=floor(sort(wblrnd(xwblA,xwblB,[1,WeibullShapes*length(ClassCOpMach)]))); 
fxShape=(xwblB/(xwblA^xwblB)).*(xShape.^(xwblB-1)).*exp(-((xShape./xwblA).^(xwblB)));


% Plot of the machines by underlying distribution of fail times
figure(22); plot(y,fy,'ro',yShape,fyShape,'r--',...
    z,fz,'bo',zShape,fzShape,'b--',...
    x,fx,'go',xShape,fxShape,'g--')
hold on;
xline(NumParts3);
hold off;
legend('Class A: Exponential Machines','Exponential Shape',...
    'Class B: Weibull-1 Machines','Weibull-1 Shape',...
    'Class C: Weibull-2 Machines','Weibull-2 Shape',...
    'Simulation Cutoff','location','northeast')
title({'Visualizing the failure time distributions',...
    'for different "types" of machines'})



% ============================
% % Alternate Step #2: Use a distribution generator to get failures at different times 
% % for all N machines (nodes)
% 
% % Set A proportion of machines will have a Weibull distribution for their
% % fail times
% setA=0.41;
% 
% % Set B proportion of machines will have an Exponential distribution for their
% % fail times
% setB=0.59;
% 
% % Exponential distribution set:
%     % with known 1 failure per 4800 seconds (one hour, 20 minutes)
% theta=4800;
% y=floor(sort(exprnd(theta*ones(1,ceil(NumNodes*setB))))); % ceil() or floor() is optional
% fy=(1/theta)*exp(-y./theta);
% % For a more clear shape of the distribution:
% ExpShapes=20;
% yShape=floor(sort(exprnd(theta*ones(1,ceil(ExpShapes*NumNodes*setB))))); % ceil() or floor() is optional
% fyShape=(1/theta)*exp(-yShape./theta);


% % Weibull distribution set:
%     %Shape parameter, beta, unitless:
% zwblB=1.3; 
%     %Scale parameter, alpha, in seconds: 
%     % close def'n: at what point of time does the curve has the smallest 
%     % derivative determines "spread" of the distribution
% zwblA=14400; 
% z=floor(sort(wblrnd(zwblA,zwblB,[1,floor(NumNodes*setA)]))); 
% fz=(zwblB/(zwblA^zwblB)).*(z.^(zwblB-1)).*exp(-((z./zwblA).^(zwblB)));
% % For a more clear shape of the distribution:
% WeibullShapes=20;
% zShape=floor(sort(wblrnd(zwblA,zwblB,[1,floor(WeibullShapes*NumNodes*setA)]))); 
% fzShape=(zwblB/(zwblA^zwblB)).*(zShape.^(zwblB-1)).*exp(-((zShape./zwblA).^(zwblB)));


% % Plot of the machines by underlying distribution of fail times
% figure(22); plot(y,fy,'ro',yShape,fyShape,'r--',z,fz,'bo',zShape,fzShape,'b--')
% hold on;
% xline(NumParts3);
% hold off;
% legend('Exponential Machines','Exponential Shape','Weibull Machines','Weibull Shape','Simulation Cutoff','location','northeast')
% title({'Visualizing the failure time distributions','for different "types" of machines'})


% ============================
% Step #3: Not all machines start as brand new machines; some may have been used before
% y has fail times for Exponential Machines
% z has fail times for Weibull Machines
% Let us assume all the Weibull Machines have already been used for 4000 time units:
yUsedShift=200;
zUsedShift=4000;
xUsedShift=400;

y=y-yUsedShift;
z=z-zUsedShift;
x=x-xUsedShift;
allTimesTemp=[y z x];
allTimesTemp(allTimesTemp<1)=randi(200,size(allTimesTemp(allTimesTemp<0)));

y=allTimesTemp(1:length(y));
z=allTimesTemp(1+length(y):length(y)+length(z));
x=allTimesTemp(1+length(y)+length(z):end);


% ============================
% Step #4: If some fail times are outside simulation/operation duration
% (consider 45 minutes or 2700 seconds), drop them

% Vector of all machine failure times, regardless of distribution

AllMachineFailTimes=[y,z,x]; %Order of Class A, Class B, Class C !!!!!!!!!!!!!!!!!!!!!!!!!!
%AllMachineFailTimes=sort([y,z,x]);


% Dropping fail times that do not fall within simulation
DegradeStart3=AllMachineFailTimes(AllMachineFailTimes<=NumParts3);

%Store fail times that happen after the simulation
DegaradePostSim=AllMachineFailTimes(AllMachineFailTimes>NumParts3);



% ============================
% Step #5: Define Set of Bad Machines. 
% Each column is a machine that has gone bad and its 
% degradation start time 
%BadSet3 = sort(randperm(NumNodes,length(DegradeStart3)));
BadSet3 = find(ismember(AllMachineFailTimes,DegradeStart3));
MachineFailureTimes=[BadSet3;DegradeStart3];

% Also defining the set of Machines that didn't go bad during the
% simulated duration of operation
PostSimFailSet=1:length(AllMachineFailTimes);
PostSimFailSet(BadSet3)=[];
MachinePostSimFailTimes=[PostSimFailSet;DegaradePostSim];



% ============================
% Step #6: Differentiate machines based on their distributions (class with
% exponential distribution, class with weibull distribution, and so on)

MachinesFailsInSeq=[sortrows(MachineFailureTimes',2)', sortrows(MachinePostSimFailTimes',2)'];

% Class A: Exponential Machines
ClassAidxFailSeq=ismember(MachinesFailsInSeq(2,:),y);
ClassAMachinesAndFailTimes=[sortrows(MachinesFailsInSeq(:,ClassAidxFailSeq)',1)';ClassAOpMach]; %Last row provides us with original machine ID #s
ClassAidx=zeros(1,NumNodes);
ClassAidx(MachinesFailsInSeq(1,ClassAidxFailSeq))=1;

% Class B: Weibull Type-1 Machines
ClassBidxFailSeq=ismember(MachinesFailsInSeq(2,:),z);
ClassBMachinesAndFailTimes=[sortrows(MachinesFailsInSeq(:,ClassBidxFailSeq)',1)';ClassBOpMach]; %Last row provides us with original machine ID #s
ClassBidx=zeros(1,NumNodes);
ClassBidx(MachinesFailsInSeq(1,ClassBidxFailSeq))=1;

% Class C: Weibull Type-2 Machines
ClassCidxFailSeq=ismember(MachinesFailsInSeq(2,:),x);
ClassCMachinesAndFailTimes=[sortrows(MachinesFailsInSeq(:,ClassCidxFailSeq)',1)';ClassCOpMach]; %Last row provides us with original machine ID #s
ClassCidx=zeros(1,NumNodes);
ClassCidx(MachinesFailsInSeq(1,ClassCidxFailSeq))=1;


% ============================
% Step #7: FMEA-Derived Values for Occurrence, or likelihood of a failure mode:

% h(t) or the hazard function
% Integral of h(t) (cumulative hazard function)
% a(t) or the proportion of hazard function per failure mode
% Beta, the conditional probability that a failure mode will 
            % actually cost that much given the failure mode occurs



tiempo=1:NumParts3;

% Find an instantaneous time-dependent failure 
% each class of machines each have an instantaneous component
% failure rate, e.g. hY(t) and hZ(t), respectively

% Class A-operation Machine, Exponential instantaneous failure rate:
hY=(1/theta).*ones(1,length(tiempo+yUsedShift));
% Class B-operation Machine, Weibull-1 instantaneous failure rate:
hZ=(zwblB/zwblA)*(((tiempo+zUsedShift)./zwblA).^(zwblB-1)); 
% Class C-operation Machine, Weibull-2 instantaneous failure rate:
hX=(xwblB/xwblA)*(((tiempo+xUsedShift)./xwblA).^(xwblB-1)); 


%
% Machines in different sets/classes have different failure rates (known as the hazard 
% function), but different failure modes contribute to it. At any moment, 
% different failure modes make up different ratios of the hazard function.

% Class A-operation Machines (Exponential distribution):
numFailModesClassAy=1;
aY= ones(numFailModesClassAy,length(tiempo));

% Class B-operation Machines (Weibull-1 distribution):
numFailModesClassBz=3;
minaZ=0;
maxaZ=1;
aZ=randfixedsum(numFailModesClassBz,length(tiempo),1,minaZ,maxaZ); %Using a function from Mathworks File Exchange


% Class C-operation Machines (Weibull-2 distribution):
numFailModesClassCx=2;
minaX=0;
maxaX=1;
aX=randfixedsum(numFailModesClassCx,length(tiempo),1,minaX,maxaX); %Using a function from Mathworks File Exchange


% % For machines in Set A (Weibull distribution):
% numFailModesClassBz=3;
% minaZ=0;
% maxaZ=1;
% aZ=randfixedsum(numFailModesClassBz,length(tiempo),1,minaZ,maxaZ); %Using a function from Mathworks File Exchange

% % Likewise for machines in Set B (Exponential distribution):
% numFailModesClassAy=1;
% aY= ones(numFailModesClassAy,length(tiempo));


figure(23);
plot(tiempo,hY,'r--','DisplayName','Class A Machines - Exponential Hazard Function')
ylabel('Instantaneous Failure Rate & Their Ratios');
xlabel('Time');
hold on;
plot(tiempo,hZ,'b--','DisplayName','Class B Machines - Weibull-1 Hazard Function')
hold on;
plot(tiempo,hX,'g--','DisplayName','Class C Machines - Weibull-2 Hazard Function')


for jj=1:numFailModesClassAy
    clr2=jj/numFailModesClassAy;
    plot(tiempo,aY(jj).*hY,'LineStyle',':','Color',[clr2 0 0],'DisplayName',strcat('Class A - Exponential Ratio ', num2str(jj)));
end

for ii=1:numFailModesClassBz
    clr1=ii/numFailModesClassBz;
    plot(tiempo,aZ(ii,:).*hZ,'LineStyle',':','Color',[0 0 clr1],'DisplayName',strcat('Class B - Weibull-1 Ratio ', num2str(ii)));
%     get(legend(gca),); % get legend from current axes.
%    legend(strcat('Weibull Ratio', num2str(ii)))
end

for kk=1:numFailModesClassCx
    clr3=kk/numFailModesClassCx;
    plot(tiempo,aX(kk).*hX,'LineStyle',':','Color',[0 clr3 0],'DisplayName',strcat('Class C - Weibull-2 Ratio ', num2str(kk)));
end


legend('location','east')

hold off



% Beta represents the conditional probability that a failure will result in 
    % identified severities, given that the failure mode actually occurs
BetaY=unifrnd(.5,1,[numFailModesClassAy,1]);
BetaZ=unifrnd(.75,1,[numFailModesClassBz,1]);
BetaX=unifrnd(.6,1,[numFailModesClassCx,1]);
figure(24);
hold on
for jj=1:numFailModesClassAy
    clr2=jj/numFailModesClassAy;
    plot(tiempo,BetaY(jj).*ones(1,length(tiempo)),'LineStyle',':','Color',[clr2 0 0],'DisplayName', ...
        strcat('Beta parameter for Class A - Exponential Failure Mode: ', num2str(jj)));
end
for ii=1:numFailModesClassBz
    clr1=ii/numFailModesClassBz;
    plot(tiempo,BetaZ(ii).*ones(1,length(tiempo)),'LineStyle','-','Color',[0 0 clr1],'DisplayName', ...
        strcat('Beta parameter for Class B - Weibull-1  Failure Mode: ', num2str(ii)));
end
for kk=1:numFailModesClassCx
    clr3=kk/numFailModesClassCx;
    plot(tiempo,BetaX(kk).*ones(1,length(tiempo)),'LineStyle','--','Color',[0 clr3 0],'DisplayName', ...
        strcat('Beta parameter for Class C - Weibull-2 Failure Mode: ', num2str(kk)));
end

legend('location','east')
hold off



tiempoAy=1:NumParts3+yUsedShift;
hYTot=(1/theta).*ones(1,length(tiempoAy));
hYinteg=cumtrapz(hYTot);
hYinteg=hYinteg(yUsedShift+1:NumParts3+yUsedShift);

tiempoBz=1:NumParts3+zUsedShift;
hZTot=(zwblB/zwblA)*(((tiempoBz)./zwblA).^(zwblB-1));
hZinteg=cumtrapz(hZTot);
hZinteg=hZinteg(zUsedShift+1:NumParts3+zUsedShift);

tiempoCx=1:NumParts3+xUsedShift;
hXTot=(xwblB/xwblA)*(((tiempoCx)./xwblA).^(xwblB-1));
hXinteg=cumtrapz(hXTot);
hXinteg=hXinteg(xUsedShift+1:NumParts3+xUsedShift);


% tiempoAy=1:NumParts3+yUsedShift;
% hYTot=(1/theta).*ones(1,length(tiempoAy));
% hYinteg=cumtrapz(hYTot);
% hYinteg=hYinteg(yUsedShift+1:NumParts3+yUsedShift);
% 
% tiempoBz=1:NumParts3+zUsedShift;
% hZTot=(zwblB/zwblA)*(((tiempoBz)./zwblA).^(zwblB-1));
% hZinteg=cumtrapz(hZTot);
% hZinteg=hZinteg(zUsedShift+1:NumParts3+zUsedShift);


figure(25);
plot(tiempo,hYinteg,'r--','DisplayName','Cumulative Class A - Exponential Hazard Function')
title('Instantaneous Failure Rate & Their Ratios');
xlabel('Time');
hold on;
plot(tiempo,hZinteg,'b--','DisplayName','Cumulative Class B - Weibull-1 Hazard Function')
hold on;
plot(tiempo,hXinteg,'b--','DisplayName','Cumulative Class C - Weibull-2 Hazard Function')

hold off;
legend('location','northeast')




% ============================
% Step #8: FMEA-Derived Values for Severity, also known as Costs:


% SEVERITY - Costs of each failure mode may in the future be 
% further divided between: 
    %(cost of failed component) + (cost of down-the-line consequences of failed
    %component) + (cost of down-the-line consequences of the failure mode)

CostsAy=unifrnd(50,200,[numFailModesClassAy,1]); %Class A, Exponential Machines, 1 failure mode
CostsBz=unifrnd(150,600,[numFailModesClassBz,1]); %Class B, Weibull-1 Machines, 3 failure modes
CostsCx=unifrnd(10,100,[numFailModesClassCx,1]); %Class C, Weibull-2 Machines, 2 failure modes

figure(26);
title('Severity');
hold on
for jj=1:numFailModesClassAy
    clr2=jj/numFailModesClassAy;
    plot(tiempo,CostsAy(jj).*ones(1,length(tiempo)),'LineStyle',':','Color',[clr2 0 0],'DisplayName', ...
        strcat('Costs for Class A, Exponential Machine Failure Mode: ', num2str(jj)));
end

for ii=1:numFailModesClassBz
    clr1=ii/numFailModesClassBz;
    plot(tiempo,CostsBz(ii).*ones(1,length(tiempo)),'LineStyle','-','Color',[0 0 clr1],'DisplayName', ...
        strcat('Costs for Class B, Weibull-1 Machine Failure Mode: ', num2str(ii)));
end

for kk=1:numFailModesClassCx
    clr3=kk/numFailModesClassCx;
    plot(tiempo,CostsCx(kk).*ones(1,length(tiempo)),'LineStyle','--','Color',[0 clr3 0],'DisplayName', ...
        strcat('Costs for Class C, Weibull-2 Machine Failure Mode: ', num2str(kk)));
end
legend('location','east')
hold off



% ============================
% Step #9: Calculating Occurrence:

figure(27);
title('Occurrence, as Criticality(Cm)');
hold on
BetaYFun=BetaY.*ones(numFailModesClassAy,length(tiempo)); %1 by 2700
BetaZFun=BetaZ.*ones(numFailModesClassBz,length(tiempo)); %3 by 2700
BetaXFun=BetaX.*ones(numFailModesClassCx,length(tiempo)); %2 by 2700

hYintegFun=repmat(hYinteg,numFailModesClassAy,1);
hZintegFun=repmat(hZinteg,numFailModesClassBz,1);
hXintegFun=repmat(hXinteg,numFailModesClassCx,1);

CriticalityAy = BetaYFun.*aY.*hYintegFun;%Failure Mode Criticalities for Class A, Exponential Machines
CriticalityBz = BetaZFun.*aZ.*hZintegFun;%Failure Mode Criticalities for Class B, Weibull-1 Machines
CriticalityCx = BetaXFun.*aX.*hXintegFun;%Failure Mode Criticalities for Class C, Weibull-2 Machines

for jj=1:numFailModesClassAy
    clr2=jj/numFailModesClassAy;
    plot(tiempo,CriticalityAy(jj,:),'LineStyle',':','Color',[clr2 0 0],'DisplayName', ...
        strcat('Criticality for Class A, Exponential Machine Failure Mode: ', num2str(jj)));
end
for ii=1:numFailModesClassBz
    clr1=ii/numFailModesClassBz;
    plot(tiempo,CriticalityBz(ii,:),'LineStyle','-','Color',[0 0 clr1],'DisplayName', ...
        strcat('Criticality for Class B, Weibull-1 Machine Failure Mode: ', num2str(ii)));
end
for kk=1:numFailModesClassCx
    clr3=kk/numFailModesClassCx;
    plot(tiempo,CriticalityCx(kk,:),'LineStyle','--','Color',[0 clr3 0],'DisplayName', ...
        strcat('Criticality for Class C, Weibull-2 Machine Failure Mode: ', num2str(kk)));
end
legend('location','east')
hold off




% ============================
% Step #10: Calculating Risks:


figure(28);
title('Failre Mode Risks = Severity * Occurrence');
hold on
CostAyFun=CostsAy.*ones(numFailModesClassAy,length(tiempo));
CostBzFun=CostsBz.*ones(numFailModesClassBz,length(tiempo));
CostCxFun=CostsCx.*ones(numFailModesClassCx,length(tiempo));


%Risk of failure modes of a Class A, Exponential Machine
RiskAyFun = CostAyFun.*CriticalityAy;
%Risk of failure modes of a Class B, Weibull-1 Machine
RiskBzFun = CostBzFun.*CriticalityBz;
%Risk of failure modes of a Class C, Weibull-2 Machine
RiskCxFun = CostCxFun.*CriticalityCx;

for jj=1:numFailModesClassAy
    clr2=jj/numFailModesClassAy;
    plot(tiempo,RiskAyFun(jj,:),'LineStyle',':','Color',[clr2 0 0],'DisplayName', ...
        strcat('Inherent Risk for Class A, Exponential Machine Failure Mode: ', num2str(jj)));
end
for ii=1:numFailModesClassBz
    clr1=ii/numFailModesClassBz;
    plot(tiempo,RiskBzFun(ii,:),'LineStyle','-','Color',[0 0 clr1],'DisplayName', ...
        strcat('Inherent Risk for Class B, Weibull-1 Machine Failure Mode: ', num2str(ii)));
end
for kk=1:numFailModesClassCx
    clr3=kk/numFailModesClassCx;
    plot(tiempo,RiskCxFun(kk,:),'LineStyle','--','Color',[0 clr3 0],'DisplayName', ...
        strcat('Inherent Risk for Class C, Weibull-2 Machine Failure Mode: ', num2str(kk)));
end

legend('location','east')
hold off



ComponentRiskAy=sum(RiskAyFun,1);
ComponentRiskBz=sum(RiskBzFun,1);
ComponentRiskCx=sum(RiskCxFun,1);

figure(29);
plot(tiempo,ComponentRiskAy,'LineStyle',':','Color',[1 0 0],'DisplayName', ...
        'Inherent Risk for Class A, Exponential Machines (sum of failure modes)');
title('Inherent Component (Machine) Risk');
hold on;
plot(tiempo,ComponentRiskBz,'LineStyle','-','Color',[0 0 1],'DisplayName', ...
    'Inherent Risk for Class B, Weibull-1 Machines (sum of failure modes)');
hold on;
plot(tiempo,ComponentRiskCx,'LineStyle','--','Color',[0 1 0],'DisplayName', ...
    'Inherent Risk for Class C, Weibull-2 Machines (sum of failure modes)');
legend('location','east')
hold off


figure(30);
SystemRisk=sum(ClassAidx)*ComponentRiskAy + sum(ClassBidx)*ComponentRiskBz + +sum(ClassCidx)*ComponentRiskCx;
plot(tiempo,SystemRisk,'k--');
hold on;
title('Total System Risk (all the machines)');
%legend('location','east')
NodePresence=nan(1,NumNodes);
for ii=1:NumNodes
    NodePresence(1,ii)=length(find(horzcat(UniquePaths{:})==ii));
end
NodePresencePreNorm=NodePresence/max(NodePresence);
NodePresenceNorm=NodePresencePreNorm + (1-mean(NodePresencePreNorm));
%NodePresenceNorm=(NodePresence-min(NodePresence))/(max(NodePresence)-min(NodePresence))

WeightedSystemRisk=sum(NodePresenceNorm.*ClassAidx)*ComponentRiskAy + ...
    sum(NodePresenceNorm.*ClassBidx)*ComponentRiskBz + sum(NodePresenceNorm.*ClassCidx)*ComponentRiskCx;
plot(tiempo,WeightedSystemRisk,'LineStyle',':','Color',[0.25 0.25 0.25]);
legend('System (Scenario) Risk','Weighted System Risk','location','northwest');
hold off;




PathsMap = zeros(length(UniquePaths),NumNodes); 
for upi = 1:length(UniquePaths)
    PathsMap(upi,UniquePaths{upi}) = 1; 
end
PathRisk=zeros(length(UniquePaths),length(SystemRisk));



for kk=1:length(UniquePaths)
    PathSetAy=sum(PathsMap(kk,:).*ClassAidx);
    PathSetBz=sum(PathsMap(kk,:).*ClassBidx);
    PathSetCx=sum(PathsMap(kk,:).*ClassCidx);
    PathRisk(kk,:)=PathSetAy*ComponentRiskAy+PathSetBz*ComponentRiskBz +PathSetCx*ComponentRiskCx ;
end
PathNames = {};


figure(31);
%xline(NumParts3);
hold on
for ppi = 1:length(UniquePaths) %for each of the unique paths 
    plot(tiempo,PathRisk(ppi,:));% +unifrnd(0,10),'.'); 
    PathNames{ppi} = sprintf('Path #%i Risk',ppi);
end
xlabel('Time (s)');
legend(PathNames{:},'location','eastoutside')
hold off


%% =============== Part 3: Simulation & Visualization  ================


load('ExperimentSet.mat');

allExpNums=arrayfun(@(x) x.Input.ExpNum,Experiment);
 % Experiment ID
ExpNum=floor(allExpNums+1);

 % Simulator ID (adds to number of structures in ExperimentSet.mat
SimNum=numel(Experiment)+1;

% Manufacturing System Simulator
[Exp3]=FunMfgSimulator_v02(UniquePaths, NumParts3, ExpNum, ...
     GoodLim, BadSet3, ProdVal3, BadLevel3, DegradeStart3)


% Save to ExperimentSet.mat
    % Experiment(SimNum)=Exp3;
    % save(sprintf('ExperimentSet.mat'),'Experiment')

 % Generate Visualizations for the Simulator
 FunMakePlots_v02(Exp3.Output,3);

%% =============== Part 4: Risk Eval: Genetic Algorithm  ================

confMatGA=Exp3.Output.confMatGA;
%confMatFmM=Exp3.Output.confMatFmM;
GAoutput=Exp3.Output.B3mean;
%Example for Genetic Algorithm:
tpGA=confMatGA(1);
fnGA=confMatGA(2);
fpGA=confMatGA(3);
tnGA=confMatGA(4);
fprGA=fpGA/(tnGA+fpGA);
tprGA=tpGA/(tpGA+fnGA);
fnrGA=1-tprGA;

% For example, we have WeightedSystemRisk(end=2700)=576.
% That is we have a projected risk of 576 at the end of operation/duration
% But if we implement a batch algorithm, we get values for:
        % true positive (tp), false positive (fp), false negative (fn)
% We take actions for tp, fp
% We incur continuing costs for fn
% We do not do anything for tn

% tp: lets say it takes 5 bucks to fix machine A (tpACost), and 10 bucks to
        %fix machine B (tpBCost). The subset of machines that are true
        %positive and of type A is tpACount. Likewise, for machines of type
        %B, there is tpBCount. 
        % tpCostTotal = tpACount*tpACost + tpBCount*tpBCost
% fp: similar but for fp

GoodSet=1:10;
GoodSet(BadSet3)=[];
thresh=0.6; 
tpAy_GA=sum((GAoutput(BadSet3)<=thresh).*ClassAidx(BadSet3));
tpBz_GA=sum((GAoutput(BadSet3)<=thresh).*ClassBidx(BadSet3));
tpCx_GA=sum((GAoutput(BadSet3)<=thresh).*ClassCidx(BadSet3));

fpAy_GA=sum((GAoutput(GoodSet)<=thresh).*ClassAidx(GoodSet));
fpBz_GA=sum((GAoutput(GoodSet)<=thresh).*ClassBidx(GoodSet));
fpCx_GA=sum((GAoutput(GoodSet)<=thresh).*ClassCidx(GoodSet));

fnAy_GA=sum((GAoutput(BadSet3)>thresh).*ClassAidx(BadSet3));
fnBz_GA=sum((GAoutput(BadSet3)>thresh).*ClassBidx(BadSet3));
fnCx_GA=sum((GAoutput(BadSet3)>thresh).*ClassCidx(BadSet3));

tnAy_GA=sum((GAoutput(GoodSet)>thresh).*ClassAidx(GoodSet));
tnBz_GA=sum((GAoutput(GoodSet)>thresh).*ClassBidx(GoodSet));
tnCx_GA=sum((GAoutput(GoodSet)>thresh).*ClassCidx(GoodSet));

% Get the following from domain-level experts
    % For the exponential machines:
CostCheckAy=20;
CostFixAy=100;
CorrMaintAy=(CostCheckAy+CostFixAy)*tpAy_GA;
ContFailAy= fnAy_GA*sum(CostsAy.*mean(aY,2));
UnnecessaryMaintAy=(CostCheckAy)*fpAy_GA;
    % For the weibull-1 machines:
CostCheckBz=10;
CostFixBz=50;
CorrMaintBz=(CostCheckBz+CostFixBz)*tpBz_GA;
ContFailBz= fnBz_GA*sum(CostsBz.*mean(aZ,2));
UnnecessaryMaintBz=(CostCheckBz)*fpBz_GA;
    % For the weibull-2 machines:
CostCheckCx=30;
CostFixCx=30;
CorrMaintCx=(CostCheckCx+CostFixCx)*tpCx_GA;
ContFailCx= fnCx_GA*sum(CostsCx.*mean(aX,2));
UnnecessaryMaintCx=(CostCheckCx)*fpCx_GA;


%PostAlg_RISK_scenario = costTP*tpr + costFN*fnr + costFP*fpr;
GAPostAlg_RISK_scenario = (CorrMaintAy + CorrMaintBz + CorrMaintCx) + ...
                         (ContFailAy + ContFailBz + ContFailCx) + ...
                         (UnnecessaryMaintAy + UnnecessaryMaintBz + UnnecessaryMaintCx)
GANoAlg_RISK_scenario = mean(WeightedSystemRisk(end-20:end)) % avg the last few ones
GAValue_Alg = GANoAlg_RISK_scenario-GAPostAlg_RISK_scenario

% Focusing on false negatives, false positives of an applied algorithm
GAAlgorithmsRisk =(ContFailAy + ContFailBz + ContFailCx) + ...
                         (UnnecessaryMaintAy + UnnecessaryMaintBz + UnnecessaryMaintCx)

%% =============== Part 5: Risk Eval: Optimization Solver  ================

confMatFmM=Exp3.Output.confMatFmM;
FmMoutput=Exp3.Output.B2';
%Example for Genetic Algorithm:
tpFmM=confMatFmM(1);
fnFmM=confMatFmM(2);
fpFmM=confMatFmM(3);
tnFmM=confMatFmM(4);
fprFmM=fpFmM/(tnFmM+fpFmM);
tprFmM=tpFmM/(tpFmM+fnFmM);
fnrFmM=1-tprFmM;

% For example, we have WeightedSystemRisk(end=2700)=576.
% That is we have a projected risk of 576 at the end of operation/duration
% But if we implement a batch algorithm, we get values for:
        % true positive (tp), false positive (fp), false negative (fn)
% We take actions for tp, fp
% We incur continuing costs for fn
% We do not do anything for tn

% tp: lets say it takes 5 bucks to fix machine A (tpACost), and 10 bucks to
        %fix machine B (tpBCost). The subset of machines that are true
        %positive and of type A is tpACount. Likewise, for machines of type
        %B, there is tpBCount. 
        % tpCostTotal = tpACount*tpACost + tpBCount*tpBCost
% fp: similar but for fp

GoodSet=1:10;
GoodSet(BadSet3)=[];
thresh=0.6; 

tpAy_FmM=sum((FmMoutput(BadSet3)<=thresh).*ClassAidx(BadSet3));
tpBz_FmM=sum((FmMoutput(BadSet3)<=thresh).*ClassBidx(BadSet3));
tpCx_FmM=sum((FmMoutput(BadSet3)<=thresh).*ClassCidx(BadSet3));

fpAy_FmM=sum((FmMoutput(GoodSet)<=thresh).*ClassAidx(GoodSet));
fpBz_FmM=sum((FmMoutput(GoodSet)<=thresh).*ClassBidx(GoodSet));
fpCx_FmM=sum((FmMoutput(GoodSet)<=thresh).*ClassCidx(GoodSet));

fnAy_FmM=sum((FmMoutput(BadSet3)>thresh).*ClassAidx(BadSet3));
fnBz_FmM=sum((FmMoutput(BadSet3)>thresh).*ClassBidx(BadSet3));
fnCx_FmM=sum((FmMoutput(BadSet3)>thresh).*ClassCidx(BadSet3));

tnAy_FmM=sum((FmMoutput(GoodSet)>thresh).*ClassAidx(GoodSet));
tnBz_FmM=sum((FmMoutput(GoodSet)>thresh).*ClassBidx(GoodSet));
tnCx_FmM=sum((FmMoutput(GoodSet)>thresh).*ClassCidx(GoodSet));

    % For the exponential machines:
CostCheckAy=20;
CostFixAy=100;
FmMCorrMaintAy=(CostCheckAy+CostFixAy)*tpAy_FmM;
FmMContFailAy= fnAy_FmM*sum(CostsAy.*mean(aY,2));
FmMUnnecessaryMaintAy=(CostCheckAy)*fpAy_FmM;
% Get the following from domain-level experts
    % For the weibull-1 machines:
CostCheckBz=10;
CostFixBz=50;
FmMCorrMaintBz=(CostCheckBz+CostFixBz)*tpBz_FmM;
FmMContFailBz= fnBz_FmM*sum(CostsBz.*mean(aZ,2));
FmMUnnecessaryMaintBz=(CostCheckBz)*fpBz_FmM;
% Get the following from domain-level experts
    % For the weibull-2 machines:
CostCheckCx=30;
CostFixCx=30;
FmMCorrMaintCx=(CostCheckCx+CostFixCx)*tpCx_FmM;
FmMContFailCx= fnCx_FmM*sum(CostsCx.*mean(aX,2));
FmMUnnecessaryMaintCx=(CostCheckCx)*fpCx_FmM;



%PostAlg_RISK_scenario = costTP*tpr + costFN*fnr + costFP*fpr;
FmMPostAlg_RISK_scenario = (FmMCorrMaintAy + FmMCorrMaintBz + FmMCorrMaintCx) + ...
                         (FmMContFailAy + FmMContFailBz + FmMContFailCx) + ...
                         (FmMUnnecessaryMaintAy + FmMUnnecessaryMaintBz + FmMUnnecessaryMaintCx)
FmMNoAlg_RISK_scenario = mean(WeightedSystemRisk(end-10:end)) % avg the last few ones
FmMValue_Alg = FmMNoAlg_RISK_scenario-FmMPostAlg_RISK_scenario

% Focusing on false negatives, false positives of an applied algorithm
FmMAlgorithmsRisk = (FmMContFailAy + FmMContFailBz + FmMContFailCx) + ... 
                         (FmMUnnecessaryMaintAy + FmMUnnecessaryMaintBz + FmMUnnecessaryMaintCx)


%% ==Part 6: Risk Eval for Online Algorithm PQCI (Variable: ProdContNorm)==


PartQuality=Exp3.Output.PartQuality;
PartPath=Exp3.Output.PartPath;
%StartValue=Exp3.Output.PartPath;
StartValue=0; %temporary, until I redo another example where I run 
                % FunMfgSimulator_v02.m all over again.
Win = max(10,floor(NumParts3/100)); %Win=10 in the example with NumParts=500

%Moving Window of:
PWN = nan(Win,NumNodes);%Part Quality With Node
PWON = nan(Win,NumNodes);%Part Qualty Without Node
%Average of Windows
QWithN = nan(NumParts3,NumNodes);%Average Quality With Node
QWithoutN = nan(NumParts3,NumNodes);%Average Quality Without Node

ProdContThresh=-0.1;

ProdCont=nan(NumParts3,NumNodes);
ProdContNorm=nan(NumParts3,NumNodes);
tpProdCont=nan(NumParts3,1);
tnProdCont=nan(NumParts3,1);
fpProdCont=nan(NumParts3,1);
fnProdCont=nan(NumParts3,1);
fprProdCont=nan(NumParts3,1);
tprProdCont=nan(NumParts3,1);
fnrProdCont=nan(NumParts3,1);

tpAy_PC=nan(NumParts3,1);
tpBz_PC=nan(NumParts3,1);
tpCx_PC=nan(NumParts3,1);

fpAy_PC=nan(NumParts3,1);
fpBz_PC=nan(NumParts3,1);
fpCx_PC=nan(NumParts3,1);

fnAy_PC=nan(NumParts3,1);
fnBz_PC=nan(NumParts3,1);
fnCx_PC=nan(NumParts3,1);

tnAy_PC=nan(NumParts3,1);
tnBz_PC=nan(NumParts3,1);
tnCx_PC=nan(NumParts3,1);

PartQualityUB = StartValue+(length(UniquePaths{1})*ProdVal3);
if length(BadSet3)<length(UniquePaths{1})
    PartQualityLB = StartValue + (length(BadSet3)*BadLevel3) + ((length(UniquePaths{1}) - length(BadSet3))*ProdVal3);
else
    PartQualityLB = StartValue+length(UniquePaths{1})*BadLevel3;
end
ProdContRange = PartQualityUB - PartQualityLB;

NodeNames={};
for Ni = 1:NumNodes  
    NodeNames{Ni} = sprintf('Node #%i',Ni);
end

GoodSet=1:10;
GoodSet(BadSet3)=[];

% Get the following from domain-level experts
    % For the exponential machines:
CostCheckAy=20;
CostFixAy=100;
PCCorrMaintAy=nan(NumParts3,1);
PCContFailAy= nan(NumParts3,1);
PCUnnecessaryMaintAy=nan(NumParts3,1);

    % For the weibull-1 machines:
CostCheckBz=10;
CostFixBz=50;
PCCorrMaintBz=nan(NumParts3,1);
PCContFailBz= nan(NumParts3,1);
PCUnnecessaryMaintBz=nan(NumParts3,1);

    % For the weibull-2 machines:
CostCheckCx=30;
CostFixCx=30;
PCCorrMaintCx=nan(NumParts3,1);
PCContFailCx= nan(NumParts3,1);
PCUnnecessaryMaintCx=nan(NumParts3,1);


PCPostAlg_RISK_scenario=nan(NumParts3,1);
PCValue_Alg=nan(NumParts3,1);
PCAlgorithmsRisk=nan(NumParts3,1);


figure(32);
hold on;
for pp = 1:NumParts3
    %Place Part Quality in Appropriate Windowed Bins
    for Ni = 1:NumNodes
        %If Part Passed this Node (machining stage)
        if any(PartPath{pp}==Ni)
            PWN(:,Ni) = [PWN(2:end,Ni); PartQuality(pp)]; %moving window, row is the window of (10) products and column is each machine/node
        else
            %If Part Did Not Pass Through This Node
            PWON(:,Ni) = [PWON(2:end,Ni); PartQuality(pp)];
        end
    end
    QWithN(pp,:)   = nanmean(PWN);  %Average of each PWN column; each column characterizes a machine/node
    QWithoutN(pp,:)= nanmean(PWON); %Average of each PWON column; each column characterizes a machine/node
    ProdCont(pp,:)=QWithN(pp,:)-QWithoutN(pp,:);
    ProdContNorm(pp,:)=ProdCont(pp,:)/abs(ProdContRange);
    %plot(1:pp,ProdContNorm(1:pp,:))
    %pause(2);
    %clf;


    tpProdCont(pp)=sum(ProdContNorm(pp,BadSet3)<=ProdContThresh);
    fpProdCont(pp)=sum(ProdContNorm(pp,GoodSet)<=ProdContThresh);
    fnProdCont(pp)=sum(ProdContNorm(pp,BadSet3)>ProdContThresh);
    tnProdCont(pp)=sum(ProdContNorm(pp,GoodSet)>ProdContThresh);

    fprProdCont(pp)=fpProdCont(pp)/(tnProdCont(pp)+fpProdCont(pp));
    tprProdCont(pp)=tpProdCont(pp)/(tpProdCont(pp)+fnProdCont(pp));
    fnrProdCont(pp)=1-tprProdCont(pp);
    
    tpAy_PC(pp)=sum((ProdContNorm(pp,BadSet3)<=ProdContThresh).*ClassAidx(BadSet3));
    tpBz_PC(pp)=sum((ProdContNorm(pp,BadSet3)<=ProdContThresh).*ClassBidx(BadSet3));
    tpCx_PC(pp)=sum((ProdContNorm(pp,BadSet3)<=ProdContThresh).*ClassCidx(BadSet3));
    
    fpAy_PC(pp)=sum((ProdContNorm(pp,GoodSet)<=ProdContThresh).*ClassAidx(GoodSet));
    fpBz_PC(pp)=sum((ProdContNorm(pp,GoodSet)<=ProdContThresh).*ClassBidx(GoodSet));
    fpCx_PC(pp)=sum((ProdContNorm(pp,GoodSet)<=ProdContThresh).*ClassCidx(GoodSet));

    fnAy_PC(pp)=sum((ProdContNorm(pp,BadSet3)>ProdContThresh).*ClassAidx(BadSet3));
    fnBz_PC(pp)=sum((ProdContNorm(pp,BadSet3)>ProdContThresh).*ClassBidx(BadSet3));
    fnCx_PC(pp)=sum((ProdContNorm(pp,BadSet3)>ProdContThresh).*ClassCidx(BadSet3));

    tnAy_PC(pp)=sum((ProdContNorm(pp,GoodSet)>ProdContThresh).*ClassAidx(GoodSet));
    tnBz_PC(pp)=sum((ProdContNorm(pp,GoodSet)>ProdContThresh).*ClassBidx(GoodSet));
    tnCx_PC(pp)=sum((ProdContNorm(pp,GoodSet)>ProdContThresh).*ClassCidx(GoodSet));
    

    
            % For the exponential machines:
    PCCorrMaintAy(pp)=(CostCheckAy+CostFixAy)*tpAy_PC(pp);
    PCContFailAy(pp)= fnAy_PC(pp)*sum(CostsAy.*mean(aY,2));
    PCUnnecessaryMaintAy(pp)=(CostCheckAy)*fpAy_PC(pp);
            % For the weibull-1 machines:
    PCCorrMaintBz(pp)=(CostCheckBz+CostFixBz)*tpBz_PC(pp);
    PCContFailBz(pp)= fnBz_PC(pp)*sum(CostsBz.*mean(aZ,2));
    PCUnnecessaryMaintBz(pp)=(CostCheckBz)*fpBz_PC(pp);
            % For the weibull-2 machines:
    PCCorrMaintCx(pp)=(CostCheckCx+CostFixCx)*tpCx_PC(pp);
    PCContFailCx(pp)= fnCx_PC(pp)*sum(CostsCx.*mean(aX,2));
    PCUnnecessaryMaintCx(pp)=(CostCheckCx)*fpCx_PC(pp);

   
    %PostAlg_RISK_scenario = costTP*tpr + costFN*fnr + costFP*fpr;
    PCPostAlg_RISK_scenario(pp) = (PCCorrMaintAy(pp) + PCCorrMaintBz(pp) + PCCorrMaintCx(pp)) + ...
                             (PCContFailAy(pp) + PCContFailBz(pp) + PCContFailCx(pp)) + ...
                             (PCUnnecessaryMaintAy(pp) + PCUnnecessaryMaintBz(pp) + PCUnnecessaryMaintCx(pp));
    PCNoAlg_RISK_scenario = WeightedSystemRisk(pp); % inherent system risk
    PCValue_Alg(pp) = PCNoAlg_RISK_scenario-PCPostAlg_RISK_scenario(pp);

    % Focusing on false negatives, false positives of an applied algorithm
    PCAlgorithmsRisk(pp) =(PCContFailAy(pp) + PCContFailBz(pp) + PCContFailCx(pp)) + ...
                             (PCUnnecessaryMaintAy(pp) + PCUnnecessaryMaintBz(pp) + PCUnnecessaryMaintCx(pp));
end
% PQCI=ProdContNorm;

plot(1:pp,ProdContNorm(1:pp,:))
yline(ProdContThresh);
xlabel('Time (s)');
title({
    ['PQCI Values of Each Machine,' ] 
    ['Against the Threshold'] 
    });

legend(NodeNames{:},'Threshold','location','eastoutside');
hold off;

confVec=[tpAy_PC, tpBz_PC, tpCx_PC, ...
    fpAy_PC, fpBz_PC, fpCx_PC,...
    fnAy_PC, fnBz_PC, fnCx_PC,...
    tnAy_PC, tnBz_PC, tnCx_PC];

figure(33);
hold on;
plot(1:pp,confVec(1:pp,:))
legend('tp A','tp B','tp C',...
    'fp A', 'fp B','fp C',...
    'fn A', 'fn B', 'fn C',...
    'tn A', 'tn B','tn C', 'location','eastoutside')
hold off;

RiskValMat=[WeightedSystemRisk', PCPostAlg_RISK_scenario, PCValue_Alg, PCAlgorithmsRisk];
figure(34);
hold on;
plot(1:pp,RiskValMat(1:pp,:));
legend('Inherent System Risk','Post Algorithm Application','Value of Applied Algorithm',...
    'Algorithmic Risk Perspective (fp, fn)','location','eastoutside');
hold off;

%% Part 7: Risk Eval for Online Algorithm EPPCI (Variable: AddBadValueNorm)

PartQuality=Exp3.Output.PartQuality;
PartPathi=Exp3.Output.PartPathi;
PartPath=Exp3.Output.PartPath;

LagStore = 10;
NumPaths = length(UniquePaths); 
PathQuality = nan(LagStore,NumPaths); 
PathPerformance  = nan(NumParts3,NumNodes);
AddBadValueRange=PartQualityUB - PartQualityLB;
PathsMap = zeros(length(UniquePaths),NumNodes); 
for upi = 1:length(UniquePaths)
    PathsMap(upi,UniquePaths{upi}) = 1; 
end

Win = max(10,floor(NumParts3/100)); %Win=10 in the example with NumParts=500
PWN = nan(Win,NumNodes);%Part Quality With Node

AddBadValue=nan(NumParts3,NumNodes);
AddBadValueNorm=nan(NumParts3,NumNodes);
BadValueLikelihoodThresh= -0.1;

tpABV=nan(NumParts3,1);
tnABV=nan(NumParts3,1);
fpABV=nan(NumParts3,1);
fnABV=nan(NumParts3,1);

fprABV=nan(NumParts3,1);
tprABV=nan(NumParts3,1);
fnrABV=nan(NumParts3,1);

tpAy_ABV=nan(NumParts3,1);
tpBz_ABV=nan(NumParts3,1);
tpCx_ABV=nan(NumParts3,1);

fpAy_ABV=nan(NumParts3,1);
fpBz_ABV=nan(NumParts3,1);
fpCx_ABV=nan(NumParts3,1);

fnAy_ABV=nan(NumParts3,1);
fnBz_ABV=nan(NumParts3,1);
fnCx_ABV=nan(NumParts3,1);

tnAy_ABV=nan(NumParts3,1);
tnBz_ABV=nan(NumParts3,1);
tnCx_ABV=nan(NumParts3,1);


NodeNames={};
for Ni = 1:NumNodes  
    NodeNames{Ni} = sprintf('Node #%i',Ni);
end

GoodSet=1:10;
GoodSet(BadSet3)=[];

% Get the following from domain-level experts
    % For the exponential machines:
CostCheckAy=20;
CostFixAy=100;
ABVCorrMaintAy=nan(NumParts3,1);
ABVContFailAy= nan(NumParts3,1);
ABVUnnecessaryMaintAy=nan(NumParts3,1);
    % For the weibull-1 machines:
CostCheckBz=10;
CostFixBz=50;
ABVCorrMaintBz=nan(NumParts3,1);
ABVContFailBz= nan(NumParts3,1);
ABVUnnecessaryMaintBz=nan(NumParts3,1);
    % For the weibull-2 machines:
CostCheckCx=30;
CostFixCx=30;
ABVCorrMaintCx=nan(NumParts3,1);
ABVContFailCx= nan(NumParts3,1);
ABVUnnecessaryMaintCx=nan(NumParts3,1);


ABVPostAlg_RISK_scenario=nan(NumParts3,1);
ABVValue_Alg=nan(NumParts3,1);
ABVAlgorithmsRisk=nan(NumParts3,1);


figure(35);
hold on;
for pp = 1:NumParts3
    %Log Part Quality By the Part's Unique Path
    PathQuality(:,PartPathi(pp)) = ...
            [PathQuality(2:end,PartPathi(pp)); PartQuality(pp)];
    %Max Performance the Node Produces
    PathPerfPreprocess=diag(nanmean(PathQuality))*PathsMap;
    ind=PathPerfPreprocess==0;
    PathPerfPreprocess(ind)=nan;
    %Max Part Quality each machine produces when it is in use (in paths):
    PathPerformance(pp,:) = max(PathPerfPreprocess);
    for Ni = 1:NumNodes
        %If Part Passed this Node (machining stage)
        if any(PartPath{pp}==Ni)
            PWN(:,Ni) = [PWN(2:end,Ni); PartQuality(pp)]; %moving window, row is the window of (10) products and column is each machine/node
        end
    end

    AddBadValue(pp,:)=(PathPerformance(pp,:)+nanmean(PWN))/2;
    AddBadValueOffset=AddBadValue(pp,:)-mean(AddBadValue(pp,:));
    AddBadValueNorm(pp,:)=AddBadValueOffset/abs(AddBadValueRange);

    tpABV(pp)=sum(AddBadValueNorm(pp,BadSet3)<=BadValueLikelihoodThresh);
    fpABV(pp)=sum(AddBadValueNorm(pp,GoodSet)<=BadValueLikelihoodThresh);
    fnABV(pp)=sum(AddBadValueNorm(pp,BadSet3)>BadValueLikelihoodThresh);
    tnABV(pp)=sum(AddBadValueNorm(pp,GoodSet)>BadValueLikelihoodThresh);

    fprABV(pp)=fpABV(pp)/(tnABV(pp)+fpABV(pp));
    tprABV(pp)=tpABV(pp)/(tpABV(pp)+fnABV(pp));
    fnrABV(pp)=1-tprABV(pp);
    
    
    tpAy_ABV(pp)=sum((AddBadValueNorm(pp,BadSet3)<=BadValueLikelihoodThresh).*ClassAidx(BadSet3));
    tpBz_ABV(pp)=sum((AddBadValueNorm(pp,BadSet3)<=BadValueLikelihoodThresh).*ClassBidx(BadSet3));
    tpCx_ABV(pp)=sum((AddBadValueNorm(pp,BadSet3)<=BadValueLikelihoodThresh).*ClassCidx(BadSet3));
    
    fpAy_ABV(pp)=sum((AddBadValueNorm(pp,GoodSet)<=BadValueLikelihoodThresh).*ClassAidx(GoodSet));
    fpBz_ABV(pp)=sum((AddBadValueNorm(pp,GoodSet)<=BadValueLikelihoodThresh).*ClassBidx(GoodSet));
    fpCx_ABV(pp)=sum((AddBadValueNorm(pp,GoodSet)<=BadValueLikelihoodThresh).*ClassCidx(GoodSet));
  
    fnAy_ABV(pp)=sum((AddBadValueNorm(pp,BadSet3)>BadValueLikelihoodThresh).*ClassAidx(BadSet3));
    fnBz_ABV(pp)=sum((AddBadValueNorm(pp,BadSet3)>BadValueLikelihoodThresh).*ClassBidx(BadSet3));
    fnCx_ABV(pp)=sum((AddBadValueNorm(pp,BadSet3)>BadValueLikelihoodThresh).*ClassCidx(BadSet3));

    tnAy_ABV(pp)=sum((AddBadValueNorm(pp,GoodSet)>BadValueLikelihoodThresh).*ClassAidx(GoodSet));
    tnBz_ABV(pp)=sum((AddBadValueNorm(pp,GoodSet)>BadValueLikelihoodThresh).*ClassBidx(GoodSet));
    tnCx_ABV(pp)=sum((AddBadValueNorm(pp,GoodSet)>BadValueLikelihoodThresh).*ClassCidx(GoodSet));

    
            % For the exponential machines:
    ABVCorrMaintAy(pp)=(CostCheckAy+CostFixAy)*tpAy_ABV(pp);
    ABVContFailAy(pp)= fnAy_ABV(pp)*sum(CostsAy.*mean(aY,2));
    ABVUnnecessaryMaintAy(pp)=(CostCheckAy)*fpAy_ABV(pp);

                % For the weibull-1 machines:
    ABVCorrMaintBz(pp)=(CostCheckBz+CostFixBz)*tpBz_ABV(pp);
    ABVContFailBz(pp)= fnBz_ABV(pp)*sum(CostsBz.*mean(aZ,2));
    ABVUnnecessaryMaintBz(pp)=(CostCheckBz)*fpBz_ABV(pp);
    
                % For the weibull-2 machines:
    ABVCorrMaintCx(pp)=(CostCheckCx+CostFixCx)*tpCx_ABV(pp);
    ABVContFailCx(pp)= fnCx_ABV(pp)*sum(CostsCx.*mean(aX,2));
    ABVUnnecessaryMaintCx(pp)=(CostCheckCx)*fpCx_ABV(pp);

    %PostAlg_RISK_scenario = costTP*tpr + costFN*fnr + costFP*fpr;
    ABVPostAlg_RISK_scenario(pp) = (ABVCorrMaintAy(pp) + ABVCorrMaintBz(pp) + ABVCorrMaintCx(pp)) + ...
                             (ABVContFailAy(pp) + ABVContFailBz(pp) + ABVContFailCx(pp)) + ...
                             (ABVUnnecessaryMaintAy(pp) + ABVUnnecessaryMaintBz(pp) + ABVUnnecessaryMaintCx(pp));
    ABVNoAlg_RISK_scenario = WeightedSystemRisk(pp); % inherent system risk
    ABVValue_Alg(pp) = ABVNoAlg_RISK_scenario-ABVPostAlg_RISK_scenario(pp);

    % Focusing on false negatives, false positives of an applied algorithm
    ABVAlgorithmsRisk(pp) =(ABVContFailAy(pp) + ABVContFailBz(pp) + + ABVContFailCx(pp)) + ...
                             (ABVUnnecessaryMaintAy(pp) + ABVUnnecessaryMaintBz(pp) + ABVUnnecessaryMaintCx(pp));
end 
% EPPCI=AddBadValueNorm;

plot(1:pp,AddBadValueNorm(1:pp,:))
yline(BadValueLikelihoodThresh);
xlabel('Time (s)');
title({
    ['EPPCI Values of Each Machine,' ] 
    ['Against the Threshold'] 
    });

legend(NodeNames{:},'Threshold','location','eastoutside');
hold off;

confVecABV=[tpAy_ABV, tpBz_ABV, tpCx_ABV,...
    fpAy_ABV, fpBz_ABV, fpCx_ABV,...
    fnAy_ABV,fnBz_ABV, fnCx_ABV,...
    tnAy_ABV, tnBz_ABV, tnCx_ABV];

figure(36);
hold on;
plot(1:pp,confVecABV(1:pp,:))
legend('tp A','tp B','tp C',...
    'fp A', 'fp B','fp C',...
    'fn A', 'fn B','fn C', ...
    'tn A', 'tn B', 'tn C', 'location','eastoutside')
hold off;

RiskValMatABV=[WeightedSystemRisk', ABVPostAlg_RISK_scenario, ABVValue_Alg, ABVAlgorithmsRisk];
figure(37);
hold on;
plot(1:pp,RiskValMatABV(1:pp,:));
legend('Inherent System Risk','Post Algorithm Application',...
    'Value of Applied EPPCI Algorithm','Algorithmic Risk Perspective (fp, fn)','location','eastoutside');
hold off;



