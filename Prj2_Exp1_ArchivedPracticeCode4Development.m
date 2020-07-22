
%% Archived pratice code for different parts that were in development in this Project's Experiment


 %% =============== Part X: No Fail Time Probability Distributions  ================
% This is the older model, MfgAnalyzerV1, with a constant 'DegradeStart'


% Define input parameters 
BadSet=[2 5];
BadLevel = 1;
DegradeStart=10; 
NumParts=120;  
GoodLim=0;

% Manufacturing System Simulator
[Exp]=MfgAnalyzerV1(UniquePaths, NumParts, ExpNum, ...
    GoodLim, BadSet, BadLevel, DegradeStart)

% Generate Visualizations for the Simulator
MakePlots(Exp.Output,1);

%%  Some testing done to create MfgAnalyzerV2 & MakePlots2  
BadShift = ProdVal3*ones(length(BadSet3)+1,NumNodes);
for i=1:length(BadSet3)
    BadShift(i+1,BadSet3(1:i)) = BadLevel3;
end
StartValue=0; 


%2700 parts going through in 45 minutes
for pp = 1:NumParts3
    PartPath{pp} = UniquePaths{randi(length(UniquePaths))}; 
    if pp <= DegradeStart3(1)
        Shift=BadShift(1,:);
        PartQuality(pp) = StartValue + sum(Shift(PartPath{pp}))-rand*.01;
    else
        ShiftFinder=find(DegradeStart3<pp);
        ShiftFinder=ShiftFinder(end);
        Shift=BadShift(ShiftFinder+1,:);
        PartQuality(pp) = StartValue + sum(Shift(PartPath{pp}))-rand*.01;

    end
    PartOutTime(pp) = pp+rand*.2;
end           


[PartOutTime, ti] = sort(PartOutTime); %ti is the index of sorted parts
%Using the same sort index, here we make sure the part quality and path 
%values refer to the same product:
PartQuality = PartQuality(ti); 
PartPath = PartPath(ti);

%ID the Part Paths
PartPathi = nan(NumParts3,1); %An array of nan's, a value slot for each part
for pp = 1:NumParts3
    upi = 1;
    while ~isequal(PartPath{pp},UniquePaths{upi})
        upi= upi+1;
    end
    PartPathi(pp) = upi;
end %The ID indicates which one of the unique paths the part has gone through 

pause;

% Temporal Visualization Section:
winsize = 15;   
WinAvgs=meanfilt2(PartQuality,winsize);
figure(21);
plot(PartOutTime,WinAvgs,'-','linewidth',1.5) 
xlabel('Time (s)') %Cycles
ylabel('Part Quality') %Avg Quality
hold on
PathNames = {};
for upi = 1:length(UniquePaths) %for each of the unique paths 
    onP = PartPathi == upi; %binary, length= number of parts
    %Show Raw Data
    ph = plot(PartOutTime(onP),PartQuality(onP),'.'); 
    PathNames{upi} = sprintf('Path #%i Quality',upi);
end
legend('Average Quality',PathNames{:},'location','eastoutside')




%% This section checks the output tp, fp, tn, fn values for 
%  the two original Q-learning-related techniques
% Check Script

BadSet=BadSet3; ProdVal=ProdVal3; BadShift = ProdVal*ones(length(BadSet)+1,NumNodes);
GoodSet=1:length(BadShift(1,:));
GoodSet(BadSet)=[];
AddBadValue=Exp3.Output.AddBadValue; 
BadLevel=BadLevel3; 
StartValue=0;
PartQualityUB = StartValue+(length(UniquePaths{1})*ProdVal);
if length(BadSet)<length(UniquePaths{1})
    PartQualityLB = StartValue + (length(BadSet)*BadLevel) + ((length(UniquePaths{1}) - length(BadSet))*ProdVal);
else
    PartQualityLB = StartValue+length(UniquePaths{1})*BadLevel;
end
 AddBadValueRange=PartQualityUB - PartQualityLB;
AddBadValueOffset=AddBadValue-mean(AddBadValue);
AddBadValueNorm=AddBadValueOffset/abs(AddBadValueRange);
 
BadValueLikelihoodThresh= -0.1;

tpAddBadValue=sum(AddBadValueNorm(BadSet)<=BadValueLikelihoodThresh);
tnAddBadValue=sum(AddBadValueNorm(GoodSet)>BadValueLikelihoodThresh);
fpAddBadValue=sum(AddBadValueNorm(GoodSet)<=BadValueLikelihoodThresh);
fnAddBadValue=sum(AddBadValueNorm(BadSet)>BadValueLikelihoodThresh);
fprAddBadValue=fpAddBadValue/(tnAddBadValue+fpAddBadValue) 
tprAddBadValue=tpAddBadValue/(tpAddBadValue+fnAddBadValue)

%% THE FOLLOWING  IS PRACTICE CODE FOR DEVELOPING RISK METRICS:
% This is a slightly older version, I'm going to use it to see how commits/pull requests work for M-files

clear all;clc; close all;
load('ExperimentSet.mat');
% Risk experiment for Experiments 166:195
expVect=166:195; %For Experiments 88:117, 118:165
RiskExp=[Experiment([expVect]).Output];
numExp=length(RiskExp); %30
MaxBadNodes =3; %4 in Experiments 118:165

BadSet=nan(numExp,MaxBadNodes);
confMatGA=nan(numExp*2,2);
confMatFmM=nan(numExp*2,2);
confMatAddBadValue=nan(numExp*2,2);
confMatProdCont=nan(numExp*2,2);

for ii=1:numExp
    tempBadSet=RiskExp(ii).BadSet;
    BadSet(ii,1:length(tempBadSet))=RiskExp(ii).BadSet;
    confMatGA(2*ii-1:2*ii,1:2)=RiskExp(ii).confMatGA;
    confMatFmM(2*ii-1:2*ii,1:2)=RiskExp(ii).confMatFmM;
    confMatAddBadValue(2*ii-1:2*ii,1:2)=RiskExp(ii).confMatAddBadValue;
    confMatProdCont(2*ii-1:2*ii,1:2)=RiskExp(ii).confMatProdCont;
end
% Example:
%tpGA=confMatGA(1)
%fnGA=confMatGA(2)
%fpGA=confMatGA(3)
%tnGA=confMatGA(4)

% For GA:
% tp = correctly identified bad machines
tpvect=confMatProdCont(1:2:end,1); %vector of true positives for each scenario
% fp = false fault alarms; overuse of inspection/maintenance
fpvect=confMatProdCont(1:2:end,2); %vector of false positives for each scenario
% fn = missed fault alarms; bad product delivery to customer and overuse of
% bad machine
fnvect=confMatProdCont(2:2:end,1); 

numBadVect=sum(~isnan(BadSet(:,:))')';
costTP=2;
costFP=3;
costFN=4;
llhTP=tpvect; %likelihood/frequency of a TP
llhFP=fpvect; %likelihood/frequency of a FP
llhFN=fnvect; %likelihood/frequency of a FN


NoAlg_RISK_scenario=costFN*(llhTP+llhFN); %risk of each scenario without diagnostic algorithm
%(llhTP_llhFN) is the likelihood of the actual fault scenario, can be
%changed to reflect that 

PostAlg_RISK_scenario=costTP*llhTP+costFN*llhFN; %risk of each fault scenario

%Assumption is that costTP<costFN
Value_Alg=NoAlg_RISK_scenario-PostAlg_RISK_scenario;


RISK_algo_scenario=costFP*llhFP+costFN*llhFN; 
%risk of using algorithm for fault isolation on the scenario


probScenario=nan(size(numBadVect));

% Assume probabilities of different numbers of machine going bad:
probScenario(numBadVect==1)=0.02;
probScenario(numBadVect==2)=0.01;
probScenario(numBadVect==3)=0.001;

RISK_algorithm=sum(probScenario.*RISK_algo_scenario);


costConditions=nan(MaxBadNodes+1,1);
costConditions(1)=mean(RISK_algo_scenario(numBadVect==1));
costConditions(2)=mean(RISK_algo_scenario(numBadVect==2));
costConditions(3)=mean(RISK_algo_scenario(numBadVect==3));
costConditions(4)=0;
probConditions=nan(MaxBadNodes+1,1);
probConditions(1)=0.2;
probConditions(2)=0.1;
probConditions(3)=0.01;
probConditions(4)=0.69;
RISK_algorithm_conditions=sum(probConditions.*costConditions);




  
