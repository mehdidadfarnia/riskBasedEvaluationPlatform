function FunMakePlots_v02(Output, c)

%% Code to Make multiple Plots from different experiments in one run:
% for c = 1:SimNum 
%     MakePlots(Experiment(c).Output,c);
% end

%% Initialization & Inputs
k=c-1;

PQCI=Output.PQCI;
ProdContThresh=Output.ProdContThresh;
EPPCI=Output.EPPCI;
BadValueLikelihoodThresh=Output.BadValueLikelihoodThresh;
BadLevel=Output.BadLevel;
ProdVal=Output.ProdVal;
PartPath=Output.PartPath;
PartQuality=Output.PartQuality;
PartPathi=Output.PartPathi;
PartOutTime=Output.PartOutTime;
UniquePaths=Output.UniquePaths;
PathsMap=Output.PathsMap;
PathRank=Output.PathRank;
RankDeficiency=Output.RankDeficiency;
DeficiencyPerNumNodes=Output.DeficiencyPerNumNodes;
WinAvgs=Output.WinAvgs;
QWithN=Output.QWithN;
NodeNames=Output.NodeNames;
QDiff=Output.QDiff;
GoodLim=Output.GoodLim;
PWN=Output.PWN;
PWON=Output.PWON;
PathPerformance=Output.PathPerformance;
PathProb=Output.PathProb;
ProduceBadPart=Output.ProduceBadPart;
AddBadValue=Output.AddBadValue;
NodeConnectivity=Output.NodeConnectivity;
BadConnectivity=Output.BadConnectivity;
BadSetUse=Output.BadSetUse;
BadSetUseRatio=Output.BadSetUseRatio;
TotalBadConnectivityRatio=Output.TotalBadConnectivityRatio;
ConnectivityOverlap=Output.ConnectivityOverlap;
BadNodePresenceRatio=Output.BadNodePresenceRatio;
BadSet=Output.BadSet;
NodePresence=Output.NodePresence;
NumNodes=Output.NumNodes;
TotalPathsWithBadMachinesRatio=Output.TotalPathsWithBadMachinesRatio;
B2=Output.B2;
B3mean=Output.B3mean;
B3std=Output.B3std;
BadShift=Output.BadShift;




%% Figure 1: Temporal Production Quality Analysis
f1=figure((k*6)+1);
clf
plot(PartOutTime,WinAvgs,'-','linewidth',1.5) 
winsize=15;
% title({'Windowed Time Stats'
%         sprintf('^{(%i part lagging window)}',winsize)})
xlabel('Time (s)') %Cycles
ylabel('Part Quality') %Avg Quality
hold on
PathNames = {};
% Plot quality of each part produced through time. The for-loop does this
% by selecting all parts that have gone through each unique machining path
% and plotting their quality vs their production time. Each data plot
% % represents a part's quality and is designated by which unique path it
% went through to get produced. 
for upi = 1:length(UniquePaths) %for each of the unique paths 
    onP = PartPathi == upi; %binary, length= number of parts
    %Show Raw Data
    ph = plot(PartOutTime(onP),PartQuality(onP),'.'); 
    PathNames{upi} = sprintf('Path #%i Quality',upi);
end
legend('Average Quality',PathNames{:},'location','eastoutside')
%%%%

%% Figure 2: Machine/Node Quality Contribution Analysis 
%For each Node/Machine, we show a new average part quality value for 
%the last window of times (i.e. 10 times) that a node has been utilized,
%as this window shifts through time until all products are through the
%system. This can be an indicator of the node/machine quality.
f2=figure((k*6)+2);
subplot(2,3,1:2)
plot(PartOutTime,QWithN)
title('Average Path Quality by Node Usage')
xlabel('Time')
legend(NodeNames{:},'location','eastoutside')

%For each Node/Machine, we take the difference between the average node 
%quality (QWithN, plotted above) and the average part quality value for 
%the last window of times (i.e. 10 times) that a node has not been utilized,
%as this window shifts through time until all products are through the
%system. 
%This difference be an indicator of the node/machine contribution to quality.
subplot(2,3,4:5)
plot(PartOutTime,QDiff)
title('Realitive Conditional Node Contribution')
xlabel('Time')
ylabel('Contribution(%)')
legend(NodeNames{:},'location','eastoutside')

%Production Acceptability is defined here as the average part quality for
%the last window (i.e. 10 times) of times that a node has been utilized for
%a part, subtracted by the selected tolerance limit for a good part 
%(GoodLim). The difference in value can be an indicator of the 
%node/machine's contribution to product acceptability. 
subplot(2,3,3)
barh(nanmean(PWN)-GoodLim)
ylabel('Node #')
title('Production Acceptablity')
xlabel('Acceptablity(%)')

%Production Contribution is defined here as the average part quality for
%the last window (i.e. 10 times) of times that a node has been utilized for
%a part, subtracted by the average part quality for the last window of
%times that that node has not been utilized for a part. The difference in
%value can be an indicator of the node/machine's quality contribution. 
subplot(2,3,6)
barh(nanmean(PWN)-nanmean(PWON)) %Also variable ProdCont in code;
ylabel('Node #')
title('Production Contribution')
xlabel('Contribution(%)')


%% Figure 3: Part Path & Machine/Node Connectivity Analysis

f3=figure((k*6)+3);
subplot(2,3,[1 2])
plot(PartOutTime,PathPerformance)
legend(NodeNames{:},'location','eastoutside')
title({'Best Avg Part Quality That Each Node Produces'
    sprintf('Across All Paths, Irrespective of Use in Path')})

    %Path quality is the average of part qualities for parts that traverse the path (in a
    %window of time).
    %Each node/machine is either used in a path or it is not. If it is used, it
    %gets a score equal to that path's quailty. If not, it gets a 0. 
    %The maximum path quality for each node/machine is the maximum of these 
    %scores across all paths, and shows the best part quality performance that 
    %each node/machine produces. 
subplot(2,3,[4 5])
plot(PartOutTime,PathProb)
title('Probablity Node will Produce Good Part')
legend(NodeNames{:},'location','eastoutside')
    %Path acceptability (whether or not it is a good part) is the avg part acceptability for parts that
    %traverse the path (in a window of time). 
    %Each node/machine is either used in a path or it is not. If it is used, it
    %gets a score equal to that path's acceptability. If not, it gets a 0.
    %Max product acceptability for each node/machine is the maximum of these 
    %scores across all paths, and shows the best part acceptance that 
    %each node/machine produces. 

subplot(2,3,3)
barh(ProduceBadPart)
title({'Likelihood of Involved Node/Machine to'
    sprintf('Produce a Bad Part')})

subplot(2,3,6)
barh(AddBadValue)
title({'Likelihood of Involved Node/Machine to'
    sprintf('Add Bad Part Value')})


%% Figure 4: Q Learning-Related Methods, PQCI & EPPCI
f4=figure((k*6)+4);
subplot(1,3,1);
barh(PQCI)
hold on;
xline(ProdContThresh);
hold off;
ylabel('Machine #')
xlabel({'Normalized Production'
    sprintf('Contribution to Bad Parts')});
legend('Part Quality Contribution Indicator','Evaluation Threshold','location','northoutside');


subplot(1,3,3);
barh(EPPCI);
hold on;
xline(BadValueLikelihoodThresh);
hold off;
ylabel('Machine #')
xlabel({'Normalized Likelihood of'
    sprintf('Machine to Add Bad Value Part')});
legend('Estimated Part-Path Contibution Indicator','Evaluation Threshold','location','northoutside');


%% Figure 5: Node Connectivity
f5=figure((k*6)+5);surf(NodeConnectivity);view(2); colorbar;




%% Figure 6: Frequency of Machine/Node Usage

f6=figure((k*6)+6);
for jj=1:NumNodes
    hold on;
    if sum(jj==BadSet)==0
        hB=bar(jj,NodePresence(jj),'b');
    else
        hB=bar(jj,NodePresence(jj),'r');
    end
end
nColors=2;
hBLG = bar(nan(2,nColors)); %the bar object array for the legend
hBLG(1).FaceColor=[0 0 1]; %Blue
hBLG(2).FaceColor=[1 0 0]; %Red
hLG=legend(hBLG,'Good Machine','Bad Set', 'location', 'northeast');
title('Number of Paths That Machine/Node Appears In') %THIS is assuming that each path may include a machine/node only once (no repeat machines in a path)
xlabel('Machine/Node #')


%% Figure 7: Faulty Node Determination/Classification
f7=figure((k*6)+7);
clf
subplot(211)
bar(1:NumNodes,[B2,B3mean'])
legend('FmM','GA','location','eastoutside')
title('Constrained Solver Solutions: Visualization 1')
xlabel('Node #')
hold on
%Capture standard deviation of genetic algorithm
er = errorbar(1:NumNodes,B3mean',B3std,'om','LineStyle', 'none','DisplayName', 'GA std');
hold off

newBS=BadShift;
newBS(newBS==ProdVal)=1;
newBS(newBS==BadLevel)=0;

subplot(212)
plot(1:NumNodes,newBS(end,:),'b-',1:NumNodes,B2,'r--',1:NumNodes,B3mean,'k:')
legend('Real','FmM','GA','location','eastoutside')
title('Constrained Solver Solutions: Visualization 2')
xlabel('Node #')
hold on
%Capture standard deviation of genetic algorithm
er = errorbar(1:NumNodes,B3mean',B3std,'ok','LineStyle', 'none','DisplayName', 'GA std');
hold off


