function [Experiment] = FunMfgSimulator_v02(UniquePaths, NumParts, ExpNum, GoodLim, BadSet, ProdVal, BadLevel, DegradeStart)

% Goal of this function is to analyze machine contribution to final quality
% of a product part, given that the products go through different paths of
% machines before becoming a finished product. 

%% Sections of this script:
 % 1. Input Data
 % 2. Rank of Machining Paths
 % 3. Temporal Production Quality Analysis
 % 4. Machine/Node Quality Contribution Analysis
 % 5. Part Path & Machine/Node Connectivity Analysis
 % 6. Node Connectivity
 % 7. Frequency of Machine/Node Usage
 % 8. Faulty Node Determination/Classification

%% Temporary Section: Potential Function Inputs
% Number of Parts = NumParts
    % Number of parts going through the total multi-stage manufacturing system
% Unique Paths = UniquePaths
    % Available Part Paths
% BadSet
    %The set of nodes/machines that degrade; defines BadShift in code
%BadLevel
    %How bad the nodes in BadSet get when they degrade; defines BadShift in
    %code
% ExpNum
    % Result save #
% Good quality limit = GoodLim
    % Tolerance limit to distinguish quality of a part
    % Think in terms of tolerances (part quality > 0 is good)
%DegradeStart


Input.ExpNum=ExpNum;
Input.UniquePaths=UniquePaths;
Input.NumParts=NumParts;
Input.GoodLim=GoodLim;
Input.BadSet=BadSet;
Input.BadLevel=BadLevel;
Input.DegradeStart=DegradeStart;
Input.timestamp=datetime;
Input.ProdVal=ProdVal;

%Discuss later:    
% BaseLine Error
% BadShift
%% Input Data
%This section defines:
%ExpNum = Result Folder Save #
%NumParts = Number of parts
%NumNodes = Number of machines/nodes
%UniquePaths = The different machining paths
%PartPath = Path for each part
%PartQuality = Quality for each part
%PartPathi = ID of the path from UniquePaths that the part goes through
%PartOutTime = Time stamp for when part reaches end of its path
%Set parameters 'BL' and 'BadShift' that determine part quality

%Assume a number of machines/nodes and the products have
%to go through a subset of the machines to reach the end and be a finalized product.


%Starting value of raw material going into production, could be an input
StartValue=0; 

%Find the total number N of machines/nodes in the multi-stage manufacturing system:
NumNodes = max(horzcat(UniquePaths{:})); 


%Badshift gives us different degrees of machine fails
BadShift = ProdVal*ones(length(BadSet)+1,NumNodes);
for i=1:length(BadSet)
    BadShift(i+1,BadSet(1:i)) = BadLevel;
end


%Here we define a machining path from the above list (of unique paths) 
%for each part.
for pp = 1:NumParts
    PartPath{pp} = UniquePaths{randi(length(UniquePaths))}; 
         %randi: specifies integer range between 1 and the # of unique
         %paths, uniformly distributed         
    
    % Defining part quality for each part. The bad shift starts after 
    % the part # indicated by DegradeStart, after which quality degrades at 
    % the stations selected above.
    % DegradeStart defines when degradation starts (by part #); i.e. a value of
    % 200 means that bad shift starts at the 200th product going through the
    % system
    if pp <= DegradeStart(1)
        Shift=BadShift(1,:);
        PartQuality(pp) = StartValue + sum(Shift(PartPath{pp}))-rand*.01;
    else
        ShiftFinder=find(DegradeStart<pp);
        ShiftFinder=ShiftFinder(end);
        Shift=BadShift(ShiftFinder+1,:);
        PartQuality(pp) = StartValue + sum(Shift(PartPath{pp}))-rand*.01;

    end
    

    %Parts output time: they come out pretty uniformly, per unit of time
    %given some small randomness
    PartOutTime(pp) = pp+rand*.2;
        %rand is a random number between 0 and 1, uniformly distributed
end           

%Sort Parts by Production Time
[PartOutTime, ti] = sort(PartOutTime); %ti is the index of sorted parts
%Using the same sort index, here we make sure the part quality and path 
%values refer to the same product:
PartQuality = PartQuality(ti); 
PartPath = PartPath(ti);

%ID the Part Paths
PartPathi = nan(NumParts,1); %An array of nan's, a value slot for each part
for pp = 1:NumParts
    upi = 1;
    while ~isequal(PartPath{pp},UniquePaths{upi})
        upi= upi+1;
    end
    PartPathi(pp) = upi;
end %The ID indicates which one of the unique paths the part has gone through 


winsize = 15;   
WinAvgs=FunMeanFilt_v02(PartQuality,winsize);

%Another way to define the input:

% % NumParts = 10000;
% % UniquePaths = {};
% % PartQuality = nan(NumParts,1);
% % PartOutTime = nan(NumParts,1);
% % PartPath = cell(NumParts,1);
% % PB = [0 1 0 1 0];
% % for pp = 1:NumParts
% %     PartPath{pp} = randi(5,1,3);
% %     damageAt = PB(PartPath{pp}) > rand(size(PartPath{pp}));
% %     if any(damageAt)&&pp>200
% %     PartQuality(pp) = 0;
% %     else
% %         PartQuality(pp) = 1;
% %     end
% %     PartOutTime(pp) = pp+rand*.2;
% %     %     %Store Unique Paths
% %     UniquePaths = HashPaths(UniquePaths,PartPath{pp});
% % end

Output.WinAvgs=WinAvgs;
Output.PartPath=PartPath;
Output.PartQuality=PartQuality;
Output.PartPathi=PartPathi;
Output.PartOutTime=PartOutTime;
Output.UniquePaths=UniquePaths;
Output.StartValue=StartValue;

%% Rank Of Machining Paths
%First, we create a matrix to represent the subset of machines/nodes (from
%the total available set of machines/nodes in the plant) traversed by each
%unique path.
%Properties of this matrix (rank) may tell us about the information we can
%derive about each node/machine during analysis. 

%Creates matrix to compare the number of unique machining paths vs the
%number of machines/nodes in the plant
PathsMap = zeros(length(UniquePaths),NumNodes); 

%With the for-loop we fill up the matrix with a logical binary 
%to show which machines the path goes through.
%Each row represents a unique path. A logical 1 is used for the subset of 
%nodes/machines traversed in the unique path. 
for upi = 1:length(UniquePaths)
    PathsMap(upi,UniquePaths{upi}) = 1; 
end

%Finds the rank of the matrix
PathRank = rank(PathsMap)
%Finds the rank deficiency - how far our matrix is from a full-rank matrix
RankDeficiency = PathRank - sum(any(PathsMap))

DeficiencyPerNumNodes=RankDeficiency/NumNodes;

Output.PathsMap=PathsMap;
Output.PathRank=PathRank;
Output.RankDeficiency=RankDeficiency;
Output.DeficiencyPerNumNodes=DeficiencyPerNumNodes;

%% Temporal Production Quality Analysis
%This section depicts the quality of products that go through the
%manufacturing system through time. This is compared against the average
%quality of products output from the system. 
%Recall: each part may go through a different sequence/path of
%machines/nodes.
%Think of production quality values in this section in terms of part 
%tolerances; good quality value limits are defined in the next section.

% % f1=figure(1);
% % clf;
% % %Time Window Averaging (for comparing part quality against average quality)
% % plot(PartOutTime,WinAvgs,'-','linewidth',1.5) 
% % title({'Windowed Time Stats'
% %         sprintf('^{(%i part lagging window)}',winsize)})
% % xlabel('Cycles')
% % ylabel('Avg Quality')
% % hold on
% % PathNames = {};
% % % Plot quality of each part produced through time. The for-loop does this
% % % by selecting all parts that have gone through each unique machining path
% % % and plotting their quality vs their production time. Each data plot
% % % represents a part's quality and is designated by which unique path it
% % % went through to get produced. 
% % for upi = 1:length(UniquePaths) %for each of the unique paths 
% %     onP = PartPathi == upi; %binary, length= number of parts
% %     %Show Raw Data
% %     ph = plot(PartOutTime(onP),PartQuality(onP),'.'); 
% %     PathNames{upi} = sprintf('Path #%i Quality',upi);
% % end
% % 
% % legend('Average Quality',PathNames{:},'location','eastoutside')

%Save Figure
%figuresdir = strcat('Results\',sprintf('Result%i',ExpNum),'\');
%filename=sprintf('Quality Per Product (with Path Infor) and Avg Quality, %i',ExpNum);
%saveas(f1,strcat(figuresdir,filename),'fig');


%% Machine/Node Quality Contribution Analysis (PQCI)
%Recall assumption that a number of machines/nodes and the parts have
%to go through a subset of the machines to reach the end and be a finalized product.
%This portion of the script depicts the indications of added quality from each
%node/machine to the final part. 


%Window to track Added Value From Node
Win = max(10,floor(NumParts/100)); %Win=10 in the example with NumParts=500

%Moving Window of:
PWN = nan(Win,NumNodes);%Part Quality With Node
PWON = nan(Win,NumNodes);%Part Qualty Without Node
%Average of Windows
QWithN = nan(NumParts,NumNodes);%Average Quality With Node
QWithoutN = nan(NumParts,NumNodes);%Average Quality Without Node

    
%Loop Parts Produced Through Time:
    %In the example of 10 machines where each part goes through a path
    %with only 4 machines, the for-loop indicates that for each part, only 4 of
    %the 10 PWN columns and 6 of 10 PWON get shifted (up) and all 10 of these shifts
    %are updated with the same part quality value. 
    %Then the average of each column (node) for the average quality (window)
    %for each part's row in QWN and QWON. 
for pp = 1:NumParts
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
    
%Average Quality With/Without Node:
    %After each part 'pp' updates a subset of the PWN columns and a subset of the
    %PWON columns with its part quality, the new averages of the PWN and
    %PWON columns give each Node/Machine a new average part quality value for 
    %the last window of times (i.e. 10 times) 
    %that a node has been utilized (for QWithN) or not utilized (for QWithoutN)
    QWithN(pp,:)   = nanmean(PWN);  %Average of each PWN column; each column characterizes a machine/node
    QWithoutN(pp,:)= nanmean(PWON); %Average of each PWON column; each column characterizes a machine/node
end

%Naming for each node to use in figure legends
NodeNames={};
for Ni = 1:NumNodes  
    NodeNames{Ni} = sprintf('Node #%i',Ni);
end

%For each Node/Machine, we show a new average part quality value for 
%the last window of times (i.e. 10 times) that a node has been utilized,
%as this window shifts through time until all products are through the
%system. This can be an indicator of the node/machine quality.
QDiff=QWithN -QWithoutN;

ProdCont= nanmean(PWN)-nanmean(PWON); %Contribution of each node to end-of-line product quality

%To figure out what the Product Contribution value bounds are, we look at
%maximum and minimum allowable values of product quality through the
%machining line:

PartQualityUB = StartValue+(length(UniquePaths{1})*ProdVal);
if length(BadSet)<length(UniquePaths{1})
    PartQualityLB = StartValue + (length(BadSet)*BadLevel) + ((length(UniquePaths{1}) - length(BadSet))*ProdVal);
else
    PartQualityLB = StartValue+length(UniquePaths{1})*BadLevel;
end

ProdContRange = PartQualityUB - PartQualityLB;
%ProdContLB = PartQualityLB - PartQualityUB; 
ProdContNorm = ProdCont/abs(ProdContRange); %range is from -1 to 1, where -1 is for node that produces only bad quality and 1 is the opposite)

%ProdContNorm2=(ProdContNorm./-2)+0.5; %range from 0 to 1, where 1 is for node that produces bad quality and 0 is the opposite.

ProdContThresh=-0.1; %%%%%%%%%%Check with michael/tim on this one!?
GoodSet=1:length(BadShift(1,:));
GoodSet(BadSet)=[];

tpProdCont=sum(ProdContNorm(BadSet)<=ProdContThresh);
tnProdCont=sum(ProdContNorm(GoodSet)>ProdContThresh);
fpProdCont=sum(ProdContNorm(GoodSet)<=ProdContThresh);
fnProdCont=sum(ProdContNorm(BadSet)>ProdContThresh);

Output.confMatProdCont = [tpProdCont, fpProdCont; fnProdCont, tnProdCont]; %Confusion Matrix

fprProdCont=fpProdCont/(tnProdCont+fpProdCont); 
tprProdCont=tpProdCont/(tpProdCont+fnProdCont);
accurProdCont=(tpProdCont+tnProdCont)/(tpProdCont+tnProdCont+fpProdCont+fnProdCont);
precisionProdCont=tpProdCont/(fpProdCont+tpProdCont);


Output.BadLevel=BadLevel;
Output.ProdVal=ProdVal;
Output.fprProdCont=fprProdCont;
Output.tprProdCont=tprProdCont;
Output.accurProdCont=accurProdCont;
Output.precisionProdCont=precisionProdCont;

Output.PQCI=ProdContNorm;
Output.ProdContThresh=ProdContThresh;
Output.QWithN=QWithN;
Output.NodeNames=NodeNames;
Output.QDiff=QDiff;
Output.GoodLim=GoodLim;
Output.PWN=PWN;
Output.PWON=PWON;
%Save Figure
%figuresdir = strcat('Results\',sprintf('Result%i',ExpNum),'\');
%filename=sprintf('Per Node Analysis %i',ExpNum);
%saveas(f2,strcat(figuresdir,filename),'fig');

%% Part Path & Machine/Node Connectivity Analysis (EPPCI)
%This portion of the script looks at connectivity between the nodes, their paths,
%and their contribution to part quality degradation.

LagStore = 10;
%Negative Voting - Counter Evidence

NumPaths = length(UniquePaths); 

%New parameters used to analyze the different machining paths to product parts:
PathQuality = nan(LagStore,NumPaths); 
    %For each part, matrix will update (column shift upwards) with values of part quality, 
    %for the column corresponding to the unique path that the part goes through
PathIsGood = nan(LagStore,NumPaths); 
    %For each part, matrix will update (column shift upwards) with a binary value of part 
    %acceptance (depending on if part quality > tolerance==GoodLim), 
    %for the column corresponding to the unique path that the part goes through
PathPerformance  = nan(NumParts,NumNodes);
    %For each part (row in this matrix), the average of PathQuality's columns (each
    %corresponding to a unique path that the part goes through - one of
    %these columns updated from the current part) is calculated and used to set up
    %a square, diagonalized matrix (dimension:NumPaths by Numpaths). This is then multiplied by the
    %PathsMap matrix (in the script section calculating the rank,
    %dimension: NumPaths by NumNodes) for a matrix with dimensions:NumPaths
    %by NumNodes.
    %The resulting matrix is similar to the PathsMap matrix with each
    %column representing the machines/nodes and each row representing a unique
    %path with a non-zero value for the subset of nodes/machines traversed in the
    %path. Instead of a logical binary 1 for this non-zero value, each row's (path's)values is
    %replaced by the average value of part quality that traversed that
    %particular path, averaged taken in the window captured by PathQuality.
    %This non-zero value can be seen as the path's quality, updated as
    %parts go through it. 
    %A row with the maximum value of each column (machine/node), no matter
    %which path they are in, tells us the max part quality performance that
    %the node/machine produces. This row is passed onto PathPerformance as
    %it's updated by each part. 
PathProb = nan(NumParts,NumNodes);
    %Similar to PathPerformance, but using PathIsGood instead of
    %PathQuality. At the end, a row with binary values of each column
    %(machine/node), no matter which path they are included in, tells us the max part
    %acceptability that the node/machine produces. This row is passed onto 
    %PathProb as it's updated by each part. 
   
% For-loop updates PathIsGood & PathQuality per product and uses them to
% calculate PathPerformance and PathProb.
for pp = 1:NumParts
    %Log Part Quality By the Part's Unique Path
    PathQuality(:,PartPathi(pp)) = ...
            [PathQuality(2:end,PartPathi(pp)); PartQuality(pp)];
    %Likewise, Log Good Parts By the Part's Unique Path
    PathIsGood(:,PartPathi(pp)) = ...
            [PathIsGood(2:end,PartPathi(pp)); PartQuality(pp) > GoodLim];
        
    %Max Performance the Node Produces
    PathPerfPreprocess=diag(nanmean(PathQuality))*PathsMap;
    ind=PathPerfPreprocess==0;
    PathPerfPreprocess(ind)=nan;
    %Max Part Quality each machine produces when it is in use (in paths):
    PathPerformance(pp,:) = max(PathPerfPreprocess);
    %Max Probablity the Node Produces a Good Part
    PathProb(pp,:) = max(diag(nanmean(PathIsGood))*PathsMap);
end 
 
%This metric ProduceBadPart is not good at finding a faulty node if all 
%the parts come out good but are of passable quality:
ProduceBadPart=(1-PathProb(pp,:)).*(nanmean(PWN<GoodLim));

%AddBadValue is used for EPPCI in the MSEC 2020 paper:
AddBadValue=(PathPerformance(pp,:)+nanmean(PWN))/2;


AddBadValueRange=PartQualityUB - PartQualityLB;
AddBadValueOffset=AddBadValue-mean(AddBadValue);

AddBadValueNorm=AddBadValueOffset/abs(AddBadValueRange);
BadValueLikelihoodThresh= -0.1;

tpAddBadValue=sum(AddBadValueNorm(BadSet)<=BadValueLikelihoodThresh);
tnAddBadValue=sum(AddBadValueNorm(GoodSet)>BadValueLikelihoodThresh);
fpAddBadValue=sum(AddBadValueNorm(GoodSet)<=BadValueLikelihoodThresh);
fnAddBadValue=sum(AddBadValueNorm(BadSet)>BadValueLikelihoodThresh);

Output.confMatAddBadValue = [tpAddBadValue, fpAddBadValue; ...
    fnAddBadValue, tnAddBadValue]; %Confusion Matrix

fprAddBadValue=fpAddBadValue/(tnAddBadValue+fpAddBadValue); 
tprAddBadValue=tpAddBadValue/(tpAddBadValue+fnAddBadValue);
accurAddBadValue=(tpAddBadValue+tnAddBadValue)/(tpAddBadValue+tnAddBadValue+fpAddBadValue+fnAddBadValue);
precisionAddBadValue=tpAddBadValue/(fpAddBadValue+tpAddBadValue);

Output.fprAddBadValue=fprAddBadValue;
Output.tprAddBadValue=tprAddBadValue;
Output.accurAddBadValue=accurAddBadValue;
Output.precisionAddBadValue=precisionAddBadValue;

Output.EPPCI=AddBadValueNorm;
Output.BadValueLikelihoodThresh=BadValueLikelihoodThresh;


Output.PathQuality=PathQuality;
Output.PathPerformance=PathPerformance;
Output.PathProb=PathProb;
Output.ProduceBadPart=ProduceBadPart;
Output.AddBadValue=AddBadValue;

%% Node Connectivity
%How often a 'Bad' Machine/Node is connected to other nodes (assuming
%selection of the unique path is random & uniform
NodeConnectivity=zeros(NumNodes,NumNodes); %Or zeros() instead of nan()? decide later
for jj = 1:length(UniquePaths)
   jjPath=UniquePaths{jj}; %The for-loop iterates through each path stored in UniquePaths
   nn=length(UniquePaths{jj}); %number of machines per path
   %nEdgesNum=nn*(nn-1)/2; %so for 4 nodes, we have 6 edges
   for qq = 1:(nn-1) %for 1 to 3
       for mm=qq:(nn-1) %for qq to 3
           NodeConnectivity(jjPath(qq),jjPath(mm+1))=NodeConnectivity(jjPath(qq),jjPath(mm+1))+1;
       end
   end
end
NodeConnectivity=NodeConnectivity+NodeConnectivity'; %order doesn't matter
NodeConnectivity = triu(NodeConnectivity,1);

BadConnectivityRow=NodeConnectivity(BadSet,:);
BadConnectivityCol=NodeConnectivity(:,BadSet);
BadConnectivity=BadConnectivityCol+BadConnectivityRow';
%Number of connections to other machines that each machine in the BadSet
%made. This should be equal to (# of unique paths it appeared) * (# of
%machines in each path - 1)
BadSetUse=sum(BadConnectivity); 
%Connections of each bad set machine as a ratio of total connections in the
%given unique paths
BadSetUseRatio=BadSetUse./(length(UniquePaths)*length(UniquePaths{1})*(length(UniquePaths{1})-1)/2);
%Connections of all bad set machine (subtracted by connections between bad set machines)
%as a ratio of total connections in the given unique paths
TotalBadConnectivityRatio=(sum(BadSetUse)-sum(sum(NodeConnectivity(BadSet,BadSet))))/(length(UniquePaths)*length(UniquePaths{1})*(length(UniquePaths{1})-1)/2); 
ConnectivityOverlap=sum(BadSetUseRatio)-TotalBadConnectivityRatio;

% % f4=figure(4);surf(NodeConnectivity);view(2); colorbar;


Output.NodeConnectivity=NodeConnectivity;
Output.BadConnectivity=BadConnectivity;
Output.BadSetUse=BadSetUse;
Output.BadSetUseRatio=BadSetUseRatio;
Output.TotalBadConnectivityRatio=TotalBadConnectivityRatio;
Output.ConnectivityOverlap=ConnectivityOverlap;
%Save Figure
%figuresdir = strcat('Results\',sprintf('Result%i',ExpNum),'\');
%filename=sprintf('Node Connectivity %i',ExpNum);
%saveas(f4,strcat(figuresdir,filename),'fig');

%% Frequency of Machine/Node Usage
NodePresence=nan(1,NumNodes);
for ii=1:NumNodes
    NodePresence(1,ii)=length(find(horzcat(UniquePaths{:})==ii));
end
cleanPaths=zeros(length(UniquePaths),1);
for kk = 1:length(UniquePaths)
    for rr=1:length(BadSet) %For all bad machines in bad set
        if sum(UniquePaths{kk}==BadSet(rr))==0 %If there are no bad machine in that path
            cleanPaths(kk)=cleanPaths(kk)+0; 
        else
            cleanPaths(kk)=cleanPaths(kk)+1;
        end
    end
end
TotalPathsWithBadMachinesRatio=(length(UniquePaths)-sum(cleanPaths==0))/length(UniquePaths);



BadNodePresenceRatio=NodePresence(BadSet)./NumNodes;


Output.BadNodePresenceRatio=BadNodePresenceRatio;
Output.BadSet=BadSet;
Output.NodePresence=NodePresence;
Output.NumNodes=NumNodes;
Output.TotalPathsWithBadMachinesRatio=TotalPathsWithBadMachinesRatio;

%% Faulty Node Determination/Classification
% This section is for comparing methods that determine (classify/predict) which
% nodes/machines are the ones that cause product defects (our original goal)

X = nanmean(PathQuality)';
loopCount = 5; %Number of times we run the algorithms below;
            %set to 5 if not testing

B3Store=nan(NumNodes, loopCount);
%Cost function
cost = @(b) [sum(abs((PathsMap*b) - (X))) ] ; 
%NEW: replaced (1-X) with (X)
%OLD COMMENT: (1-X) can be thought of the quality gap for each path, for all paths. It
%is the difference in ideal path quality (=1) vs the true path qualities
%stored in 'X'

b0 = randn(NumNodes,1)*.01; %Initialize
%  [B1,B1_unc] = regress(X,PathsMap);


% WHY CANT THIS BE DONE WITH FMINCON INSTEAD OF FMINIMAX
%Algorithms (constrained solvers) to solve the cost function:
%previously fminimiax
B2 = fmincon(cost,b0,[],[],[],[],zeros(size(b0)),ones(size(b0)));%zeros and ones are lower bounds and upper bounds, no other constraints
%Running the genetic algorithm a few times
for bb = 1:loopCount
    [B3,~,~,~,pop] = ga(@(x)cost(x'),NumNodes,[],[],[],[],zeros(NumNodes,1),ones(NumNodes,1));
    B3Store(:,bb)=B3';
end

%2 different visualizations to compare the solutions/classifications: a bar graph and a
%line plot. The line plot also includes the nodes/machines that cause part
%defects since we defined it with the 'Bl' and 'BadShift' parameters that
%were defined earlier in the script.


B3mean=nanmean(B3Store'); %GA mean per node
B3std=nanstd(B3Store'); %GA standard deve per node


Output.B3Store=B3Store;
Output.B3=B3;
Output.B2=B2;
Output.B3mean=B3mean;
Output.B3std=B3std;
Output.BadShift=BadShift;


%% Evaluating Predictions
thresh=0.6; 


tpGA=sum(B3mean(BadSet)<=thresh);
fpGA=sum(B3mean(GoodSet)<=thresh);
fnGA=sum(B3mean(BadSet)>thresh);
tnGA=sum(B3mean(GoodSet)>thresh);
Output.confMatGA = [tpGA, fpGA; fnGA, tnGA]; %Confusion Matrix

BS=BadShift(end,:);
errGA=immse(BS,B3mean); %MSE between actual and GA prediction; 
errGAbad=immse(BS(BadSet),B3mean(BadSet)); %MSE between actual and GA prediction of bad set
errGAgood=immse(BS(GoodSet),B3mean(GoodSet)); %MSE between actual and GA prediction of good set
fprGA=fpGA/(tnGA+fpGA); 
tprGA=tpGA/(tpGA+fnGA);
accurGA=(tpGA+tnGA)/(tpGA+tnGA+fpGA+fnGA);
precisionGA=tpGA/(fpGA+tpGA);
 
tpFmM=sum(B2(BadSet)<=thresh);
tnFmM=sum(B2(GoodSet)>thresh);
fpFmM=sum(B2(GoodSet)<=thresh);
fnFmM=sum(B2(BadSet)>thresh);
Output.confMatFmM = [tpFmM, fpFmM; fnFmM, tnFmM]; %Confusion Matrix

errFmM=immse(BS,B2'); %MSE between actual and minimax prediction; accuracy
errFmMbad=immse(BS(BadSet),B2(BadSet)');
errFmMgood=immse(BS(GoodSet),B2(GoodSet)');
fprFmM=fpFmM/(tnFmM+fpFmM); 
tprFmM=tpFmM/(tpFmM+fnFmM);
accurFmM=(tpFmM+tnFmM)/(tpFmM+tnFmM+fpFmM+fnFmM);
precisionFmM=tpFmM/(fpFmM+tpFmM);

Output.thresh=thresh;

Output.errGA=errGA;
Output.errGAbad=errGAbad;
Output.errGAgood=errGAgood;
Output.fprGA=fprGA;
Output.tprGA=tprGA;
Output.accurGA=accurGA;
Output.precisionGA=precisionGA;

Output.errFmM=errFmM;
Output.errFmMbad=errFmMbad;
Output.errFmMgood=errFmMgood;
Output.fprFmM=fprFmM;
Output.tprFmM=tprFmM;
Output.accurFmM=accurFmM;
Output.precisionFmM=precisionFmM;
%BadSetMeansRatioGA=mean(B3mean(BadSet))/mean(BadShift(BadSet)); %Similar to 'recall' or 'true positive rate'
%BadSetMeansRatioFmM=mean(B2(BadSet)')/mean(BadShift(BadSet)); %Similar to 'recall' or 'true positive rate'

Experiment.Input=Input;
Experiment.Output=Output;
%Save total data from script
%save(sprintf('Experiment.mat'),'Experiment')
%save(sprintf('Experiment %g .mat',ExpNum),'Input','Output')