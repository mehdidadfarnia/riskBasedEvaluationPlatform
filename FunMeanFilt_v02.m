function xf = FunMeanFilt_v02(x,win)
%Mean Filter Data
%  xf = meanfilt1(x,win)
%Mean filter data (x) useing a window size (win)
%Default window size = 10
%** LAGGING Window ** 
if nargin < 2
    win =10;
end
%size Data
[o, v] = size(x);

if o == 1
    x = x';
    [o, v] = size(x);
end

%Allocate
xf = nan(o,v);

for i = 1:o
    si = max(1,i-win);
    
    xf(i,:) = mean(x(si:i,:));
end

% 
% %Half Window
% hwin = ceil(win/2);
% 
% %Initalize Scrolling Sum
% rs = sum(x(1:hwin,:))*2+x(hwin+1,:);
% 
% for i = 1:o
%     %Upper Index
%     ui = i+hwin;
%     if ui > o
%         ui = 2*o - ui;
%     end
%     
%     %Lower index
%     li = i-hwin;
%     if li < 0
%         li = -li;
%     elseif li == 0
%         li = 1;
%     end
%     %Scrolling Sum
%     rs = rs + x(ui,:) - x(li,:);
%     
%     %Scrolling mean
%     xf(i,:) = rs/win;
% end
% 
%  % %Fix Missing Values
% xf(xf-xf~=0) = mean(x(x-x==0));