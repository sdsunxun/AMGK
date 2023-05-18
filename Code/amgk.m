 %V3
 %20220908
 %归一化曲率，并自带。
function [cout,marked_img,cd] = amgk(varargin)

% Parse the input parameters
[I,High,Low,Gap_size,EP,Th] = parse_inputs(varargin{:});
if size(I,3)==3
    I=rgb2gray(I); % Transform RGB image to a Gray one. 
end
tic
% Detect edges
BW = edge(I,'canny',[Low,High]);
time_for_detecting_edge=toc;
%BW = edge(I,'canny_old',[0,0.35]);
tic
% Extract curves from the edge-image
[curve,curve_start,curve_end,curve_mode,curve_num,TJ,img1] = extract_curve(BW,Gap_size);  
time_for_extracting_curve=toc;
[sizex sizey] = size(I);
if size(curve{1})>0
    % Detect corners on the extracted edges
    tic
    [corner_out index cd2] = getcorner(curve,curve_mode,curve_start,curve_num,sizex,sizey,Th); 
    
    % Update the T-junctions
    [corner_final cd3] = Refine_TJunctions(corner_out,TJ,cd2,curve, curve_num, curve_start, curve_end, curve_mode,EP);
time_for_detecting_corner=toc;
   
%     img=I; % to show corners on the input image
img=~BW;%to show corners on the BW image
size(corner_final,1);
%     for i=1:size(corner_final,1)
%         img=mark(img,corner_final(i,1),corner_final(i,2),5);
%     end
% 
%     marked_img=img;
%     figure(); imshow(marked_img);
    cout = corner_final;
    cd = cd3;
%     figure; imshow(img); hold on; plot(cout(:,2),cout(:,1),'*r');
else
    cout = [];
    marked_img = [];
    cd = [];
end

here = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [corners index cd] = getcorner(curve,curve_mode,curve_start,curve_num,sizex,sizey,Th);

corners = [];
cor = []; % candidate corners
cd = [];


CLen = [4 6 8];
T_angle = acos((Th-1)*2)/pi()*180;

% CLen = [6 8 10];%%%for gray level
% Th=-0.92;
for i=1:curve_num;
    C3=[];
    C = [];
    SC=[]; % 有符号的曲率
    xs=curve{i}(:,1);
    ys=curve{i}(:,2);
    curveLen = size(xs,1);    
    W=10;
    if size(xs,1)>10 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!待修改！！！！！！！！！！
    if curve_mode(i,:)=='loop'
        xs1=[xs(curveLen-W+1:curveLen);xs;xs(1:W)];
        ys1=[ys(curveLen-W+1:curveLen);ys;ys(1:W)];
    else %expand the ends to gaussian window
        xs1=[ones(W,1)*2*xs(1)-xs(W+1:-1:2);xs;ones(W,1)*2*xs(curveLen)-xs(curveLen-1:-1:curveLen-W)];
        ys1=[ones(W,1)*2*ys(1)-ys(W+1:-1:2);ys;ones(W,1)*2*ys(curveLen)-ys(curveLen-1:-1:curveLen-W)];
    end   
    xs = xs1;
    ys = ys1;
    L = curveLen+2*W;   
    for j = 1:3
        chordLen = CLen(1,j);
        k=kcosine(xs,ys,chordLen,L,curve_mode(i,:),W);
        C3 =[C3,k];
    end

    c1 = C3(:,1);
    c2 = C3(:,2);
    c3 = C3(:,3);
    
    C = (c1+c2+c3)/3;
    
    SC=C;
    C=abs(C);
    
    L = curveLen;
    xs = xs(W+1:L+W);
    ys = ys(W+1:L+W);

    % Find curvature local maxima as corner candidates
    extremum=[];
    N=size(C,1);
    n=0;
    Search=1;
        
    for j=1:N-1
        if (C(j+1)-C(j))*Search>0
            n=n+1;
            extremum(n)=j;  % In extremum, odd points are minima and even points are maxima
            Search=-Search; % minima: when K starts to go up; maxima: when K starts to go down 
        end
    end
    if mod(size(extremum,2),2)==0 %to make odd number of extrema
        n=n+1;
        extremum(n)=N;
    end

    %%% accumulate candidate corners
    n = size(extremum,2);    
    for j = 1:n
        cor = [cor; curve{i}(extremum(j),:)];
    end   
    %%%
    
    n = size(extremum,2);
    flag = ones(size(extremum));
     
    % Compare each maxima with its contour average
    for j=2:2:n % if the maxima is less than local minima, remove it as flase corner
        if (C(extremum(j)) > Th) 
                flag(j)=0;
                [a,index1]=min(C(extremum(j):-1:extremum(j-1)));
                [a,index2]=min(C(extremum(j):extremum(j+1)));
                ROS=C(extremum(j)-index1+1:extremum(j)+index2-1);
                mean_k=mean(ROS);
                K_thre(j)=mean(ROS)+abs(mean(ROS))*0.05;     
                if C(extremum(j))<K_thre(j)
                    flag(j)=1;
                end
            end
    end
    extremum = extremum(2:2:n); % only maxima are corners, not minima
    flag = flag(2:2:n);
    extremum = extremum(find(flag==0));    
        
    % Check corner angle to remove false corners due to boundary noise and trivial details
    %fl = 0;
    %if fl
        flag=0;
        smoothed_curve=[xs,ys]; 
        while sum(flag==0)>0
            n=size(extremum,2);
            flag=ones(size(extremum)); 
            for j=1:n % second argument of curve_tangent function is always the position of the extrema in the first argument
                %which is array of points between two exterama
                if j==1 & j==n
                    ang=curve_tangent(smoothed_curve(1:L,:),extremum(j));
                elseif j==1 
                    ang=curve_tangent(smoothed_curve(1:extremum(j+1),:),extremum(j));
                elseif j==n
                    ang=curve_tangent(smoothed_curve(extremum(j-1):L,:),extremum(j)-extremum(j-1)+1);
                else
                    ang=curve_tangent(smoothed_curve(extremum(j-1):extremum(j+1),:),extremum(j)-extremum(j-1)+1);
                end     
                if ang>T_angle & ang<(360-T_angle) % if angle is between T_angle = 162 and (360-T_angle) = 198
                    flag(j)=0;  
                end
            end
             
            if size(extremum,2)==0
                extremum=[];            
            else
                extremum=extremum(find(flag~=0));             
            end
        end   
    extremum=extremum(find(extremum>0 & extremum<=curveLen)); % find corners which are not endpoints of the curve             
    index{i} = extremum';
    n = size(extremum,2);
    for j = 1:n
        corners = [corners; curve{i}(extremum(j),:)];
        cd = [cd; C(extremum(j))];
    end    
   
    fl = 1;
    if fl
     if curve_mode(i,:)=='loop'       
        if n>1
            compare_corner=corners-ones(size(corners,1),1)*curve_start(i,:);
            compare_corner=compare_corner.^2;
            compare_corner=compare_corner(:,1)+compare_corner(:,2);
            if min(compare_corner)>100      % Add end points far from detected corners, i.e. outside of 5 by 5 neighbor                                                
            if C(1) > Th
                 corners = [corners; curve_start(i,:)];
                 cd = [cd;5];
            end
                             
            end
        end
     end
    end  
    end
    here = 1;
end

here = 1;


        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cd = kcosine(xs,ys,chordLen,curveLen,curve_mode,W);


ax = [xs(W+1:curveLen-W)-xs(W-chordLen+1:curveLen-W-chordLen)];
ay = [ys(W+1:curveLen-W)-ys(W-chordLen+1:curveLen-W-chordLen)];
bx = [xs(W+1:curveLen-W)-xs(W+chordLen+1:curveLen-W+chordLen)];
by = [ys(W+1:curveLen-W)-ys(W+chordLen+1:curveLen-W+chordLen)];
Cd = (ax.*bx+ay.*by)./( (ax.^2+ay.^2).*(bx.^2+by.^2) ).^0.5;

Cd=(1+Cd)./2;

c_sign=-sign(ax.*by-ay.*bx);

Cd=Cd.*c_sign;

% figure()
% plot(Cd,'k-') 
% hold on

% figure()
% plot(Cd,'k:.') 
% hold on
for l = 1:3   %Mr zhong's method for weighted average in order to prevent the effect caused by noise or other ~
    for k = 2:curveLen-2*W-1
        Cd(k) = 0.25*Cd(k-1)+0.5*Cd(k)+0.25*Cd(k+1);
    end
    if curve_mode=='loop'
        Cd(1) = 0.25*Cd(curveLen-2*W)+0.5*Cd(1)+0.25*Cd(2);
        Cd(curveLen-2*W) = 0.25*Cd(curveLen-2*W-1)+0.5*Cd(curveLen-2*W)+0.25*Cd(1);
    end
end
% plot(Cd,'r.')
% plot(Cd,'r') 
% % axis([0,curveLen,-1,0.2])
here=1;

%%%%%%%%%%%55
function [xse yse] = enlarge(xs,ys,CL,curve_mode);
%CL = chord length
L = size(xs,1);
if curve_mode=='loop' % wrap around the curve by CL pixles at both ends
    xse = [xs(L-CL+1:L);xs;xs(1:CL)];
    yse = [ys(L-CL+1:L);ys;ys(1:CL)];
else % extend each line curve by CL pixels at both ends
    xse = [ones(CL,1)*2*xs(1)-xs(CL+1:-1:2);xs;ones(CL,1)*2*xs(L)-xs(L-1:-1:L-CL)];
    yse = [ones(CL,1)*2*ys(1)-ys(CL+1:-1:2);ys;ones(CL,1)*2*ys(L)-ys(L-1:-1:L-CL)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
function [xs ys W] = smoothing(x,y,L,curve_mode,gau,W);

if L>W
    if curve_mode=='loop' % wrap around the curve by W pixles at both ends
        x1 = [x(L-W+1:L);x;x(1:W)];
        y1 = [y(L-W+1:L);y;y(1:W)];
    else % extend each line curve by W pixels at both ends
        x1 = [ones(W,1)*2*x(1)-x(W+1:-1:2);x;ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
        y1 = [ones(W,1)*2*y(1)-y(W+1:-1:2);y;ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
    end
    
    xx=conv(x1,gau);
    xs=xx(2*W+1:L+2*W);
    yy=conv(y1,gau);
    ys=yy(2*W+1:L+2*W);    
else
    xs = [];
    ys = [];    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% extract curves from input edge-image
function [curve,curve_start,curve_end,curve_mode,cur_num,TJ,img]=extract_curve(BW,Gap_size)
%   Function to extract curves from binary edge map, if the endpoint of a
%   contour is nearly connected to another endpoint, fill the gap and continue
%   the extraction. The default gap size is 1 pixel
[L,W]=size(BW);
BW1=zeros(L+2*Gap_size,W+2*Gap_size);
BW_edge=zeros(L,W);
BW1(Gap_size+1:Gap_size+L,Gap_size+1:Gap_size+W)=BW;
[r,c]=find(BW1==1); %returns indices of non-zero elements
cur_num=0;

while size(r,1)>0 %when number of rows > 0
    point=[r(1),c(1)];
    cur=point;
    BW1(point(1),point(2))=0; %make the pixel black
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1); 
                               %find if any pixel around the current point is an edge pixel
    while size(I,1)>0 %if number of row > 0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [min_dist,index]=min(dist);
        p=point+[I(index),J(index)];
        point = p-Gap_size-1; % next is the current point
        cur=[cur;point]; %add point to curve 
        BW1(point(1),point(2))=0;%make the pixel black
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
                                %find if any pixel around the current point 
                                %is an edge pixel
    end
    
    % Extract edge towards another direction
    point=[r(1),c(1)];
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [min_dist,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[point;cur];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
        
    if size(cur,1)>(size(BW,1)+size(BW,2))/25 % for 512 by 512 image, choose curve if its length > 40
        cur_num=cur_num+1;                    % One can change this value to control the length of the extracted edges
        curve{cur_num}=cur-Gap_size;
    end
    [r,c]=find(BW1==1);
    
end

for i=1:cur_num
    curve_start(i,:)=curve{i}(1,:);
    curve_end(i,:)=curve{i}(size(curve{i},1),:);
    if (curve_start(i,1)-curve_end(i,1))^2+...
        (curve_start(i,2)-curve_end(i,2))^2<=25  %if curve's ends are within sqrt(32) pixels
        curve_mode(i,:)='loop';
    else
        curve_mode(i,:)='line';
    end
    BW_edge(curve{i}(:,1)+(curve{i}(:,2)-1)*L)=1;
end
%%%
if cur_num>0
    TJ = find_TJunctions(curve, cur_num, Gap_size+1); % if a contour goes just outsize of ends, i.e., outside of gapsize we note a T-junction there
else
    curve{1} = [];
    curve_start = [];
    curve_end = [];
    curve_mode = [];
    cur_num = [];
    TJ = [];    
end
%%%
img=~BW_edge;



% find T-junctions within (gap by gap) neighborhood, where gap = Gap_size +
% 1; edges were continued (see edge_extraction procedure) when ends are within (Gap_size by Gap_size)
function TJ = find_TJunctions(curve, cur_num, gap); % finds T-Junctions in planar curves
TJ = [];
for i = 1:cur_num
    cur = curve{i};
    szi = size(cur,1);
    for j = 1:cur_num
        if i ~= j
            temp_cur = curve{j};
            compare_send = temp_cur - ones(size(temp_cur, 1),1)* cur(1,:);
            compare_send = compare_send.^2;
            compare_send = compare_send(:,1)+compare_send(:,2);
            if min(compare_send)<=gap*gap       % Add curve strat-points as T-junctions using a (gap by gap) neighborhood
                TJ = [TJ; cur(1,:)];
            end
            
            compare_eend = temp_cur - ones(size(temp_cur, 1),1)* cur(szi,:);
            compare_eend = compare_eend.^2;
            compare_eend = compare_eend(:,1)+compare_eend(:,2);
            if min(compare_eend)<=gap*gap       % Add end-points T-junctions using a (gap by gap) neighborhood
                TJ = [TJ; cur(szi,:)];
            end
        end
    end
end
%%%

% Compare T-junctions with obtained corners and add T-junctions to corners
% which are far away (outside a 5 by 5 neighborhood) from detected corners
function [corner_final c3] = Refine_TJunctions(corner_out,TJ,c2,curve, curve_num, curve_start, curve_end, curve_mode,EP);
%corner_final = corner_out;
c3=c2;

%%%%% add end points
if EP
    corner_num = size(corner_out,1);
    for i=1:curve_num
            if size(curve{i},1)>0 & curve_mode(i,:)=='line'

                % Start point compare with detected corners
                compare_corner=corner_out-ones(size(corner_out,1),1)*curve_start(i,:);
                compare_corner=compare_corner.^2;
                compare_corner=compare_corner(:,1)+compare_corner(:,2);
                if min(compare_corner)>100       % Add end points far from detected corners 
                    corner_num=corner_num+1;
                    corner_out(corner_num,:)=curve_start(i,:);
                    c3 = [c3;8];
                end

                % End point compare with detected corners
                compare_corner=corner_out-ones(size(corner_out,1),1)*curve_end(i,:);
                compare_corner=compare_corner.^2;
                compare_corner=compare_corner(:,1)+compare_corner(:,2);
                if min(compare_corner)>100
                    corner_num=corner_num+1;
                    corner_out(corner_num,:)=curve_end(i,:);
                    c3 = [c3;9];
                end
            end
    end
end
%%%%%%%%%%%%%%%5

%%%%%Add T-junctions
corner_final = corner_out;
for i=1:size(TJ,1)
    % T-junctions compared with detected corners
    if size(corner_final)>0
        compare_corner=corner_final-ones(size(corner_final,1),1)*TJ(i,:);
        compare_corner=compare_corner.^2;
        compare_corner=compare_corner(:,1)+compare_corner(:,2);
        if min(compare_corner)>100       % Add end points far from detected corners, i.e. outside of 5 by 5 neighbor
            corner_final = [corner_final; TJ(i,:)];
            c3 = [c3;10];
        end
    else
        corner_final = [corner_final; TJ(i,:)];
        c3 = [c3;10];
    end
end


corner_last=[];
corner_final2=corner_final;
flag(1:size(corner_final2,1),1)=1;
for k=1:size(corner_final2,1)
    if(flag(k)==0)
        continue;
    end
    comp_corner=corner_final2(k,:);
    compare_corner=corner_final2-ones(size(corner_final2,1),1)*comp_corner;
    compare_corner=compare_corner.^2;
    compare_corner=compare_corner(:,1)+compare_corner(:,2);
    mini_index=find(compare_corner<32&flag==1&compare_corner>0);
    if(size(mini_index,1)>0)
        for m=1:size(mini_index)
            corner2=corner_final2(mini_index(m),:);
            valid_corner=[floor((comp_corner(1,1)+corner2(1,1))/2),floor((comp_corner(1,2)+corner2(1,2))/2),];
            flag(k)=0;
            flag(mini_index(m))=0;
            corner_last=[corner_last;valid_corner];
        end
    else
        corner_last=[corner_last;comp_corner];
    end
end
corner_final=corner_last;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show corners into the output images or into the edge-image
function img1=mark(img,x,y,w)
[M,N,C]=size(img);
img1=img;
if isa(img,'logical')
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<1);
    img1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:)=...
        img(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
else
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<128)*255;
    img1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:)=...
        img(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
function ang=curve_tangent(cur,center) % center is always the position of the corresponding extrema in cur

for i=1:2
    if i==1
        curve=cur(center:-1:1,:);
    else
        curve=cur(center:size(cur,1),:);
    end
    L=size(curve,1);
    
    if L>3
        if sum(curve(1,:)~=curve(L,:))~=0 % if not collinear
            M=ceil(L/2);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M,1);
            y2=curve(M,2);
            x3=curve(L,1);
            y3=curve(L,2);
        else
            M1=ceil(L/3);
            M2=ceil(2*L/3);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M1,1);
            y2=curve(M1,2);
            x3=curve(M2,1);
            y3=curve(M2,2);
        end
        
        if abs((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2))<1e-8  % straight line
            tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
        else
            % Fit a circle 
            x0 = 1/2*(-y1*x2^2+y3*x2^2-y3*y1^2-y3*x1^2-y2*y3^2+x3^2*y1+y2*y1^2-y2*x3^2-y2^2*y1+y2*x1^2+y3^2*y1+y2^2*y3)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            y0 = -1/2*(x1^2*x2-x1^2*x3+y1^2*x2-y1^2*x3+x1*x3^2-x1*x2^2-x3^2*x2-y3^2*x2+x3*y2^2+x1*y3^2-x1*y2^2+x3*x2^2)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            % R = (x0-x1)^2+(y0-y1)^2;

            radius_direction=angle(complex(x0-x1,y0-y1));
            if radius_direction<0
                radius_direction = 2*pi-abs(radius_direction);
            end
            
            adjacent_direction=angle(complex(x2-x1,y2-y1));
            
            if adjacent_direction<0
                adjacent_direction = 2*pi-abs(adjacent_direction);
            end
            
            tangent_direction=sign(sin(adjacent_direction-radius_direction))*pi/2+radius_direction;
            if tangent_direction<0
                tangent_direction = 2*pi-abs(tangent_direction);
            elseif tangent_direction>2*pi
                tangent_direction = tangent_direction-2*pi;
            end
        end
    
    else % very short line
        tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
    end
    direction(i)=tangent_direction*180/pi;
end
ang=abs(direction(1)-direction(2));
%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parses the inputs into input_parameters
function [I,Hg,Lo,Gap_size,EP,Th] = parse_inputs(varargin);

error(nargchk(0,5,nargin));
Para=[3,1,0.35,0,0.06]; %Default experience value;
if nargin>=2
    I=varargin{1};
    for i=2:nargin
        if size(varargin{i},1)>0
            Para(i-1)=varargin{i};
        end
    end
end

if nargin==1
    I=varargin{1};
end
    
if nargin==0 | size(I,1)==0
    [fname,dire]=uigetfile('*.bmp;*.jpg;*.gif','Open the image to be detected');
    I=imread([dire,fname]);
end

Gap_size = Para(1);
EP = Para(2);
Hg = Para(3); % high edge detection threshold
Lo = Para(4); % low edge detection threshold
Th=Para(5);
%%%%%%%%%%%%%%%%%%%%%%%%

function [G W] = makeGFilter(sig);

GaussianDieOff = .0001; 
pw = 1:100;

ssq = sig*sig;
W = max(find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff));
if isempty(W)
    W = 1;  
end
t = (-W:W);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq); 
G=gau/sum(gau);
