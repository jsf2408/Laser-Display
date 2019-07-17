clear;clc;close all;

warning('off','all')
warning

contourZero = zeros(1,11);
contourLength = 0:0.1:1;
contourLengthNegative = 1:-0.1:0;
contourOne = ones(1,11);

shapeX = [contourLength, contourOne, contourLengthNegative, contourZero];
shapeY = [contourZero, contourLength, contourOne, contourLengthNegative];

contours = struct('X',{},'Y',{},'m',{},'c',{},'backX',{},'backY',{});

contours(1).X = shapeX;
contours(1).Y = shapeY;

contours(2).X = shapeX+3;
contours(2).Y = shapeY;

contours(3).X = shapeX;
contours(3).Y = shapeY+3;

contours(4).X = shapeX+1.5;
contours(4).Y = shapeY+1.5;

figure
grid on
hold on

for i = 1:length(contours)
    m = [];
    c = [];
    for j = 1:length(contours(i).X)
        if j == length(contours(i).X)
            k = 1;
        else
            k = j;
        end
        x1 = contours(i).X(k);
        x2 = contours(i).X(k+1);
        y1 = contours(i).Y(k);
        y2 = contours(i).Y(k+1);
        coefficients = polyfit([x1, x2], [y1, y2], 1);
        m = [m, coefficients(1)];
        c = [c, coefficients(2)];
    end
    contours(i).m = m;
    contours(i).c = c;
end

allDist = [];
allStartX = [];
allStartY = [];
allEndX = [];
allEndY = [];
for g = 1:length(contours)
    minDist = [];
    endX = [];
    endY = [];
    startX = [];
    startY = [];
    for h = 1:length(contours)
        dist = [];
        x = [];
        y = [];
        x2 = [];
        y2 = [];
        for i = 1:length(contours(g).X)
            for j = 1:length(contours(h).X)
                d = sqrt((contours(h).X(j)-contours(g).X(i))^2+(contours(h).Y(j)-contours(g).Y(i))^2);
                coefficients = polyfit([contours(h).X(j), contours(g).X(i)], [contours(h).Y(j), contours(g).Y(i)], 1);
                m = coefficients(1);
                c = coefficients(2);
                save = 1;
%                 for k = 1:length(contours)
%                     for L = 1:length(contours(k).m)
%                         intercept = round((c-contours(k).c(L))/(contours(k).m(L)-m),3);
%                         v1 = round(contours(k).X(L),3);
%                         if L == length(contours(k).X)
%                             v2 = round(contours(k).X(1),3);
%                         else
%                             v2 = round(contours(k).X(L+1),3);
%                         end
%                         if ((v1<=intercept) && (intercept<=v2)) || ((v2<=intercept) && (intercept<=v1))
%                             if k ~= g && k ~= h
%                                 save = 0;
%                             end
%                         end
%                         %if k == 4 && ((g == 2 && h == 3)||(g == 3 && h == 2))
%                             fprintf('%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',k,L,g,h,intercept,v1,v2,save)
%                         %end
%                     end
%                 end
                if save == 1
                    dist = [dist, d];
                    x = [x,contours(g).X(i)];
                    y = [y,contours(g).Y(i)];
                    x2 = [x2,contours(h).X(j)];
                    y2 = [y2,contours(h).Y(j)];
                end
            end
        end
        [min_dist,idx] = min(dist(:));
        minDist = [minDist,min_dist];
        startX = [startX, x(idx)];
        startY = [startY, y(idx)];
        endX = [endX, x2(idx)];
        endY = [endY, y2(idx)];
    end
    allDist = [allDist; minDist];
    allStartX = [allStartX; startX];
    allStartY = [allStartY; startY];
    allEndX = [allEndX; endX];
    allEndY = [allEndY; endY];
end

allDist(~allDist) = inf;

[m,n] = size(allDist);
start = 1;
path = [];
path = [path,start];
for i = 1:m
    allDist(start,:) = inf;
    index = allDist(:,start);
    [val,start] = min(index((index > 0)));
    path = [path,start];
end
for i = 1:length(contours)
    plot(contours(i).X,contours(i).Y)
end

original = contours;

lineX = [];
lineY = [];
for i = 1:length(path)-1
    a = path(i);
    b = path(i+1);
    
    [~,sX] = find(contours(a).X==allStartX(a,b));
    [~,sY] = find(contours(a).Y==allStartY(a,b));
    
    index = intersect(sX,sY);
    index = index(1);

    contours(a).X = circshift(contours(a).X, -index);
    contours(a).Y = circshift(contours(a).Y, -index);    
end
for i = 1:length(path)-1
    a = path(i);
    b = path(i+1);
    
    [s,eX] = find(contours(b).X==allEndX(a,b));
    [d,eY] = find(contours(b).Y==allEndY(a,b));
    
    index1 = intersect(eX,eY);
    index1 = index1(1);
    
    contours(b).backX = (contours(b).X(index1:-1:1));
    contours(b).backY = (contours(b).Y(index1:-1:1));
end
figure
hold on
for i=1:length(path)-1
    plot(contours(path(i)).X,contours(path(i)).Y)
    plot(contours(path(i)).backX,contours(path(i)).backY)
end

figure
for i=1:length(path)-1
    lineX = [lineX,contours(path(i)).backX];
    lineY = [lineY,contours(path(i)).backY];
    lineX = [lineX,contours(path(i)).X];
    lineY = [lineY,contours(path(i)).Y];
end
plot(lineX,lineY)
