clc
clear 
close all

% import a image file and convert to a .wav for laser display
image = imread('MH_test.jpg');

figure
imshow(image)

% convert to grascale 
I = rgb2gray(image);

figure
imshow(I)

% do a contor plot of the image
% use a single level to start with
figure;
[C, h] = imcontour(I,1);

% split out to individual contours
index = 1;
i = 1;
while index < size(C,2)
    contour(i).level = C(1,index); %#ok<SAGROW>
    
    % read in the cords
    contour(i).x = C(1,index+(1:C(2,index))); %#ok<SAGROW>
    contour(i).y = C(2,index+(1:C(2,index)))*-1; %#ok<SAGROW>
           
    i = i + 1;
    index = index + C(2,index) + 1;
end

max_size = max(size(image,1),size(image,2));
for n = 1:numel(contour)
    % move and scale between +- 1, maintain Aspect ratio
    % take size from original image
    contour(n).x = (contour(n).x / (0.5*max_size)) - 1;
    contour(n).y = (contour(n).y / (0.5*max_size)) + 1;

    % caculate the length, maybe also throw out by lenght, this is sort of
    % draw time i guess
    contour(n).length = sum(sqrt(diff(contour(n).x).^2 + diff(contour(n).y).^2));
    
    % start and end points of the contour
    contour(n).start = [contour(n).x(1),contour(n).y(1)];
    contour(n).end = [contour(n).x(2),contour(n).y(2)];
 end

% remove any that have zero lenght
contour([contour.length] == 0) = [];

% plot individual contors
figure 
hold all
xlim([-1,1])
ylim([-1,1])
for n = 1:numel(contour)
    plot(contour(n).x,contour(n).y)
end


% conect the contors using greedy algotithum
% somthing more fancy here could do a better job
% also this also only considers joining up from the start and end
% could be some improment to be had by backtracing, or joining part way

% start at with the longest
[ ~, index] = max([contour.length]);
uesed = false(numel(contour),1);
sorted = zeros(numel(contour),1);
for n = 1:numel(contour)
    
    uesed(index) = true;
    sorted(n) = index;
    
    % find the diatance from this contors end to the start of all contours    
    starts = cell2mat({contour.start}');
        
    distance = inf(numel(contour),1);
    distance(~uesed) = sqrt(sum((starts(~uesed) - contour(index).end).^2,2));
    
    [ ~, index] = min(distance);
end

% make sure we have used all the contors, but only onece each
if ~all(uesed) || numel(unique(sorted)) ~= numel(sorted)
    error('sorting cock up')
end

% convert into single X Y data set
x_data = zeros(numel([contour.x])+1,1);
y_data = zeros(numel([contour.x])+1,1);
index = 1;
for n = 1:numel(contour)
    contour_size = numel([contour(sorted(n)).x]);
    
    x_data(index+(0:contour_size-1)) = contour(sorted(n)).x;
    y_data(index+(0:contour_size-1)) = contour(sorted(n)).y;
    
    index = index + contour_size;
end
% make joined up
x_data(end+1) = x_data(1);
y_data(end+1) = y_data(1);


% plot joined up contors
figure 
plot(x_data,y_data)
xlim([-1,1])
ylim([-1,1])

figure
subplot(2,1,1)
plot(x_data)
subplot(2,1,2)
plot(y_data)

plot_time = 1/0.2;

%Sampling Frequency
Fs = 1/(plot_time/numel(x_data));


% do a low pass of 300hz
cutoff_freq = 30;

low_passed(1,:) = [0,0];
for i = 2:length(x_data)
    sample = [x_data(i),y_data(i)];
    
    rc = 1.0/(2*pi*cutoff_freq);
    
    alpha = (1/Fs)/((1/Fs)+rc);
    
    alpha = max(alpha,0);
    alpha = min(alpha,1);
    
    low_passed(i,:) = low_passed(i-1,:) + ((sample - low_passed(i-1,:)) * alpha);
end

figure
subplot(2,1,1)
plot(low_passed(:,1))
    
subplot(2,1,2)
plot(low_passed(:,2))  

figure
plot(low_passed(:,1),low_passed(:,2))
    
num_repeat = 50;

audiowrite('test.wav',repmat(myoutput,num_repeat,1),round(Fs))
    
    
    
