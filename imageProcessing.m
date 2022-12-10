clc
clear all
close all

% Read image
folder = "fig/trial5/";
files = dir(folder+"*.jpg");
for i = 1:length(files)
% % % %     figure()
    I = imread(folder+files(i).name);
% % % %     subplot(1,2,1)
% % % %     imshow(I)
    
    % Search for gradients
    Igray = rgb2gray(I);
    Igrad = imgradient(Igray);
    
%     max_grad = max(max(Igrad));
%     I2 = imbinarize(Igrad,max_grad/10);
%     subplot(1,2,2)
%     imshow(I2)
%     hold on
%     [y,x] = find(I2);
%     ellipse_t = fit_ellipse(x,y,gca);
% % % % %     subplot(1,2,2)
    [y,x] = find(Igray);
    a = (max(x)-min(x))/2;
    b = (max(y)-min(y))/2;
    xc = (max(x)+min(x))/2;
    xcs(i) = xc;
    yc = (max(y)+min(y))/2;
    ycs(i) = yc;
    I3 = imbinarize(Igrad);
% % % % %     imshow(I3)
% % % % %     hold on
    xs = [xc,max(x),min(x),xc,xc];
    ys = [yc,yc,yc,max(y),min(y)];
% % % % %     plot(xs,ys,'rx','MarkerSize',10)
%     plot(max(x),yc,'x')
%     plot(min(x),yc,'x')
%     plot(xc,max(y),'x')
%     plot(xc,min(y),'x')
end
