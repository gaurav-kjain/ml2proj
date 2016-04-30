BASE_DIR='d:/Study/MTech/Y2016/Project/sax_12/sax_12/';
%BASE_DIR = 'D:\Study\MTech\Y2016\Project\6\study\sax_10';

files = dir( fullfile(BASE_DIR,'*.dcm') );
files = {files.name}'; 
%Get height and width of images
fne = fullfile(BASE_DIR,files{1});
dicinfo= dicominfo(fne);

nrow=dicinfo.Height;
ncol=dicinfo.Width;
totframes=30;
pspace=dicinfo.PixelSpacing;

X = repmat(uint8(0), [nrow ncol 1 totframes]);

for i=1:numel(files)
    fname = fullfile(BASE_DIR,files{i});
    X(:,:,1,i) = medfilt2(dicomread(fname));
end

Ystd=zeros(nrow,ncol);
for irow=1:nrow
   for icol=1:ncol
       Ystd(irow,icol)=std2(X(irow,icol,1,:));
   end
end
stdmap=uint8(Ystd);

%stdmap=im2double(Ystd)/255;
%getting 80% of histogram value
%stdmap=double(stdmap)/double(max(stdmap(:)));
h=imhist(stdmap);
level = prctile(h,81);


%level  = graythresh(stdmap);
%tobinimg = X(:,:,1,12);
binimg = stdmap>level;
tobinimg = X(:,:,1,12);%systole image
segmtimg = uint8(tobinimg).*uint8(binimg);
%segmtimg=segmtimg/255;

%figure()
%imshow(stdmap);
%imshow(binimg);
%imshow(segmtimg);

hp=imhist(segmtimg);
hp(1)=0;

%figure,
T=otsuthresh(hp);
g=im2bw(tobinimg, T);

BW2 = imfill(binimg,'holes');
BW3 = bwareaopen(BW2, 200);


figure()
subplot(3, 3, 1);
imshow(uint8(tobinimg));
subplot(3, 3, 2);
imshow(segmtimg);
subplot(3, 3, 3);
imshow(binimg);
subplot(3, 3, 4);
imshow(g);
subplot(3, 3, 5);
imshow(BW2);
subplot(3, 3, 6);
imshow(BW3);

Ysys=zeros(nrow,ncol);
for irow=1:nrow
   for icol=1:ncol
       if(g(irow,icol)==1 && BW3(irow,icol)==1)
           Ysys(irow,icol)=127;%blood area
       end
       if(BW3(irow,icol)==1 && g(irow,icol)==0)
               Ysys(irow,icol)=255;%non blood area
       end
   end
end

figure()
imshow(uint8(Ysys));
bfin=bwperim(uint8(Ysys));
figure()
imshow(bfin);


%center of image
csysh=nrow/2;
csysw=ncol/2;
approxC=110;
%total range of ROI to consider
pixsysh=ceil( double(approxC) * double(pspace(1)));
pixsysw=ceil( double(approxC) * double(pspace(2)));

if(rem(pixsysh,2) ~= 0)
    pixsysh=pixsysh+1;
end

if(rem(pixsysw,2) ~= 0)
    pixsysw=pixsysw+1;
end


%Draw ROI
rectangle = int32([(csysw-pixsysw/2) (csysh-pixsysh/2) pixsysw pixsysh]);
RGB = insertShape(tobinimg, 'rectangle', rectangle, 'LineWidth', 5); 
%imshow(uint8(RGB));

ROIImgBin   = imcrop(binimg,rectangle);
ROIImg      = imcrop(tobinimg,rectangle);

for irow=1:pixsysh
   for icol=1:pixsysw
      ROIImg(irow,icol) =(ROIImgBin(irow,icol)*ROIImgBin(irow,icol))/(255*255);
   end
end

figure()
subplot(2, 2, 1);
imshow(stdmap);

subplot(2, 2, 2);
imshow(binimg);

subplot(2, 2, 3);
imshow(ROIImg);

subplot(2, 2, 4);
imshow(ROIImgBin);

    
%features extraction
[Gx,Gy] = imgradientxy(img);
figure()
subplot(3, 3, 1);
imshow(Gx);
subplot(3, 3, 2);
imshow(Gy);

[Gxx,Gxy] = imgradientxy(Gx);

subplot(3, 3, 3);
imshow(Gxx);
subplot(3, 3, 4);
imshow(Gxy);

%multiply binary image with 10th image

%element wise product
Lx=double(img).*double(Gx);
subplot(3, 3, 5);
imshow(Lx);

[Gmag, Gdir] = imgradient(Gx, Gy);
subplot(3, 3, 6);
imshow(uint8(Gmag));

img_gradx_max=imregionalmax(Gmag);
subplot(3, 3, 7);
imshow(img_gradx_max);

