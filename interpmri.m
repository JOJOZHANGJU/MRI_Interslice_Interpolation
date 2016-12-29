function[] = interpmri(imDir,interpnum,outDir)
% INTERPMRI can be used to intsert MR images between two images.
% dir is the directory that contains the MR images
% internum is the number of images to be insert between images.
if nargin==1
    interpnum = 2;
elseif nargin==2
    % output dir
    outDir = fullfile(imDir,'interp');
end

if ~exist(imDir,'dir')
    error('The MR image directory does not exist');
end

if ~exist(outDir,'dir')
    mkdir(outDir);
end

disp(['Interp MR images,input:',imDir,',output:',outDir,'.']);

% get image file list
imList = dir(fullfile(imDir,'*.bmp'));
imNum = length(imList);
imAll = cell(1,imNum);
if imNum==0
    disp('No *.bmp images');
    return;
end

disp(['input/output=',num2str(imNum),'/',num2str((imNum-1)*(1+interpnum)+1)]);

% load all images
for i=1:imNum
    fname = fullfile(imDir,imList(i).name);
    im = imread(fname);
    imAll{i} = im;
end

[IE,JE] = size(im);
prefix = getPrefix(imList(i).name);

% interp
num = 0;
X = [1 (interpnum+2)];
XI = 1:(interpnum+2);
for k=1:imNum-1
    curim = imAll{k};
    nxtim = imAll{k+1};
    Y = zeros(2,IE,JE);
    Y(1,:,:) = double(curim);
    %Y(1,:,:) = uint8(curim);
    Y(2,:,:) = double(nxtim);
    %Y(2,:,:) = uint8(nxtim);
    IMI = zeros(interpnum,IE,JE);
%     YI = zeros(2+interpnum,IE,JE);
%     for i=1:IE
%         yt = zeros(2,JE);
%         yt(:,:) = Y(:,i,:);
%         YI(:,i,:) = interp1(X,yt,XI);
%     end
    %YI = interp1(X,Y,XI);
    YI = interp1(X,Y,XI,'spline');
    IMI(:,:,:) = YI(2:end-1,:,:);
%     for i=1:IE
%         for j=1:JE
%             Y  = double([curim(i,j),nxtim(i,j)]);
% %             YI = interp1(X,Y,XI);
%             
% 
%             IMI(:,i,j) = YI(2:end-1);
%         end
%     end
    
    % save the images
    num = num+1;
    saveImage(curim,outDir,prefix,num);
    for i=1:interpnum
        num = num+1;
        im = zeros(IE,JE,'uint8');
        im(:,:) = IMI(i,:,:);
        saveImage(im,outDir,prefix,num);
    end
    disp(['instert between ',num2str(k),' and ',num2str(k+1)]);
end
num = num+1;
saveImage(nxtim,outDir,prefix,num);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save the images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[] = saveImage(im,dir,prefix,num)
    im = uint8(im);
    ind = [strrep(sprintf('%3d',num),' ','0'),'.tif'];
    fname = fullfile(dir,[prefix,'_',ind]);
    imwrite(im,fname,'TIFF');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the prefix of an image file, exclude the number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[prefix] = getPrefix(fn)
    prefix = '';
    ins1 = strfind(fn,'_');
    if isempty(fn)
        return;
    end
    
    ins2 = strfind(fn,'.');
    if ~isempty(ins2) && ins2(end)>ins1(end)
        num = fn(ins1(end):ins2(end)-1);
        if ~isempty(num2str(num))
            prefix = fn(1:ins1(end)-1);
        end
    end
end