%%   User-fiendly code based on the following article:
%$£  Segmented Fluorescence Correlation Spectroscopy (FCS) on a commercial laser scanning microscope
%$£  Sci Rep 14, 17555 (2024). https://doi.org/10.1038/s41598-024-68317-7
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£ 
%$£                                                                     %$£
%$£                Elisa Longo, Silvia Scalisi, Luca Lanzanò            %$£
%$£      University of Catania - Department of Physics and Astronomy    %$£
%$£         Istituto Italiano di Tecnologia - Nanoscopy Department      %$£
%$£                      User-Friendly Version (30-07-2024)             %$£
%$£                                                                     %$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£

function [Iav, Hav]=SegmentedFCS_v0(filename)


% Processes FCS data acquired by SLOW line/raster scan or FAST line scan in a commercial microscope.
% Each line (SLOW) or Column (FAST) is divided in segments. Each segment is processed to produce ACF. 
% The segments are selected and the corresponding ACFs averaged.
%
% INPUT FILE:
% The input file must be .tif 
% X-segmentation:
% - single image (1 chanel or 2 channel) 
% NOTE: If the file has 2 channels, the script will calculate the cross-correlation
% Y-segmentation:
% - image or stack of images (e.g. scanning FCS data are normally exported
% as a stack of multiple images)
%  NOTE: the user has to specify if the data has 1 channel or 2 channels
%
%
% OUTPUT FILES:
% - ACF function exported as .txt 
% - Mask exported as .tif
% - Normalized image exported as .tif

answerfile = questdlg('Select segmented FCS mode:', ...
	'Options', ...
	'X-Segmentation (slow scan FCS)','Y-Segmentation (fast scan FCS)','X-Segmentation (slow scan FCS)');

switch answerfile
    case 'X-Segmentation (slow scan FCS)'                      
 % SLOW SCAN 1 channel or 2  Channels (cross)
%open image (Y lines data acquired as xt or xy) and reshape into a matrix D
%of size Y per T (=x pixels)
%open files ... 
[filename,pathname, filterindex] = uigetfile({'*.tif'},'Please select an image ');
filenamefull = [pathname, filename];   
filenamenoext = [pathname, filename(1:end-4)];
D0=ISFCS_readfiletif(filenamefull);
% A(isnan(A))=0;
% T=size(A,3);
Z=size(D0,3);
D=D0(1:1:end,:,1);
D2=D;
if Z> 1
    D2=D0(1:1:end,:,2);
end        
subseg=8;
pxdef=0.47 ;
Nch=1;

% set temporal resolution (Leica sp8 at 1Hz, is 0.0305 ms)
prompt = {'Pixel time (ms)'}; 
dlg_title = 'Settings'  ; 
num_lines = 1;
def = {num2str(0.0305)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
dtpixel=str2num(answer{1});
dt=dtpixel/1000; % dwell time in s
    
    case 'Y-Segmentation (fast scan FCS)'
% LINE SCAN 1 Channel 
%open image stack (N lines acquired as xt mode) and reshape into a matrix D of size X per T (=N lines)
[stackname,pathname, filterindex] = uigetfile({'*.tif'},'Select an image stack ', 'MultiSelect', 'on');
if iscell(stackname) % this is for multiple images, not yet modified for 2 ch
stacknamefull{1} = [pathname, stackname{1}];  
D=ISFCS_readfiletif(stacknamefull{1});
D=permute(D, [2 1]);
for i=2:length(stackname)
stacknamefull{i} = [pathname, stackname{i}];  
A1=ISFCS_readfiletif(stacknamefull{i});
A1=permute(A1, [2 1]);
D=cat(2,D,A1);
end
filenamenoext = [pathname, stackname{i}(1:end-4)];
X=size(D,1);
T=size(D,2);
D2=D;

else   % this is OK: for stacks
    
filenamefull = [pathname, stackname];
filenamenoext = [pathname, stackname(1:end-4)];
fileext=stackname(end-2:end);
if fileext=='tif'
A0=ISFCS_readfiletif(filenamefull);  
elseif fileext=='avi' % not used
    A=ISFCS_AVI2Matrix(filenamefull);
end
A0=permute(A0, [2 1 3]);


% set temporal resolution (Leica sp8 at 1Hz, is 0.0305 ms)
prompt = {'Line time (ms)', 'Number of channels (1 or 2)'}; 
dlg_title = 'Settings'  ; 
num_lines = 1;
def = {  num2str(0.55), num2str(1) };
answer = inputdlg(prompt,dlg_title,num_lines,def);
dtline=str2num(answer{1});
Nch=str2num(answer{2});
dt=dtline/1000;        

if Nch==2
A=A0(:,:,1:2:end);
A2=A0(:,:,2:2:end);
else
A=A0;    
A2=A0;
end
Nframes=size(A,3);
X=size(A,1);
YF=size(A,2);
T=YF*Nframes;
D=reshape(A,X,T,1);   
D2=reshape(A2,X,T,1);  
clear A A2;
end             
        
subseg=32;
pxdef=100 ;

end

Y=size(D,1);
T=size(D,2);  

% load background if available, for G0 correction
[filenamepatt,pathnamepatt, filterindex] = uigetfile({'*.mat'},'Load background file?');
if filenamepatt ~= 0
filenamefullpatt = [pathnamepatt, filenamepatt];   
load( [filenamefullpatt], 'Iav', 'Hav' ) ;
Ibkgd=Iav;
Gbkgdt=Hav;
else
    Ibkgd=0;
    Gbkgdt=0;
end   



tag='';
segmin=1;

maxseg=256;
Lseg=ceil(Y/maxseg);

Lseg=1;
nseg=floor(Y/Lseg);
 

% figure
% imagesc(imresize(double(D),0.125))

Lseg=1;
Nlog=120;

% figure
% plot(1:nseg,I)
% figure; imagesc(Inorm');
% figure; imagesc(I');

% selection of segments and calculation of final ACF

Thrg1=1.1;
Type='Low';
colmin=1;
minseguser=1;
maxseguser=nseg;
figcheck=1;
while figcheck==1 
clear    Hseg2 Htemp Hav
colmax=subseg;
Tseg=floor(T/subseg);
% data segmentation and ACF calculation (data can be Nlines * Npx in SLOW  or Npx * Nlines in FAST)
I=zeros(subseg,nseg) ;
Inormx=I;
Inormy=I;
for i=1:nseg
    for j=1:subseg
Aseg=D(  1+ Lseg*(i-1):Lseg*(i), 1+ Tseg*(j-1):Tseg*(j)  );
Aseg2=D2(  1+ Lseg*(i-1):Lseg*(i), 1+ Tseg*(j-1):Tseg*(j)  );
% I(i)=mean(Aseg,'all');
I(j,i)=mean(Aseg,'all');
[Hlog, taulog, Hlin, tau]=ISFCS_imgdata_X(Aseg,Aseg2,dt,0,Nlog,0) ;
Hseg2(j,i,:) = mean(Hlog,1);
    end
if max(I(:,i)) ~= min(I(:,i))    
Inormx(:,i)=(I(:,i)-min(I(:,i)))./(max(I(:,i))-min(I(:,i)));  % normalized at each line (0,1)    
else
Inormx(:,i)=0;    
end  
% Hseg(i,:) = squeeze(mean(Hseg2(i,:,:),2));

end
for j=1:subseg
if max(I(j,:)) ~= min(I(j,:))
    Inormy(j,:)=(I(j,:)-min(I(j,:)))./(max(I(j,:))-min(I(j,:))); 
else
    Inormy(j,:)=0;
end
end

switch answerfile
    case 'X-Segmentation (slow scan FCS)' 
Inorm=Inormx;
    case 'Y-Segmentation (fast scan FCS)' 
Inorm=Inormy;       
end    
    
switch answerfile
    case 'X-Segmentation (slow scan FCS)' 
[mask]=FCS_Threshold(Inorm',Thrg1, Type,colmin,colmax,minseguser,maxseguser,3); % masks matrix: maxseg per X!
maskseg=mask';
imagenorm=Inorm';
    case 'Y-Segmentation (fast scan FCS)' 
[mask]=FCS_Threshold(Inorm,Thrg1, Type,minseguser,maxseguser,colmin,colmax,3);       
maskseg=mask;  
imagenorm=Inorm;
end    

for kk=1:length(Hlog)
Htemp =  Hseg2(:,:,kk);
Hav(kk) =  mean(Htemp(maskseg>0));
end
Iav=mean(I(maskseg>0));
G0=Hav(1);
factor = (1 - (Gbkgdt/G0)*(Ibkgd/Iav)^2) / (1 - Ibkgd/Iav )^2  ;  % correction factor in presence of background 
% seg1=segmin;
% seg2=maxseg;

% Havcorr = factor .* mean(Hseg(seg1:seg2,:),1);
Havcorr = factor .*Hav ;
subplot(3,1,1)
title('Normalized segment intensity')
subplot(3,1,2)
title('Selected segments')
subplot(3,1,3)
semilogx(taulog, Hav,'-', taulog, Havcorr,'-.');
title('Average ACF')


prompt2 = {' Threshold (0-1):','Low or High','First column', 'Last column','First Row', 'Last Row','Segments'}; 
dlg_title2 = 'Input parameters'; 
num_lines = 1;
def2 = {num2str(Thrg1),Type,num2str(colmin),num2str(colmax),num2str(minseguser),num2str(maxseguser),num2str(subseg)};
answer2 = inputdlg(prompt2,dlg_title2,num_lines,def2);
figcheck=~isempty(answer2); 
if figcheck==1
Thrg1=str2num(answer2{1});
Type=(answer2{2});
colmin=str2num(answer2{3});
colmax=str2num(answer2{4});
minseguser=str2num(answer2{5});
maxseguser=str2num(answer2{6});
subseg=str2num(answer2{7});
end
end


answer = questdlg('Save data?');
if strcmp(answer,'Yes')==1
prompt = {'File name'}; 
dlg_title = ''  ; 
num_lines = 1;
def = {filenamenoext };
answerfile = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answerfile) == 0
filenamenoext=(answerfile{1});
end
    
ACFout1=cat(2, taulog' , Hav');
dlmwrite([filenamenoext,'_ACF',tag,'_',Type,'.txt'],ACFout1,'delimiter',';','precision',4);
if filenamepatt ~= 0
ACFout1corr=cat(2, taulog' , Havcorr');
dlmwrite([filenamenoext,'_ACF_corr',tag,'_',Type,'.txt'],ACFout1corr,'delimiter',';','precision',4);
end
%save data in Matlab
save([filenamenoext,'.mat'] ,  'Iav' , 'Hav' );


%save images
%export images Inorm', the selection masks
A1=uint16(1000*imagenorm);
outputFileName = [filenamenoext,'-imagenorm.tif'];
imwrite(A1, outputFileName);
A2=uint16(mask);
outputFileName = [filenamenoext,'-mask.tif'];
imwrite(A2, outputFileName);

end



end


%FUNCTIONS


% function read
function A=ISFCS_readfiletif(fname)
% fname = [filename, '.tif'];
info = imfinfo(fname);
nslice = numel(info);
A=imread(fname, 1);  % read tif stack data
for k = 2:nslice
    B = imread(fname, k);
    A=cat(3,A,B);
end
end



%function for FCS calc
function [Hlog, taulog, Hlin, tau]=ISFCS_imgdata(A,dt,plot,Nlog,d)
% given an input matrix A of dimension Y,T containing intensity corrsponding to 1)spatial
% position and 2)time (intensity 'carpet') yields:
% Hlog: matrix of dimension Y,circa-Nlog containing the ACF sampled in log
% scale for each pixel (ACF carpet); taulog: the values of correlation
% time;
% d is the distance between column: d=0 ACF, d>0 pCF

dimimg=size(A);
Y=dimimg(1);
T=dimimg(2);
H=zeros(Y,T);

 for k=1:Y-d       % calc ACF curve via FFT for each position k
 H(k,:)=ifft(fft(A(k,:)).*conj(fft(A(k+d,:))))/((mean(A(k,:)))*(mean(A(k+d,:)))*T)-1; 
 end

Hlin=H(:,2:floor(T/2+1));
tau=dt:dt:(T/2)*dt;

% get log scale and rebin data in time-log scale
tlog0=[1:29];
nplg=(Nlog)/log10(T/2);
tlog=cat (2,tlog0, (10.^(1/nplg)).^[nplg*log10(30):nplg*log10(T/2)-1] ) ;
tindx=uint64(round(tlog));
Hlogav=ones(Y,29);
for j=1:29
  Hlogav(:,j)= Hlin(:, tindx(j) ) ;
end
for i=30:length(tindx)-1
Hlogav=cat(2, Hlogav, mean(Hlin(:, tindx(i-1):tindx(i+1) ),2 )  ) ;
end
Hlogav=cat(2, Hlogav, mean(Hlin(:, tindx(length(tindx)-1):tindx(length(tindx))),2 )   ) ;

taulog=dt*tlog;
Hlog=Hlogav;

if plot==1
% figure
% plot(log10(tau),Hlin);
% axis([log10(dt) log10(T*dt/2) 0.0 0.1]);
figure
semilogx(tlog,Hlogav)
end

end

function [Hlog, taulog, Hlin, tau]=ISFCS_imgdata_X(A,A2,dt,plot,Nlog,d)  % check
% given 2 input matrix A, A2 of dimension Y,T containing intensity corrsponding to 1)spatial
% position and 2)time (intensity 'carpet') yields:
% Hlog: matrix of dimension Y,circa-Nlog containing the ACF sampled in log
% scale for each pixel (ACF carpet); taulog: the values of correlation
% time;
% d is the distance between column: d=0 ACF, d>0 pCF

dimimg=size(A);
Y=dimimg(1);
T=dimimg(2);
H=zeros(Y,T);

 for k=1:Y-d       % calc ACF curve via FFT for each position k
 H(k,:)=ifft(fft(A(k,:)).*conj(fft(A2(k+d,:))))/((mean(A(k,:)))*(mean(A2(k+d,:)))*T)-1; 
 end

Hlin=H(:,2:floor(T/2+1));
tau=dt:dt:(T/2)*dt;

% get log scale and rebin data in time-log scale
tlog0=[1:29];
nplg=(Nlog)/log10(T/2);
tlog=cat (2,tlog0, (10.^(1/nplg)).^[nplg*log10(30):nplg*log10(T/2)-1] ) ;
tindx=uint64(round(tlog));
Hlogav=ones(Y,29);
for j=1:29
  Hlogav(:,j)= Hlin(:, tindx(j) ) ;
end
for i=30:length(tindx)-1
Hlogav=cat(2, Hlogav, mean(Hlin(:, tindx(i-1):tindx(i+1) ),2 )  ) ;
end
Hlogav=cat(2, Hlogav, mean(Hlin(:, tindx(length(tindx)-1):tindx(length(tindx))),2 )   ) ;

taulog=dt*tlog;
Hlog=Hlogav;

if plot==1
% figure
% plot(log10(tau),Hlin);
% axis([log10(dt) log10(T*dt/2) 0.0 0.1]);
figure
semilogx(tlog,Hlogav)
end

end

% function threshold
function [B,C]=ISFCS_Threshold(A,thr1,thr2,colmin,colmax)
B1=A;
C=A;
%condition on intensity
B1(A<=thr1)=0;
B1(A>thr1)=1;
B=1-B1;
C(A<thr2)=0;
C(A>=thr2)=1;
%add condition on columns
if colmin<=colmax
    if colmin > 1
    B(:,1:colmin)=0;
    C(:,1:colmin)=0;
    end
    if colmax < size(A,2)
    B(:,colmax:end)=0;
    C(:,colmax:end)=0;
    end
else
    B(:,colmax:colmin)=0;
    C(:,colmax:colmin)=0;	
end
% figure
subplot(2,2,3)
imagesc(A,[min(nonzeros(A)),max(nonzeros(A))])
% axis image
subplot(2,2,1)
imagesc(B)
colormap(hot)
% axis image
subplot(2,2,2)
imagesc(C)
colormap(hot)
% axis image

end

function [param, fval]=ISFCS_Fit_single_gauss2D(x,y,Ginf0, G00, taud0,Display)
% model: Gauss PSF,2D, free diffusion


my=min(y);
My=max(y);
G00=min(G00,My);

[param, fval]=fminsearch(@(Param) sum( ( (Param(1)+(Param(2)./(1+x./Param(3)))) -y ).^2 ),[Ginf0 G00 taud0]);

if Display==1
%     figure;
    semilogx(x,y,'s')
    hold on
    semilogx(x,(param(1)+(param(2)./(1+x./param(3)))),'--r')
    hold off
%     param(1)=param(1)./param(2);
    title(strcat('taud=',num2str(param(3),1),'  G0=',num2str(param(2),1)))
end
end


function [param, fval]=ISFCS_Fit_single_gauss2D_dist(x,y,Ginf0, G00, taud0,distw, Display)
% model: Gauss PSF,2D, free diffusion - dual foci

my=min(y);
My=max(y);
G00=min(G00,My);

[param, fval]=fminsearch(@(Param) sum( ( (Param(1)+(Param(2)./(1+x./Param(3))).* (exp(-(distw^2)./(1+x./Param(3)))) ) -y ).^2 ),[Ginf0 G00 taud0]);

if Display==1
%     figure;
    semilogx(x,y,'s')
    hold on
    semilogx(x,(param(1)+(param(2)./(1+x./param(3))).*(exp(-(distw^2)./(1+x./param(3)))) ),'--r')
    hold off
%     param(1)=param(1)./param(2);
    title(strcat('taud=',num2str(param(3),1),'  G0=',num2str(param(2),1)))
end
end




function [Idet ]=Test_DetrendData(I,detrmode, Ldet,Lav, pow, noise)
%detrend data I(y,t) along t 
% 
% 
Lavh=floor(Lav/2);
% Ism=I;
dimimg=size(I);
Y=dimimg(1);
T=dimimg(2);
% TotTime=dt*T;
% tmax=min(tmax,TotTime);
% pseg = nextpow2(tmax/dt);
% Lseg =round( 2.^(pseg-1) );
maxseg=floor(T/Ldet);
% H=zeros(Y,Lseg);
switch detrmode
    case 'polyfit'
%using polyfit
for k=1+Lavh:Y-Lavh
Linex=double(mean(I( k-Lavh:k+Lavh , : ),1));
xx=1:length(Linex) ; 
xx=xx'; 
yy=Linex';
p = polyfit(xx,yy,pow);
yf = polyval(p,xx);
% figure; plot( xx, yy) ;
% figure; plot(xx,yf1 ) ;
yf1=max(yf)-yf;
Idet(k , : ) = double( I ( k,:)) + yf1' ; 
end
 
    case 'segment'
% for k=1+Lavh:Y-Lavh
% yf = movmean(I( k , : ),Ldet) ;
% yf1=max(yf)-yf;
% Idet(k , : ) = double( I ( k,:)) + yf1 ; 
% end

%using segments
for k=1+Lavh:Y-Lavh
for i=1:maxseg
A=sum( I( k-Lavh:k+Lavh, 1+ Ldet*(i-1):Ldet*(i) ) , 1 );
Iav=mean(A);
Ism(k, 1+ Ldet*(i-1):Ldet*(i) )=double( ones(1,Ldet).*Iav/Lav );
% I2(i)=mean(A(i2,:));
  % calc cross-CF curve via FFT
% H(i,:)=ifft(fft(A(i1,:)).*conj(fft(A(i2,:))))/((mean(A(i1,:)))*(mean(A(i2,:)))*Lseg)-1; 
end
Imax=max(Ism(k,:));
Iadd=-Ism(k,:)+Imax;
switch noise
    case 'Poisson'
    Iadd= 1e12 * imnoise(Iadd*1e-12,'poisson') ;  %require im proc toolbox
    Idet(k,:)=uint16(I(k,:))+uint16(Iadd);        %require im proc toolbox
    case 'Random'
    ErrAdd = 1 ./ sqrt(mean(Iadd, 'all')) ; % rel noise to add
    Idet(k,:)=double(I(k,:))+double(Iadd .*(1 + ErrAdd.*randn(size(Iadd)) ) );  % add random noise
end
end
Idet=Idet(1+Lavh:Y-Lavh,:);

end

end




function [B]=FCS_Threshold(A,thr1,Type,colmin,colmax,minseguser,maxseguser, plotsize)
B1=zeros(size(A));
B2=zeros(size(A));
A1=ones(maxseguser-minseguser+1,colmax-colmin+1);
B1(minseguser:maxseguser,colmin:colmax)=A1;
switch Type
    case 'Low'
B2(A>thr1)=0;
B2(A<=thr1)=1;        
    case 'High'
B2(A<=thr1)=0;
B2(A>thr1)=1;           
end

B=B1.*B2;

% figure
subplot(plotsize,1,1)
imagesc(A)
% axis image
subplot(plotsize,1,2)
imagesc(B)
% axis image


end
