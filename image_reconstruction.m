%-----------------------------


%

clear;
[filename, pathname] = uigetfile( {'*.oct', 'OCT Files (*.oct)'; ...
 '*.*','All Files (*.*)'},  'Pick an  Oct File');
if isequal(filename,0) || isequal(pathname,0),
 return;
end
fpath=[pathname filename];
fid = fopen(fpath, 'r');
octdata=fread(fid,'ushort');
sta=fclose(fid);
pointsnumerA=1344;
pointsnumerB=1250;
pointsnumerC=1344;
data=zeros(1000,pointsnumerC);
for i=1:1000
    for j=1:pointsnumerC
    data(i,j)=octdata(pointsnumerC*1000*0+(i-1)*pointsnumerC+j);
    end
end
 rawdata=data(1:1000,1:pointsnumerA);
mean_data=mean(rawdata(:,1:pointsnumerA));
padzeros=zeros(1,704);
hilbert_value=zeros(1000,pointsnumerA);
padzeros_value=zeros(1000,2048);
imagetrans=zeros(1000,2048);
image_data=zeros(1000,2048);
%  a=-0.00000008;
%  b=0.00000002;
%  c=-0.00000000;% far
 a=0.000028;
 b=-0.00000001;
 c=-0.000000003;%near
win2=hanning(1344);
win=hanning(2048);
for i=1:1000
    rawdata(i,:)=rawdata(i,:);%;-mean_data;
    hilbert_value(i,1:2048)=[rawdata(i,1:pointsnumerA).*win2' padzeros];
    for ii=1:2048
        hilbert_value(i,ii)=hilbert_value(i,ii)*exp((-1i*((a+b*ii)*(ii-pointsnumerA/2)^2))+(-1i*((c)*(ii-pointsnumerA/2)^3)));
    end
    padzeros_value(i,1:2048)=(hilbert_value(i,1:2048)).*win';
    imagetrans(i,1:2048)=abs(fft(padzeros_value(i,1:2048)));
    image_data(i,1:2048)=20*log10(abs(fft(padzeros_value(i,1:2048))));   
end
figure,imshow(image_data/155);
%hold on;
%plot(imagetrans(800,20:1024));
figure;plot(image_data(800,20:1024));

var_value=var(image_data(800,1:1024));

