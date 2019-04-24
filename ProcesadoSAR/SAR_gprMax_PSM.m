function [ x_prima ,z_prima ,refl ] = SAR_gprMax_PSM( filename )
header.iterations = double(h5readatt(filename,'/', 'Iterations'));
header.dt = h5readatt(filename, '/', 'dt');
fields.time = linspace(0, (header.iterations)*(header.dt)*1E9, header.iterations);
%   Detailed explanation goes here
fields.ex = h5read(filename, strcat('/rxs/rx1/', 'Ex'))';
fields.ey = h5read(filename, strcat('/rxs/rx1/', 'Ey'))';
fields.ez = h5read(filename, strcat('/rxs/rx1/', 'Ez'))';
fields.hx = h5read(filename, strcat('/rxs/rx1/', 'Hx'))';
fields.hy = h5read(filename, strcat('/rxs/rx1/', 'Hy'))';
fields.hz = h5read(filename, strcat('/rxs/rx1/', 'Hz'))';
%Cargamos campo, tiempo y el numero de traza
field = fields.ez;
%field1=field(200:end,:)';%Quitamos el acoplo. Para la fft, tienen que estar en filas
time = 0:header.dt:header.iterations * header.dt;
prof = time(1:(end-1))*3e8/2;
index2=find(prof >= 0.17);
index=find(prof <= 2.27);
%field = field(index:end,:);%Quitamos los 15 primeros centrimetos de la medida
% media=zeros(size(field,1),1);
% media(index)=sum(field(index,:),2)/size(field,2);
%field = field - repmat(media,1,size(field,2));
traces = 0:size(field,2)-1;
%Definimos el eje de frecuencias
NFFT=2^nextpow2(length(field));%Numero de puntos de la fft optimo
Fs=1/header.dt;%Frecuencia de muestreo del GprMax
f=linspace(-Fs/2,Fs/2,NFFT);
f=f(f>=0);%Frecuencias positivas
%Hacemos fftshift que centra la frecuencia cero en dimension 1(por
%columnas)
S1=fft(field,NFFT,1);
S=S1(1:length(f),:);%Nos quedamos en el campo en frecuencias positivas
% S=S(1:3:end,:);
% f=f(1:3:end);
S=permute(S,[2 3 1]);
Sk=fft(S,[],1);
Sk=fftshift(Sk,1);
% Sk=fft(Sk,[],2);
% Sk=fftshift(Sk,2);

%x_prima=traces*0.025;
x_prima=traces*0.008;
y_prima=0;
%z_prima=0:0.01:1-0.01;
z_prima=prof(index2);%200:4:900);
% refl=field(index2,:);
% return;
ksx=2*pi/0.008;
kx=linspace(-ksx/2,ksx/2,length(x_prima));

kz=sqrt(4*(2*pi*permute(repmat(f,length(kx),1),[1,3,2])*sqrt(1)/3e8).^2-permute(repmat(kx,length(f),1),[2 3 1]).^2);
kz(kz.*kz < 0) = 0;
refl=zeros(numel(x_prima),length(y_prima),numel(z_prima));

for i=1:numel(z_prima)
    for j=1:numel(f)
        refl(:,:,i)=refl(:,:,i)+Sk(:,:,j).*exp(1j*kz(:,:,j)*z_prima(i));
    end
end
refl=ifftshift(refl,1);
refl=ifft(refl,[],1);

% refl=ifft(refl,[],2);
% refl=ifftshift(refl,2);
refl=permute(refl,[3 1 2]);

end

