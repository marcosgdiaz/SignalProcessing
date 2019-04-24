function [ x_prima ,z_prima ,refl ] = SAR_gprMax_PSM_multi(filename, lim, eps )
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
index2=find(prof >= 0.27 & prof <= 0.77);
index=find(prof <= 0.3);
%field = field(index:end,:);%Quitamos los 15 primeros centrimetos de la medida
media=zeros(size(field));
media(index,:)=field(index,:);
field = field - repmat(sum(media,2)/size(media,2),1,size(media,2));
%field = field - repmat(min(media,[],2),1,size(media,2));

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
S=permute(S,[2 3 1]);
Sk=fft(S,[],1);
Sk=fftshift(Sk,1);
% Sk=fft(Sk,[],2);
% Sk=fftshift(Sk,2);

%x_prima=traces*0.025;
x_prima=traces*0.008;
y_prima=0;
%z_prima=0.15:0.01:0.85-0.01;
z_prima=prof(index2);

ksx=2*pi/0.008;
kx=linspace(-ksx/2,ksx/2,length(x_prima));
ky=0;

%kz1=sqrt(4*(2*pi*permute(repmat(f,length(kx),1),[1,3,2])/3e8).^2-permute(repmat(kx,length(f),1),[2 3 1]).^2);
kz=sqrt(4*(2*pi*permute(repmat(f,length(kx),length(ky),1,length(eps)),[1,3,2,4]).*real(sqrt(repmat(permute(eps,[4 3 1 2]),...
    length(kx),length(ky),length(f))))/3e8).^2-permute(repmat(kx,length(f),length(ky),1,length(eps)),[2 3 1 4]).^2);


% T12=2*eps*kz1./(eps*kz1+1*k2z);
% T21=2*1*kz2./(1*kz2+eps*kz1);
% kz1(kz1.*kz1 < 0) = 0;


kz(kz.*kz < 0) = 0; %Est?n ordenados en la cuarta dimension

refl=zeros(numel(x_prima),length(y_prima),numel(z_prima));

for i=1:numel(z_prima)
    for j=1:numel(f)
        log = z_prima(i) >= lim;
        pos = find(log,1,'last');
        if pos>1
            refl(:,:,i)=refl(:,:,i)+Sk(:,:,j).*prod(exp(1j*kz(:,:,j,1:pos-1)...
                .*repmat(permute(lim(2:pos)-lim(1:pos-1),[4 3 1 2]),size(kz(:,:,1,1)))),4).*exp(1j*kz(:,:,j,pos)*(z_prima(i)-lim(pos)));
        else
            refl(:,:,i)=refl(:,:,i)+Sk(:,:,j).*exp(1j*kz(:,:,j,pos)*z_prima(i));
        end
%         if z_prima(i)<=lim
%             refl(:,:,i)=refl(:,:,i)+Sk(:,:,j).*exp(1j*kz1(:,:,j)*z_prima(i));
%         else
%             refl(:,:,i)=refl(:,:,i)+Sk(:,:,j).*exp(1j*kz1(:,:,j)*lim).*exp(1j*kz2(:,:,j)*(z_prima(i)-lim))./exp(-1j*angle(T12(:,:,j).*T21(:,:,j)));
%         end
    end
end
refl=ifftshift(refl,1);
refl=ifft(refl,[],1);
% refl=ifft(refl,[],2);
% refl=ifftshift(refl,2);
refl=permute(refl,[3 1 2]);

end



