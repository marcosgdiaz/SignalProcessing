function [x_prima ,y_prima ,imagen] = SAR_gprMax_multi(filename, lim, eps)
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
prof = time*3e8/2;
index2=find(prof >= 0.27 & prof <= 0.77);
index=find(prof <= 0.30);
%field1 = field(index:end,:);%Quitamos los 15 primeros centrimetos de la medida
media=zeros(size(field));
media(index,:)=field(index,:);
field = field - repmat(sum(media,2)/size(media,2),1,size(media,2));
traces = 0:size(field,2)-1;
%field = field(index:end,:);%Quitamos los 15 primeros centrimetos de la medida
%Definimos el eje de frecuencias
NFFT=2^nextpow2(length(field));%Numero de puntos de la fft optimo
Fs=1/header.dt;%Frecuencia de muestreo del GprMax
f=linspace(-Fs/2,Fs/2,NFFT);
f=f(f>=0);%Frecuencias positivas
%Hacemos fftshift que centra la frecuencia cero en dimension 1(por
%columnas)
FIELD1=fft(field,NFFT,1);
FIELD=FIELD1(1:length(f),:);%Nos quedamos en el campo en frecuencias positivas
%Preparamos la matriz que va a multiplicar al campo
k0=2*pi*f/3e8;
k=[k0*sqrt(eps(1)); k0*sqrt(eps(2)); k0*sqrt(eps(3))];
%Definimos las coordenadas de la imagen
x_prima=traces*0.008;%0.025;
%y_prima=0.15:0.01:0.85-0.01;
y_prima=prof(index2);

[X,Y]=meshgrid(x_prima,y_prima);%Para hacer la matriz de distancias
%Matriz de exponenciales
imagen=zeros(length(y_prima),length(x_prima), numel(k0));%Creacion de la matriz de la imagen SAR

D1=zeros(length(y_prima),length(x_prima),size(FIELD,2));
D2=D1;
D3=D1;
for m=1:size(FIELD,2)
    for i=1:length(y_prima)
        if y_prima(i) <= lim(1)
           D1(i,:,m)=sqrt((x_prima(m)-X(i,:)).^2+y_prima(i).^2);         
        elseif y_prima(i) > lim(1) && y_prima(i) <= lim(2)
           x2 = x_prima(m) - X(i,:);
           ratio=sqrt(eps(1)/(eps(2)));%Ponderacion dependiente de Er2
           xb=x2+ratio*(x2*lim(1)./y_prima(i)-x2);%Distancia x del la ruta optima
           D1(i,:,m)=sqrt(lim(1)^2+xb.^2);
           D2(i,:,m)=sqrt((y_prima(i)-lim(1)).^2+(x2-xb).^2);
        else
           x3 = x_prima(m) - X(i,:);
           num = (1 + sqrt(eps(2)/eps(3)) * (lim(2) - y_prima(i)) / (y_prima(i) - lim(1)))*x3;
           denom = (1 - sqrt(eps(2)/eps(3)) * (y_prima(i) - lim(2)) / (y_prima(i) - lim(1)) - sqrt(eps(1)/eps(3)) * (y_prima(i) - lim(2)) * (lim(1) - lim(2)) / ((y_prima(i) - lim(1)) * lim(2)));
           xb2 =  num / denom;
           xb1 = xb2 * (1 + sqrt(eps(1)/eps(2)) * (lim(1) - lim(2)) / lim(2));
           D1(i,:,m)=sqrt(lim(1)^2+xb1.^2);
           D2(i,:,m)=sqrt((lim(2)-lim(1))^2+(xb2-xb1).^2);
           D3(i,:,m)=sqrt((y_prima(i)-lim(2)).^2+(x3-xb2).^2);
        end
    end
end
       
bp=waitbar(0,'Procesado DAS proyeccion 3D'); 
for m=1:size(FIELD,2)%Se recorren las antenas
   %disp(['Medida m = ' num2str(m) ' de ' num2str(Nx*Ny)]);
       for n=1:numel(k0) 
         S=FIELD(n,:);
         imagen(:,:,n)=imagen(:,:,n)+S(m)*exp(2j*(k(1,n)*D1(:,:,m)+k(2,n)*D2(:,:,m)+k(3,n)*D3(:,:,m))); %*(4*pi*2*D*f(n)/(3e8*2*D0int(n))).
       end
   
   waitbar(m/(size(FIELD,2)),bp,sprintf('%0.2f %%',100*m/(size(FIELD,2))));
end
close(bp)
imagen = sum(imagen,3);

end
      