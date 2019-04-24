function [ x_prima ,y_prima ,imagen ] = SAR_gprMax( filename )
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
index2=find(prof >= 0.27);
index=find(prof <= 2.27);
%field1 = field(index:end,:);%Quitamos los 15 primeros centrimetos de la medida
%field = field(index:end,:);%Quitamos los 15 primeros centrimetos de la medida
% media=zeros(size(field));
% media(index,:)=field(index,:);
%field = field - repmat(sum(media,2)/size(media,2),1,size(media,2));
traces = 0:size(field,2)-1;
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
k0=2*pi*f*sqrt(1)/3e8;
%k1=k0*sqrt(2.5);
%Definimos las coordenadas de la imagen
x_prima=traces*0.008;%0.025;
%y_prima=0.15:0.01:0.85-0.01;
y_prima=prof(index2);%0:0.001:0.30;

[X,Y]=meshgrid(x_prima,y_prima);%Para hacer la matriz de distancias
%Matriz de exponenciales
exponenciales=repmat(2j,length(y_prima),length(FIELD));
imagen=zeros(length(y_prima),length(x_prima));%Creacion de la matriz de la imagen SAR
for m=1:size(FIELD,2) %Se recorren las antenas
    %Creamos la matriz distancias dependiente de la antena que se este midiendo3 
    D=sqrt(abs(X-X(1,m)).^2+(Y).^2);
    for x=1:length(x_prima)%Se recorren las columnas de la matriz imagen
%         imagen(:,x) = imagen(:,x) + exp(2*1j*(2*pi*repmat(f,length(y_prima),1).*...
%             repmat(D(:,x), 1, length(f)))/3e8)*FIELD(:,m);
        omega=repmat(k0,length(y_prima),1).*repmat(D(:,x),1,length(FIELD));%Nos quedamos con las distancias de la columna x
        imagen(:,x)=imagen(:,x)+exp(exponenciales.*omega)*FIELD(:,m);
        
    end 


end

