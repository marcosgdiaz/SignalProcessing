clear all, close all, clc
root='/Users/marcosgonzalezdiaz/Documents/Universidad/TFG/gprMax-master/Ensayos/';
A='primer_modelo_merged.out';
B='segundo_modelo_merged.out';
C='tercer_modelo_merged.out';
D='cuarto_modelo_merged.out';
E='quinto_modelo_merged.out';
F='sexto_modelo_merged.out';
G='varios_objetos_merged.out';
H='varios_objetos1_merged.out';
I='interfaz_rugosa_merged.out';
J='heterogeneo_plano_merged.out';
K='heterogeneo_rugoso_merged.out';
L='paper_TFG_hf_merged.out';
M='prueba_mortero_merged.out';
N='prueba2_merged.out';
O='tfg1_merged.out';
P='tfg2_merged.out';
Q='tfg3_merged.out';
R='tfg4_merged.out';
S='tfg5_merged.out';
T='tfg6_merged.out';
U='paper1_merged.out';


%[x, y, imagen]=SAR_gprMax_correction(strcat(root,P));
[x, y, imagen]=SAR_gprMax_PSM_multi(strcat(root,O),[0 0.42 0.52], [1 2.5 3.5]);%[0 0.14],[1 4+0.02j]); 
%[0 0.072 0.09 0.14],[1 2.19-0.02j 1.1-0.001j 4-0.2j]);%[0 0.072 0.09 0.14],[1 2.19-0.02j 1.1-0.001j 4+0.2j]);% 0.09 0.14],[1 2.19-0.02j 1.1-0.001j 4+0.2j]);
%pcolor(x*100-24,y*100-7,20*log10(abs(imagen/max(max(imagen)))))
find(y>=0.17 & y<=0.77);
% index=find(y<=100);
% media=zeros(size(imagen,1),1);
% media(index)=sum(imagen(index,:),2)/size(imagen,2);
% imagen = imagen - repmat(media,1,size(imagen,2));
pcolor(x*100-24,y(ans)*100-7,20*log10(abs(imagen(ans,:)/max(max(imagen(ans,:))))))
fontSize = 24;
set(gca,'YDir','reverse')
colormap jet;
shading interp;
xlabel('x [cm]','fontSize',fontSize);
ylabel('z [cm]','fontSize',fontSize);
h = colorbar;
h.Location='southoutside';
ylabel(h, 'Amplitud normalizada [dB]', 'fontSize',fontSize)
axis equal;
xlim([-24 24])
ylim([20 70])
caxis([-30 0]);
set(gca,'fontsize',fontSize)

