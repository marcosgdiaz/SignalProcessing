clear all; close all; clc
ImagePartition('test03.jpg',210, 600, 11500 )
%load('plano_met.mat')
%load('objetos_enterrados.mat')
%load('hormigon_1layer_75.mat')
%load('hormigon_2layer_55.mat')
%load('tierra_placaEncimayEnterrada_z150.mat')

% load('../RadarILMsens/imagen_3D_arena_plastico.mat')
% y=linspace(-0.5,0.5,41);
% x=linspace(-0.4,0.4,36);
% z=linspace(0,1.5,50);
% Image3DPartition(x,y,z,20*log10(abs(imagen/max(max(max(imagen))))),-15.5,50, 2400)
