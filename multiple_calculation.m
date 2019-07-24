% Calculation of all experimental data and writing it to one excel file 
clear all;
close all;
clc;

global y counter1 Name1 Container Experiment_number PlotAll

Experiment_number = 13; % Experiment Number, u can check it in Name1
PlotAll = 0; % to plot all the experiments enter 1

Name1 =  {'_04_06.txt','_02_07.txt','_06_07.txt','_13_07.txt','_13(2)_07.txt','_02_08.txt','_07_08.txt',...
          '_08_08.txt','_13_08.txt','_13(2)_08.txt','_14_08.txt','_28_08.txt','_08_10.txt','_08(2)_10.txt','_08(3)_10.txt','_08(4)_10.txt',...
          '_18_10.txt','_18_10(2).txt','_18_10(3).txt','_18_10(4).txt','_18_10(5).txt','_18_10(6).txt','_18_10(7).txt','_18_10(8).txt',...
          '_08_11.txt','_08(2)_11.txt','_08(3)_11.txt','_08(4)_11.txt','_08(5)_11.txt','_08(6)_11.txt','_08(7)_11.txt';...
    'ДРШ_250_4.11.04','ДРШ_250_4.11.04','ДРШ_250_4.11.04','ДРШ_250_4.04.16','ДРШ_250_4.04.16','ДРКС_200_3.03','ДРКС_200_3.03','ДРШ_250_СССР_митино',...
    'ДРШ_250_СССР_митино','ДРШ_250_СССР_митино','ДРШ_250_2M_Разряд','Osram_HBO_w&by','Osram_HBO_w&by','Osram_HBO_w&by','Osram_HBO_w&by','Osram_HBO_w&by',...
    'ДРШ_250_66','ДРШ_250_66','ДРШ_250_66','ДРШ_250_66','ДРШ_250_66','ДРШ_250_66','ДРШ_250_66','ДРШ_250_66',...
    'ДРKC_200_27','ДРKC_200_27','ДРKC_200_27','ДРKC_200_27','ДРKC_200_27','ДРKC_200_27','ДРKC_200_27'};% file names for all experiments
N = length(Name1) ;
counter1 = 0;
for y = 1:N    
clear 
%clear - except : y, counter1, Name1
clearvars variables 
close all
clc

getGlobals

%tanf
%hx
load('Vlyambda.mat') % spectral luminous efficiency function
step=5;
range = 0:1100; % range of interpolation
stitching = 450;  % the point where the connection of 2 graphics (with filtra and without) takes plase 
lyambdamin1 = 200;% 1 left boundary of integration 
lyambdamax1 = 300;% 1 right boundary of integration
lyambdamin2 = 200;% 2 left boundary of integration
lyambdamax2 = 800;% 2 right boundary of integration                        

Name = char(Name1(1,y)) ;% File name
%Name = '_02_07.txt' ;
id = fopen(['Exp',Name], 'r'); % obtain information about open files
fgetl(id);   % to skip first line

%Fspec = textscan(fid,'%f',10); another method read formatted data from text file
% N = 34; number of columns
% Fr = '';
% for i = 1:N
%     Fr = [Fr, '%f'];
% end
% Fspec = textscan(fid, Fr,'Delimiter', ' ');
%%  reading of experimental data taken by spectrometer

Fspectra = textscan(id, '%f','Delimiter', ' '); % read formatted data from text file
Fspectra = Fspectra{1}; % open a container
Flag=0;
f=0;
N=0;
while Flag < 1 % calculating the number of columns in a file
    f=f+1;
    if Fspectra(f) ~= 1 % determine the second row
    N=N+1;
    else
    Flag = 1;
    end
end

Fspectra1 = [];
i = 1;
j = 1;
while (i+N-1) <= length(Fspectra) % create a container with the data in the correct format
    Fspectra1(j,:) = Fspectra(i:i+N-1); % a data matrix
    i = i + N;
    j = j + 1;
end
Fspectra2 = cell(1,N);
for i = 1:N
    Fspectra2{i} = Fspectra1(:,i); % an actual container
end
fclose(id);
Fspectra = Fspectra2;
%% reading of experimental data 
id = fopen(['Lux',Name], 'r');% -4  % obtain information about open files (illumination)
fgetl(id);   %to skip first line
Evv = textscan(id, '%f', 'Delimiter', ' '); % read formatted data from text file
Ev = cell2mat(Evv); % convert cell array to ordinary array
fclose(id);

%Ev = load('Luxmetr.txt') ;
k = length(Ev) ;

counter1 = counter1 + 4 ;

id = fopen(['Temperature',Name], 'r');% -4  % obtain information about open files (temperature)
fgetl(id);   %to skip first line
Tc = textscan(id, '%f', 'Delimiter', ' '); % read formatted data from text file
T = cell2mat(Tc); % convert cell array to ordinary array
fclose(id);

id = fopen(['I',Name], 'r');% -4  % obtain information about open files (current)
fgetl(id);   %to skip first line
Ia = textscan(id, '%f', 'Delimiter', ' '); % read formatted data from text file
I = cell2mat(Ia); % convert cell array to ordinary array
fclose(id);

id = fopen(['U',Name], 'r');% obtain information about open files (voltage)
fgetl(id);   %to skip first line
Uv = textscan(id, '%f', 'Delimiter', ' '); % read formatted data from text file
U = cell2mat(Uv); % convert cell array to ordinary array
fclose(id);

id = fopen(['Ufans',Name], 'r');% obtain information about open files (voltage)
fgetl(id);   %to skip first line
Uf = textscan(id, '%f', 'Delimiter', ' '); % read formatted data from text file
Ufans = cell2mat(Uf); % convert cell array to ordinary array
fclose(id);

for i = 1:k 
    Power(i) = I(i)*U(i); % power calculation
end
%% Calculation of solid angles
for i = 1:k
    
diameter = 0.024; % the diameter of the receiver
l = 55*diameter; % distance from source to receiver
L = sqrt(l^2+(diameter/2)^2); % hypotenuse 

fi = acosd(l/L) ; % opening angle

% tanf % the tangent of the angle (Mathcad)
% hx  % the size of the aperture (Mathcad)

SolidAngle = 2*pi*(1-cosd(fi)); % solid angle
SolidAngleall = 2*pi*(1-cosd(90));

Fv(i) = (pi*diameter^2)*Ev(i)/4; % lighting flow in (lm) 
end
%% Main calculation
 ji=0;
VLyambda{j} = interp1(lyambda,Vlyambda,range); % interpolation for spectral sensitivity
for i = 1:k  % calculation of the energy flow  
    
 if i>1 
    ji=ji+1;
 end
Fholder1 = Fspectra{i+4+ji}; % spectrum with a filter
Fholder2 = Fspectra{i+5+ji}; % spectrum without filter

flag = 0;
f = 0;
while flag < 1 
    f=f+1;
 if Fspectra{2}(f) < stitching % calculation of the point where the value of Fspectra{2} (wavelength) reaches the value of stitching
 else flag=1;
 end 
end     

Fholder = Fholder1(1:f); 
Length=length(Fspectra{2});
Fholder(f:Length) = Fholder2(f:Length).*1;
Fspectra{i+4} = Fholder; % stitched spectrum with and without filter
clear Fdouble1
clear Fdouble2
clear Fholder

FspectraInterp{i} = interp1(Fspectra{2},Fspectra{i+4},range,'linear',0); % interpolation for each spectrum
FeVlFunction = FspectraInterp{i}.* VLyambda{j}; % multiplication by the spectral luminous efficiency function
FeVlFunction(isnan(FeVlFunction))= 0; % if nan = 0
%VQ11 = VQ(~isnan(VQ));

Fvrelative(i) = 680*trapz(FeVlFunction(400:750)); % relative units
Kc(i) = Fv(i)/Fvrelative(i); %  coupling coefficient of the lighting and energy systems 
Fesolidangle1(i) = Kc(i)*trapz(FspectraInterp{i}(lyambdamin1:lyambdamax1)); % the energy flow that falls on the receiver (in the first range)
Fesource1(i) = 2*Fesolidangle1(i)/(1-cosd(fi)); % the energy flow from a source (in the first range)
Fesolidangle2(i) = Kc(i)*trapz(FspectraInterp{i}(lyambdamin2:lyambdamax2)); % the energy flow that falls on the receiver (in the second range)
Fesource2(i) = 2*Fesolidangle2(i)/(1-cosd(fi)); % the energy flow from a source (in the second range)

FeAllsolidangle(i) = Kc(i)*trapz(FspectraInterp{i}(200:800)); % the all energy flow that falls on the receiver(200:800)nm
FeAllsource(i) = 2*FeAllsolidangle(i)/(1-cosd(fi)); % the all energy flow from a source (200:800)nm

%calculation of the minimum flow ratio for the appliance (the minimum portion is 16(mW))

A = 0.004; % one side of the receiving field
B = 0.002; % another side of the receiving field
ll = 0.0168; % distance from source to receing field 
hypotenuse = sqrt(ll^2+(A/2)^2); % hypotenuse 
fangle = acosd(ll/hypotenuse); % opening angle
fefield(i) = (FeAllsource(i)*(1-cosd(fangle))/2)*0.69; % the energy flow that falls in a solid angle
Ferecevingfield(i) = fefield(i)*A*B/(pi*A^2/4); % the energy flow that falls on the receiving field (Fe*S_circle/S_rectangle)

Standard_diviation(i) = std(Fspectra{i+4}(lyambdamin1:lyambdamax2)); % Standard_diviation
end
%% finding extremums
MaxFFesolidangle1 = max(Fesolidangle1); % finding the maximum element in an array
MaxFesource1 = max(Fesource1); % finding the maximum element in an array
MaxFesolidangle2 = max(Fesolidangle2); % finding the maximum element in an array
MaxFesource2 = max(Fesource2); % finding the maximum element in an array
MaxFeAllsource = max(FeAllsource); % finding the maximum element in an array

Procent1 = 100*Fesource1./ FeAllsource;  % finding the percentage of the flow from the all flow (200:300)nm
Procent2 =  100*Fesource2./ FeAllsource;% finding the percentage of the flow from the all flow (200:800)nm
Procent1Max = 100*MaxFesource1 / MaxFeAllsource; % finding the percentage of the flow from the all flow (200:300)nm
Procent2Max = 100*MaxFesource2 / MaxFeAllsource; % finding the percentage of the flow from the all flow (200:800)nm

index1 = find(Fesolidangle1 == max(Fesolidangle1), 1); % finding the number of the maximum element
index2 = find(Fesource1 == max(Fesource1), 1); % finding the number of the maximum element
index3 = find(Fesolidangle2 == max(Fesolidangle2), 1); % finding the number of the maximum element
index4 = find(Fesource2 == max(Fesource2), 1); % finding the number of the maximum element

MaxT = max(T); % maximum temperature (for relative T)
MaxPower = max(Power); % maximum power (for relative Power)
Efficiency = (Fesource1*100)./Power;  % (in the first range)
Efficiencyall = (Fesource2*100)./Power ; % (in the second range)
%% plotting data for one experiment with extremums (is used in a program single_calculation)
% 
% disp(Ferecevingfield); % the energy flow that falls on the receing field
% disp(Fesource1);
% disp(Fesource2);
% 
% figure 
% grid on
% hold on
% for i= 1:k % cycle to plot all the stitched graphics
%      subplot(ceil(k/4),4,i);
%      plot(Fspectra{2},Fspectra{4+i});
%      title([num2str(4+i),' № ',num2str(i)])
% end 
% 
% figure 
% subplot(2,2,1)
% grid on
% hold on
% plot(1:k,Fesource1)
% set(gca,'XTick', [1:1:k])
% xlabel('Номер измерения')
% ylabel('Мощность излучения (Вт)')
% title(['Интегральный энергетический поток от времени в диапазоне ',num2str(lyambdamin1),'-',num2str(lyambdamax1),' нм']);
% plot(1:k,Power/MaxPower,'b') 
% plot(1:k,T*1/MaxT,'r') % relative units
% plot(1:k,Efficiency,'k') % relative units
% 
% subplot(2,2,3)
% grid on
% hold on
% plot(Fspectra{2},Fspectra{4+index1}*Kc(index1)); 
% set(gca,'XTick', [0:100:1100])
% xlabel('Длина волны (нм)')
% ylabel('Мощность излучения')
% title(['Спектр оптического излучения для Fe(Max) = ',num2str(MaxFesource1),' Вт = ',num2str(Procent1Max),'% от Fe(200-800)']); 
% 
% subplot(2,2,2)
% grid on
% hold on
% plot(1:k,Fesource2)
% set(gca,'XTick', [1:1:k])
% xlabel('Номер измерения')
% ylabel('Мощность излучения (Вт)')
% title(['Интегральный энергетический поток от времени в диапазоне ',num2str(lyambdamin2),'-',num2str(lyambdamax2),' нм']);
% plot(1:k,Efficiencyall,'k') 
%  
% subplot(2,2,4)
% grid on
% hold on
% plot(Fspectra{2},Fspectra{4+index2}*Kc(index2)); 
% set(gca,'XTick', [0:100:1100])
% xlabel('Длина волны (нм)')
% ylabel('Мощность излучения')
% title(['Спектр оптического излучения для Fe(Max) = ',num2str(MaxFesource2),' Вт = ',num2str(Procent2Max),'% от Fe(200-800)']); %Fspec{i+4}

%% Writing all data to excel

header = {'№_of_measurments','I','U','Power','Ev','Tc','Fesource1(200-300)','Fesource2(200-800)',...
    'Standard_diviation','Ufans','Efficiency', 'Efficiencyall','Procent1','Procent2',Name};
indicate = 0;

    if counter1 < 13 %
     Rows_down = ['B',num2str(3)];
    Rows_down2 = ['B',num2str(2)];
     
    else 
         Rows_down = ['B',num2str(1+counter1)];
    Rows_down2 = ['B',num2str(counter1)];
       
    end 
    
measurments = [1:k];
xlswrite('data.xlsx',[measurments(:),I(:),U(:),Power(:),Ev(:),T(:),Fesource1(:),Fesource2(:)...
    ,Standard_diviation(:),Ufans(:),Efficiency(:), Efficiencyall(:), Procent1(:), Procent2(:)],'Sheet1',Rows_down);
xlswrite('data.xlsx',header,'Sheet1',Rows_down2);

counter1=counter1+k;

%% Writing all data to the сontainer
  Container{y} = [measurments(:),I(:),U(:),Power(:),Ev(:),T(:),Fesource1(:),Fesource2(:)...
 ,Standard_diviation(:),Ufans(:),Efficiency(:), Efficiencyall(:), Procent1(:), Procent2(:)]; % an actual containe

end
%% plotting data for one experiment

figure('Name',[char(Name1(2,Experiment_number)),'  Date',char(Name1(1,Experiment_number))]) % new figure

for i= 1:6
ax(i) = subplot(6,1,i); 
hold(ax(i),'on')
grid(ax(i),'on');
end 

plot(ax(1),Container{Experiment_number}(:,1),Container{Experiment_number}(:,7),'-o') 
title(ax(1),'integral energy flow in W')
ylabel(ax(1),'Energy flow (W)')
xlabel(ax(1),'measurement number')

plot(ax(2),Container{Experiment_number}(:,1),Container{Experiment_number}(:,6),'-o')
title(ax(2),'Temperature of the bulb')
ylabel(ax(2),'Temperature C')
xlabel(ax(2),'measurement number')

plot(ax(3),Container{Experiment_number}(:,1),Container{Experiment_number}(:,10),'-o')
title(ax(3),'Cooling voltage')
ylabel(ax(3),'voltage')
xlabel(ax(3),'measurement number')

plot(ax(4),Container{Experiment_number}(:,1),Container{Experiment_number}(:,11),'-o')
title(ax(4),'Efficiency for UV')
ylabel(ax(4),'%')
xlabel(ax(4),'measurement number')

plot(ax(5),Container{Experiment_number}(:,1),Container{Experiment_number}(:,2),'-o')
title(ax(5),'Current')
ylabel(ax(5),'A')
xlabel(ax(5),'measurement number')

plot(ax(6),Container{Experiment_number}(:,1),Container{Experiment_number}(:,4),'-o')
title(ax(6),'Power')
ylabel(ax(6),'W')
xlabel(ax(6),'measurement number')

fig = gcf;
%fig.PaperFormat='A4';
%fig.PaperOrientation='landscape';
%fig.PaperSize=[29.7 21];
%fig.PaperPosition=[0.63, 0.2, 28.43, 20.6];
saveas(fig,[num2str(Experiment_number),'№_',char(Name1(2,Experiment_number)),'_Figure.fig'])

%% plotting for all experimental data

% for i= 1:k % cycle to plot all the stitched graphics
%      subplot(ceil(k/4),4,i);
%      plot(Fspectra{2},Fspectra{4+i});
%      title([num2str(4+i),' № ',num2str(i)])
% end 


if PlotAll > 0 
for i = 1:12
figure('Name',[char(Name1(2,i)),'  Date',char(Name1(1,i))]) % new figure
 
for j= 1:6
ax(j) = subplot(6,1,j); 
hold(ax(j),'on')
grid(ax(j),'on');
end 

plot(ax(1),Container{i}(:,1),Container{i}(:,7),'-o') 
title(ax(1),'integral energy flow in W')
ylabel(ax(1),'Energy flow (W)')
xlabel(ax(1),'measurement number')

plot(ax(2),Container{i}(:,1),Container{i}(:,6),'-o')
title(ax(2),'Temperature of the bulb')
ylabel(ax(2),'Temperature C')
xlabel(ax(2),'measurement number')

plot(ax(3),Container{i}(:,1),Container{i}(:,10),'-o')
title(ax(3),'Cooling voltage')
ylabel(ax(3),'voltage')
xlabel(ax(3),'measurement number')

plot(ax(4),Container{i}(:,1),Container{i}(:,11),'-o')
title(ax(4),'Efficiency for UV')
ylabel(ax(4),'%')
xlabel(ax(4),'measurement number')

plot(ax(5),Container{i}(:,1),Container{i}(:,2),'-o')
title(ax(5),'Current')
ylabel(ax(5),'A')
xlabel(ax(5),'measurement number')

plot(ax(6),Container{i}(:,1),Container{i}(:,4),'-o')
title(ax(6),'Power')
ylabel(ax(6),'W')
xlabel(ax(6),'measurement number')
fig = gcf;
saveas(fig,[num2str(i),'№_',char(Name1(2,i)),'_Figure.fig'])
end
end
 %% Comparison of temperature and energy flow from a source (200:300)nm, 
figure ('Name',[char(Name1(2,Experiment_number)),'  Date',char(Name1(1,Experiment_number))])
hold 'on'
grid 'on'
plot(Container{Experiment_number}(:,6),Container{Experiment_number}(:,7),'o')
title('Energy flow from a source (200:300)nm')
ylabel('Energy flow (W)')
xlabel('Temperature C ')
a = Container{Experiment_number}(:,1); b = num2str(a); c = cellstr(b);
text(Container{Experiment_number}(:,6), Container{Experiment_number}(:,7), c)
fig = gcf;
saveas(fig,[num2str(Experiment_number),'№_',char(Name1(2,Experiment_number)),'_Energy flow&Temperature_(200&300)nm_Figure.fig.fig'])
% Comparison of temperature and all energy flow from a source (200:800)nm, I - const
figure ('Name',[char(Name1(2,Experiment_number)),'  Date',char(Name1(1,Experiment_number))])
hold 'on'
grid 'on'
plot(Container{Experiment_number}(:,6),Container{Experiment_number}(:,8),'o')
title('All energy flow from a source (200:800)nm')
ylabel('Energy flow (W)')
xlabel('Temperature C ')
a = Container{Experiment_number}(:,1); b = num2str(a); c = cellstr(b);
text(Container{Experiment_number}(:,6), Container{Experiment_number}(:,8), c)
fig = gcf;
saveas(fig,[num2str(Experiment_number),'№_',char(Name1(2,Experiment_number)),'_Energy flow&Temperature_(200&800)nm_Figure.fig'])
%% Comparison of current and energy flow from a source (200:300)nm, 
figure ('Name',[char(Name1(2,Experiment_number)),'  Date',char(Name1(1,Experiment_number))])
hold 'on'
grid 'on'
plot(Container{Experiment_number}(:,2),Container{Experiment_number}(:,7),'o')
title('Energy flow from a source (200:300)nm')
ylabel('Energy flow (W)')
xlabel('I current A ')
a = Container{Experiment_number}(:,1); b = num2str(a); c = cellstr(b);
text(Container{Experiment_number}(:,2), Container{Experiment_number}(:,7), c)
fig = gcf;
saveas(fig,[num2str(Experiment_number),'№_',char(Name1(2,Experiment_number)),'_Energy flow&I current_(200&300)nm_','Figure.fig'])
% Comparison of current and all energy flow from a source (200:800)nm, T - const
figure ('Name',[char(Name1(2,Experiment_number)),'  Date',char(Name1(1,Experiment_number))])
hold 'on'
grid 'on'
plot(Container{Experiment_number}(:,2),Container{Experiment_number}(:,8),'o')
title('All energy flow from a source (200:800)nm')
ylabel('Energy flow (W)')
xlabel('I current A ')
a = Container{Experiment_number}(:,1); b = num2str(a); c = cellstr(b);
text(Container{Experiment_number}(:,2), Container{Experiment_number}(:,8), c)
fig = gcf;
saveas(fig,[num2str(Experiment_number),'№_',char(Name1(2,Experiment_number)),'_Energy flow&I current_(200&800)nm_','Figure.fig'])
%% Power comparison
figure ('Name',[char(Name1(2,Experiment_number)),'  Date',char(Name1(1,Experiment_number))])
hold 'on'
grid 'on'
plot(Container{Experiment_number}(:,2),Container{Experiment_number}(:,8),'ro')
plot(Container{Experiment_number}(:,2),Container{Experiment_number}(:,7),'bo')
plot(Container{Experiment_number}(:,2),Container{Experiment_number}(:,4),'o')
title('Power comparison')
ylabel('Power(W)')
xlabel('I current A')
a = Container{Experiment_number}(:,1); b = num2str(a); c = cellstr(b);
text(Container{Experiment_number}(:,2), Container{Experiment_number}(:,8), c)
text(Container{Experiment_number}(:,2), Container{Experiment_number}(:,7), c)
text(Container{Experiment_number}(:,2), Container{Experiment_number}(:,4), c)
%fig = gcf;
%saveas(fig,[num2str(Experiment_number),'№_',char(Name1(2,Experiment_number)),'Power comparison'])
%% STD 
% figure ('Name',[char(Name1(2,Experiment_number)),'  Date',char(Name1(1,Experiment_number))])
% hold 'on'
% grid 'on'
% plot(Container{Experiment_number}(:,4),Container{Experiment_number}(:,9),'ro')
% title('STD')
% ylabel('STD')
% xlabel('№')
% a = Container{Experiment_number}(:,1); b = num2str(a); c = cellstr(b);
% text(Container{Experiment_number}(:,4), Container{Experiment_number}(:,9), c)

% Osram = ExpDataTo([12 13 14 15 16]);
% DRSH2m = Container{11}(:,[2,6,7]);
% DRSHMitino = ExpDataTo([8 9 10]);
% DRKS = ExpDataTo([6 7]);
% DRSH40416 = ExpDataTo([4 5]);
% DRSH41104 = ExpDataTo([1 2 3]);
% T = (min(Osram(:,2)):25:max(Osram(:,2)));
% I=(min(Osram(:,1)):0.25:max(Osram(:,1)));

% Just for OSRAM cooling comparison 
figure
hold 'on'
grid 'on'
plot(Container{13}(:,2),Container{13}(:,7),'om')
text(Container{13}(:,2),Container{13}(:,7),'0,V')
plot(Container{14}(:,2),Container{14}(:,7),'ob')
text(Container{14}(:,2),Container{14}(:,7),'14,V')
plot(Container{15}(:,2),Container{15}(:,7),'or')
text(Container{15}(:,2),Container{15}(:,7),'28,V')

title('Cooling comparison 200-300')
ylabel('Power UV 200-300(W)')
xlabel('I A')

