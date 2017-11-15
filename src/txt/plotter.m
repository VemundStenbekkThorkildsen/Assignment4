%% load files 
clc
close all
clear all
A1=load('acceptanceMCrandom1.000000.txt');
A2=load('acceptanceMCrandom2.400000.txt');

A3=load('acceptanceMCup1.000000.txt');
A4=load('acceptanceMCup2.400000.txt');

a=load('TempAcceptanceRandom.txt');
b=load('TempAcceptanceUp.txt');

E1=load('energy1.000000.txt');
E24=load('energy2.400000.txt');

M1=load('magnetization1.000000.txt');
M24=load('magnetization2.400000.txt');


%% acceptance random for t=1.0 an t=2.4
cycles1=linspace(1,length(A1),length(A1));
subplot(1,2,1)
plot(cycles1,A1)
title('T=1.0');
xlabel('Number of cycles')
ylabel('Accepted configuration')
grid on
subplot(1,2,2)
plot (cycles1,A2)
title('T=2.4')
xlabel('Number of cycles')
ylabel('Accepted configuration')
grid on




%% acceptance ordered for t=1.0 and t=2.4
cycles2=linspace(1,length(A3),length(A3));
subplot(1,2,1)
plot(cycles2,A3);
title('T=1.0');
xlabel('Number of cycles')
ylabel('Accepted configuration')
grid on
subplot(1,2,2)
plot(cycles2,A4)
title('T=2.4');
xlabel('Number of cycles')
ylabel('Accepted configuration')
grid on


%% acceptance vs temperature



figure(1)
plot(a(:,1),a(:,2));
title('Random initial matrix')
xlabel('T')
ylabel('Accepted moves')
grid on


figure(2)
plot(b(:,1),b(:,2));
title('All spins up initial matrix')
xlabel('T')
ylabel('Accepted moves')
grid on 

%% energyplots for t=1.0 and t=2.4

cyclesE1=linspace(1,length(E1),length(E1))*5000;
cyclesE24=linspace(1,length(E24),length(E24))*5000;

figure(3)

plot(cyclesE1,E1,'linewidth',2)
title('T=1.0', 'fontsize',28)
xlabel('Number of cycles', 'fontsize',28)
ylabel('E','fontsize',28)
ylim([-810,-300])
grid on

figure(4)

plot(cyclesE24,E24,'linewidth',2)
title('T=2,4','fontsize',28)
xlabel('Number of cycles','fontsize',28)
ylabel('E','fontsize',28)
 
grid on





%%
figure(5)

%M1=abs(M1);
plot(cyclesE1,M1,'linewidth',2)
title('T=1.0','fontsize',28)
xlabel('Number of cycles','fontsize',28)
ylabel('M','fontsize',28)

grid on


figure(6)

%M24=abs(M24);
plot(cyclesE24,M24,'linewidth',2);
title('T=2.4','fontsize',28);
xlabel('number of cycles','fontsize',28)
ylabel('M','fontsize',28)
grid on


%% number of appeareances per energy 
E1=load('energy1.000000.txt');
E125=load('energy1.250000.txt');
E15=load('energy1.500000.txt');
E175=load('energy1.750000.txt');
E2=load('energy2.000000.txt');
E225=load('energy2.250000.txt');
E25=load('energy2.500000.txt');
E275=load('energy2.750000.txt');
E3=load('energy3.000000.txt');


edges1 = unique(E1);
edges125 = unique(E125);
edges15 = unique(E15);
edges175 = unique(E175);
edges2 = unique(E2);
edges225 = unique(E225);
edges25 = unique(E25);
edges275 = unique(E275);
edges3 = unique(E3);


counts1 = histc(E1(:), edges1);
counts125 = histc(E125(:), edges125);
counts15 = histc(E15(:), edges15);
counts175 = histc(E175(:), edges175);
counts2 = histc(E2(:), edges2);
counts225 = histc(E225(:), edges225);
counts25 = histc(E25(:), edges25);
counts275 = histc(E275(:), edges275);
counts3 = histc(E3(:), edges3);



figure(7)
hold all
plot (edges1,counts1,'linewidth',2);
plot (edges125,counts125,'linewidth',2);
plot (edges15,counts15,'linewidth',2);
plot (edges175,counts175,'linewidth',2);
plot (edges2,counts2,'linewidth',2);
plot (edges225,counts225,'linewidth',2);
plot (edges25,counts25,'linewidth',2);
plot (edges275,counts275,'linewidth',2);
plot (edges3,counts3,'linewidth',2);
title('Random initial matrix for varying temperature','fontsize',28)
l=legend('1.0','1.25','1.5','1.75','2.0','2,25','2,5','2.75','3.0')
l.FontSize=28;
ylabel('number of occurences')
xlabel('Energy')
grid on
xlim([-800,-100])

skwedd=skewness(E225);
Skrfgqw=skewness(275)
%%
E24=load('energy2.400000.txt');
edges24 = unique(E24);
counts24 = histc(E24(:), edges24);


figure(8)
plot (edges24,counts24,'linewidth',2);
title('T=2.4')
xlabel('Energy')
ylabel('number of occurances')

meanE= mean(E24)
medaiE=median(E24)
skew=skewness(E24)


%% Energy RANDOM matrix
E24=load('CumEnergyRAND2.400000.txt');
E1=load('CumEnergyRAND1.000000.txt');
UE24=load('CumEnergyUP2.400000.txt');
UE1=load('CumEnergyUP1.000000.txt');


avg1=E1(:,2)./E1(:,1);
avg24=E24(:,2)./E24(:,1);
Uavg1=UE1(:,2)./UE1(:,1);
Uavg24=UE24(:,2)./UE24(:,1);

figure(9)
hold on
plot(E1(:,1),avg1,'linewidth',2);
plot(UE1(:,1),Uavg1,'linewidth',2);
title('Average energy when T=1.0','fontsize',28)
xlabel('Monte carlo cycles','fontsize',28);
ylabel('Average Energy','fontsize',28)
p=legend('Random Matrix','Up matrix')
ylim([-810,-200])
p.FontSize=28;
grid on

figure(10)
hold on
plot(E24(:,1),avg24,'linewidth',2);
plot(UE24(:,1),Uavg24,'linewidth',2);
title('Average energy when T=2.4','fontsize',28)
xlabel('Monte carlo cycles','fontsize',28);
ylabel('Average Energy','fontsize',28)
l=legend('Random Matrix','Up matrix')
l.FontSize=28;
grid on

%% magnetization RAND and UP matrix
UM24=load('CumMagneUP2.400000.txt');
UM1=load('CumMagneUP1.000000.txt');
M24=load('CumMagneRAND2.400000.txt');
M1=load('CumMagneRAND1.000000.txt');



avgM1=(M1(:,2)./M1(:,1));
avgM24=M24(:,2)./M24(:,1);
UavgM1=(UM1(:,2)./UM1(:,1));
UavgM24=UM24(:,2)./UM24(:,1);


figure(11)
hold on
plot(M1(:,1),avgM1,'linewidth',2);
plot(UM1(:,1),UavgM1,'linewidth',2);
title('Average absolute magnetization when T=1.0','fontsize',28)
p=legend('Random Matrix','Up matrix')
p.FontSize=28;
xlabel('Monte carlo cycles','fontsize',28);
ylabel('Average Absolute magnetization','fontsize',28)
ylim([50,410])
grid on

figure(12)
hold on
plot(M24(:,1),avgM24,'linewidth',2);
plot(UM24(:,1),UavgM24,'linewidth',2);
title('Average absolute magnetization when T=2.4','fontsize',28)
l=legend('Random Matrix','Up matrix')
l.FontSize=28;
xlabel('Monte carlo cycles','fontsize',28);
ylabel('Average Absolute magnetization','fontsize',28)
grid on


%% Magnetization up matrix
UM24=load('CumMagneUP2.400000.txt');
UM1=load('CumMagneUP1.000000.txt');



UavgM1=(UM1(:,2)./UM1(:,1));
UavgM24=UM24(:,2)./UM24(:,1);

figure(15)
plot(UM1(:,1),UavgM1,'linewidth',2);
title('Average absolute magnetization when T=1.0','fontsize',28)
xlabel('Monte carlo cycles','fontsize',28);
ylabel('Absolute magnetizm','fontsize',28)
grid on

figure(16)
plot(UM24(:,1),UavgM24,'linewidth',2);
title('Average absolute magnetization when T=2.4','fontsize',28)
xlabel('Monte carlo cycles','fontsize',28);
ylabel('Absolute magnetism','fontsize',28)
grid on


%% Energy UP matrix


UE1(1,2)=UE1(1,2)/100;
UE24(1,2)=UE24(1,2)/100;

Uavg1=UE1(:,2)./UE1(:,1);
Uavg24=UE24(:,2)./UE24(:,1);

figure(13)
plot(UE1(:,1),Uavg1,'linewidth',2);
title('Average energy when T=1.0','fontsize',28)
xlabel('Monte carlo cycles','fontsize',28);
ylabel('Energy','fontsize',28)
grid on

figure(14)
plot(UE24(:,1),Uavg24,'linewidth',2);
title('Average energy when T=2.4','fontsize',28)
xlabel('Monte carlo cycles','fontsize',28);
ylabel('Energy','fontsize',28)
grid on










