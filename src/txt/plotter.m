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

cyclesE1=linspace(1,length(E1),length(E1))*2500;
cyclesE24=linspace(1,length(E24),length(E24))*2500;

figure(3)
E1=E1/400;
plot(cyclesE1,E1,'linewidth',2)
title('T=1.0', 'fontsize',28)
xlabel('Number of cycles', 'fontsize',28)
ylabel('E','fontsize',28)
ylim([-2.1,-0.6])
grid on

figure(4)
E24=E24/400;
plot(cyclesE24,E24,'linewidth',2)
title('T=2,4','fontsize',28)
xlabel('Number of cycles','fontsize',28)
ylabel('E','fontsize',28)
 
grid on





%%
figure(5)
M1=M1/400;
M1=abs(M1);
plot(cyclesE1,M1,'linewidth',2)
title('T=1.0','fontsize',28)
xlabel('Number of cycles','fontsize',28)
ylabel('M','fontsize',28)
ylim([0.2,1.1])
grid on


figure(6)
M24=M24/400;
M24=abs(M24);
plot(cyclesE1,M24,'linewidth',2);
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























