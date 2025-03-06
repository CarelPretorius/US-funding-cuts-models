clear all; 
load calibres1.mat;

xsto = xsto3;

ix0 = 2e4; nx = 150; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);
pct_xs = prctile(xs,[2.5,50,97.5])

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    [out, aux] = obj(xs(ii,:));
    sims(ii,:) = [aux.incd(end),aux.noti(end), aux.mort(end),aux.mdr(end-4), aux.mdr(end), aux.mdriniTX(end), aux.sym];
    incd(:,ii) = aux.incd;  
    mdr(:,ii)  = aux.mdr; 
end
fprintf('\n');

sim_pct = prctile(sims',[2.5,50,97.5],2);
incd_pct = prctile(incd,[2.5,50,97.5],2);  %1997-2023
mdr_pct = prctile(mdr,[2.5,50,97.5],2);    %1997-2023

% --- Show all on a plot --------------------------------------------------
f1 = figure; 
dx = 0.1; ms = 24; fs = 14;
sizes1 = [28, 28, 28];
sizes2 = [28, 28, 28, 28];

% --- Population rates, i.e. prevalence and notifications 
subplot(1,2,1)
sim  = sim_pct(1:3,:);
dat  = [data.inc2023; data.noti2023; data.mort]';  

xx   = [1:3]-dx;
plt  = dat';
hilo = diff(plt,[],2);
lg(1,:) = scatter(xx, plt(:,2),sizes1,'filled','r'); hold on
errorbar(xx, plt(:,2), hilo(:,1), hilo(:,2), 'Linestyle', 'None','Color', 'r'); hold on

xx   = [1:3]+dx;
plt  = sim(1:3,:);
hilo = diff(plt,[],2);
lg(2,:) = scatter(xx, plt(:,2),sizes1,'filled','b');
errorbar(xx, plt(:,2), hilo(:,1), hilo(:,2), 'Linestyle', 'None', 'Color', 'b');

set(gca,'fontsize',fs,'XTick',1:3,'XTicklabel',{'Incidence','Notification','Mortality'});

xtickangle(30);
title('Country2')
ylabel('Rate per 100,000','fontsize',fs);
yl = ylim; yl(1) = 0; ylim(yl);
xlim([0.5 3.5]);
legend(lg, 'Data','Model','Location','NorthEast');



subplot(1,2,2)
sim  = sim_pct(4:end,:);
dat  = [data.mdr2019; data.mdr2023; data.mdriniTX; data.sym]';  

xx   = [1:4]-dx;
plt  = dat';
hilo = diff(plt,[],2);
lg(1,:) = scatter(xx, plt(:,2),sizes2,'filled','r'); hold on
errorbar(xx, plt(:,2), hilo(:,1), hilo(:,2), 'Linestyle', 'None','Color', 'r'); hold on

xx   = [1:4]+dx;
plt  = sim;
hilo = diff(plt,[],2);
lg(2,:) = scatter(xx, plt(:,2),sizes2,'filled','b');
errorbar(xx, plt(:,2), hilo(:,1), hilo(:,2), 'Linestyle', 'None', 'Color', 'b');

set(gca,'fontsize',fs,'XTick',1:4,'XTicklabel',{'MDR Incidence 2019','MDR Incidence 2023','2nd line TX initiation','Proportion of symptomatic'});

xtickangle(30);
title('Country2')
ylabel('Rate per 100,000','fontsize',fs);
yl = ylim; yl(1) = 0; ylim(yl);
xlim([0.5 4.5]);
%title('Population rates');
legend(lg, 'Data','Model','Location','NorthEast');






%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
years = 1997:2023;
subplot(1,2,1)
plt = incd_pct'; 
lg(1,:) = plot(years, plt(2,:),'b');    % Incidence
jbfill(years, plt(3,:), plt(1,:), 'b', 'None', 1, 0.3); hold on

ms = 24; fs = 14;
dat  = [data.inc2023];
plt2  = dat(1,:);
hilo = diff(plt2,[],2);
lg(2,:) = plot(2023, plt2(:,2), 'm.', 'markersize', ms); hold on
errorbar(2023, plt2(:,2), hilo(:,1), hilo(:,2),'m', 'Linestyle', 'None','LineWidth',1.5);
xlim([2011 2023]);

set(gca,'fontsize',fs);
ylabel('Rate per 100,000','fontsize',fs);
title('Incidence')

subplot(1,2,2)
plt = mdr_pct';
lg(3,:) = plot(years, plt(2,:),'r');  % MDR incidence
jbfill(years, plt(3,:), plt(1,:), 'r', 'None', 1, 0.3); hold on

dat  = [data.mdr2019];
plt2  = dat(1,:);
hilo = diff(plt2,[],2);
lg(4,:) = plot(2015, plt2(:,2), 'm.', 'markersize', ms); hold on
errorbar(2015, plt2(:,2), hilo(:,1), hilo(:,2),'m', 'Linestyle', 'None','LineWidth',1.5);

dat  = [data.mdr2023];
plt2  = dat(1,:);
hilo = diff(plt2,[],2);
lg(5,:) = plot(2023, plt2(:,2), 'm.', 'markersize', ms); hold on
errorbar(2023, plt2(:,2), hilo(:,1), hilo(:,2),'m', 'Linestyle', 'None','LineWidth',1.5);

set(gca,'fontsize',fs);
ylabel('Rate per 100,000','fontsize',fs);
xlim([2011 2023]);
title('MDR-Incidence')

% Save data for plotting calibration result 
caldat.model = [sim_pct(1:end-1,:); sim_pct(end,:)*100];
caldat.datain = [data.inc2023;data.noti2023;data.mort;data.mdr2019;data.mdr2023;data.mdriniTX;data.sym*100];
save ('country2_cal','caldat')
