hxlb=-0.7;
hxub=0.7;
hylb=0; 
hyub=2500;

sxlb=-0.7;
sxub=0.7;
sylb=-0.7; 
syub=0.7;

figure(10)
%1
ax1 = subplot(size(pl_Wind,2),size(pl_Wind,2),1);
histogram(pl_Wind(:,1))
xlim(ax1,[hxlb hxub])
ylim(ax1,[hylb hyub])
ax2 = subplot(size(pl_Wind,2),size(pl_Wind,2),2);
scatter(pl_Wind(:,1),pl_Wind(:,2))
xlim(ax2,[sxlb sxub])
ylim(ax2,[sylb syub])
ax3 = subplot(size(pl_Wind,2),size(pl_Wind,2),3);
scatter(pl_Wind(:,1),pl_Wind(:,3))
xlim(ax3,[sxlb sxub])
ylim(ax3,[sylb syub])
ax4 = subplot(size(pl_Wind,2),size(pl_Wind,2),4);
scatter(pl_Wind(:,1),pl_Wind(:,4))
xlim(ax4,[sxlb sxub])
ylim(ax4,[sylb syub])
ax5 = subplot(size(pl_Wind,2),size(pl_Wind,2),5);
scatter(pl_Wind(:,1),pl_Wind(:,5))
xlim(ax5,[sxlb sxub])
ylim(ax5,[sylb syub])
ax6 = subplot(size(pl_Wind,2),size(pl_Wind,2),6);
scatter(pl_Wind(:,1),pl_Wind(:,6))
xlim(ax6,[sxlb sxub])
ylim(ax6,[sylb syub])
%2
ax7 = subplot(size(pl_Wind,2),size(pl_Wind,2),7);
scatter(pl_Wind(:,2),pl_Wind(:,1))
xlim(ax7,[sxlb sxub])
ylim(ax7,[sylb syub])
ax8 = subplot(size(pl_Wind,2),size(pl_Wind,2),8);
histogram(pl_Wind(:,2))
xlim(ax8,[hxlb hxub])
ylim(ax8,[hylb hyub])
ax9 = subplot(size(pl_Wind,2),size(pl_Wind,2),9);
scatter(pl_Wind(:,2),pl_Wind(:,3))
xlim(ax9,[sxlb sxub])
ylim(ax9,[sylb syub])
ax10 = subplot(size(pl_Wind,2),size(pl_Wind,2),10);
scatter(pl_Wind(:,2),pl_Wind(:,4))
xlim(ax10,[sxlb sxub])
ylim(ax10,[sylb syub])
ax11 = subplot(size(pl_Wind,2),size(pl_Wind,2),11);
scatter(pl_Wind(:,2),pl_Wind(:,5))
xlim(ax11,[sxlb sxub])
ylim(ax11,[sylb syub])
ax12 = subplot(size(pl_Wind,2),size(pl_Wind,2),12);
scatter(pl_Wind(:,2),pl_Wind(:,6))
xlim(ax12,[sxlb sxub])
ylim(ax12,[sylb syub])
%3
ax13 = subplot(size(pl_Wind,2),size(pl_Wind,2),13);
scatter(pl_Wind(:,3),pl_Wind(:,1))
xlim(ax13,[sxlb sxub])
ylim(ax13,[sylb syub])
ax14 = subplot(size(pl_Wind,2),size(pl_Wind,2),14);
scatter(pl_Wind(:,3),pl_Wind(:,2))
xlim(ax14,[sxlb sxub])
ylim(ax14,[sylb syub])
ax15 = subplot(size(pl_Wind,2),size(pl_Wind,2),15);
histogram(pl_Wind(:,3))
xlim(ax15,[hxlb hxub])
ylim(ax15,[hylb hyub])
ax16 = subplot(size(pl_Wind,2),size(pl_Wind,2),16);
scatter(pl_Wind(:,3),pl_Wind(:,4))
xlim(ax16,[sxlb sxub])
ylim(ax16,[sylb syub])
ax17 = subplot(size(pl_Wind,2),size(pl_Wind,2),17);
scatter(pl_Wind(:,3),pl_Wind(:,5))
xlim(ax17,[sxlb sxub])
ylim(ax17,[sylb syub])
ax18 = subplot(size(pl_Wind,2),size(pl_Wind,2),18);
scatter(pl_Wind(:,3),pl_Wind(:,6))
xlim(ax18,[sxlb sxub])
ylim(ax18,[sylb syub])
%4
ax19 = subplot(size(pl_Wind,2),size(pl_Wind,2),19);
scatter(pl_Wind(:,4),pl_Wind(:,1))
xlim(ax19,[sxlb sxub])
ylim(ax19,[sylb syub])
ax20 = subplot(size(pl_Wind,2),size(pl_Wind,2),20);
scatter(pl_Wind(:,4),pl_Wind(:,2))
xlim(ax20,[sxlb sxub])
ylim(ax20,[sylb syub])
ax21 = subplot(size(pl_Wind,2),size(pl_Wind,2),21);
scatter(pl_Wind(:,4),pl_Wind(:,3))
xlim(ax21,[sxlb sxub])
ylim(ax21,[sylb syub])
ax22 = subplot(size(pl_Wind,2),size(pl_Wind,2),22);
histogram(pl_Wind(:,4))
xlim(ax22,[hxlb hxub])
ylim(ax22,[hylb hyub])
ax23 = subplot(size(pl_Wind,2),size(pl_Wind,2),23);
scatter(pl_Wind(:,4),pl_Wind(:,5))
xlim(ax23,[sxlb sxub])
ylim(ax23,[sylb syub])
ax24 = subplot(size(pl_Wind,2),size(pl_Wind,2),24);
scatter(pl_Wind(:,4),pl_Wind(:,6))
xlim(ax24,[sxlb sxub])
ylim(ax24,[sylb syub])
%5
ax25 = subplot(size(pl_Wind,2),size(pl_Wind,2),25);
scatter(pl_Wind(:,5),pl_Wind(:,1))
xlim(ax25,[sxlb sxub])
ylim(ax25,[sylb syub])
ax26 = subplot(size(pl_Wind,2),size(pl_Wind,2),26);
scatter(pl_Wind(:,5),pl_Wind(:,2))
xlim(ax26,[sxlb sxub])
ylim(ax26,[sylb syub])
ax27 = subplot(size(pl_Wind,2),size(pl_Wind,2),27);
scatter(pl_Wind(:,5),pl_Wind(:,3))
xlim(ax27,[sxlb sxub])
ylim(ax27,[sylb syub])
ax28 = subplot(size(pl_Wind,2),size(pl_Wind,2),28);
scatter(pl_Wind(:,5),pl_Wind(:,4))
xlim(ax28,[sxlb sxub])
ylim(ax28,[sylb syub])
ax29 = subplot(size(pl_Wind,2),size(pl_Wind,2),29);
histogram(pl_Wind(:,5))
xlim(ax29,[hxlb hxub])
ylim(ax29,[hylb hyub])
ax30 = subplot(size(pl_Wind,2),size(pl_Wind,2),30);
scatter(pl_Wind(:,5),pl_Wind(:,6))
xlim(ax30,[sxlb sxub])
ylim(ax30,[sylb syub])
%6
ax31 = subplot(size(pl_Wind,2),size(pl_Wind,2),31);
scatter(pl_Wind(:,6),pl_Wind(:,1))
xlim(ax31,[sxlb sxub])
ylim(ax31,[sylb syub])
ax32 = subplot(size(pl_Wind,2),size(pl_Wind,2),32);
scatter(pl_Wind(:,6),pl_Wind(:,2))
xlim(ax32,[sxlb sxub])
ylim(ax32,[sylb syub])
ax33 = subplot(size(pl_Wind,2),size(pl_Wind,2),33);
scatter(pl_Wind(:,6),pl_Wind(:,3))
xlim(ax33,[sxlb sxub])
ylim(ax33,[sylb syub])
ax34 = subplot(size(pl_Wind,2),size(pl_Wind,2),34);
scatter(pl_Wind(:,6),pl_Wind(:,4))
xlim(ax34,[sxlb sxub])
ylim(ax34,[sylb syub])
ax35 = subplot(size(pl_Wind,2),size(pl_Wind,2),35);
scatter(pl_Wind(:,6),pl_Wind(:,5))
xlim(ax35,[sxlb sxub])
ylim(ax35,[sylb syub])
ax36 = subplot(size(pl_Wind,2),size(pl_Wind,2),36);
histogram(pl_Wind(:,6))
xlim(ax36,[hxlb hxub])
ylim(ax36,[hylb hyub])

set(findall(gcf,'-property','FontSize'),'FontSize',20)