function nmse_plot(eta,NMSE_eqnorm,NMSE1,NMSE2,UpperBound_eqnorm,UpperBound1,UpperBound2)
figure
box on
plot(eta,10*log10(NMSE_eqnorm),'r-s','LineWidth',1.8,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(eta,10*log10(NMSE1),'r-','LineWidth',1.8)
plot(eta,10*log10(NMSE2),'r--','LineWidth',1.8)
plot(eta,10*log10(UpperBound_eqnorm),'b:s','LineWidth',1.8,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','y')
plot(eta,10*log10(UpperBound1),'b:','LineWidth',1.8)
plot(eta,10*log10(UpperBound2),'b-.','LineWidth',1.8)
xlim([5.99,24.001]);
ylim([-60,12]);
grid on
set(gca, 'FontName', 'Times new roman','FontSize',18)
legend({'$d_{\mathrm{max}}/ d_{\mathrm{min}} = 1$','$d_{\mathrm{max}}/ d_{\mathrm{min}} \approx 2.43$',...
    '$d_{\mathrm{max}}/ d_{\mathrm{min}} \approx 2.83$','$d_{\mathrm{max}}/ d_{\mathrm{min}} \approx 1$',...
    '$d_{\mathrm{max}}/ d_{\mathrm{min}} \approx 2.43$','$d_{\mathrm{max}}/ d_{\mathrm{min}} \approx 2.83$',...
    },'Location','northeast','NumColumns',2,'Interpreter','latex','FontSize', 18)
legend show
xlabel('$\frac{\mathbf{x}_{\mathrm{min}}}{\sigma\sqrt{2\log(N)}}$','Interpreter','latex','FontSize', 22)
ylabel('Normalized MSE [dB]','FontSize', 18)

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.138,0.162,0.23,0.125])
box(ax,'on')
plot(eta,10*log10(NMSE_eqnorm),'r-s','LineWidth',1.8,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','y','parent',ax)
hold on
plot(eta,10*log10(NMSE1),'r-','LineWidth',1.8,'parent',ax)
plot(eta,10*log10(NMSE2),'r--','LineWidth',1.8,'parent',ax)
% set(ax,'xlim',[9.8,10.2],'ylim',[-44.75,-43])
set(ax,'xlim',[9.8,10.2],'ylim',[-44.5,-42.5])
set(gca, 'FontName', 'Times new roman','FontSize',16)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

% -- Annotations --
% dim_rect = [.28 .3073 .054 .0324];
dim_rect = [.28 .313 .054 .0324];
annotation('rectangle',dim_rect,'FaceColor','blue','FaceAlpha',.2)

aa = annotation('line',[.13757 .28],[.287 0.3454]);
aa.LineWidth = 0.6;
aa.LineStyle = '-';

bb = annotation('line',[.368 .334],[.287 0.3454]);
bb.LineWidth = 0.6;
bb.LineStyle = '-';

% --
dim_ell_bound = [.192 .433 .05 .361];
annotation('ellipse',dim_ell_bound)

dim = [.3 .433 .4 .1];
str = 'Upper-bounds from Theorem 2';
cc = annotation('textbox',dim,'String',str);
cc.FontSize = 16;
cc.FontName = 'Times new roman';
cc.EdgeColor = 'none';

dd = annotation('arrow',[.242 .31],[.608 0.523]);
dd.LineWidth = 0.6;
dd.LineStyle = '-';

% --
dim_ell_emp = [.4655 .254 .02 .07];
annotation('ellipse',dim_ell_emp)

dim = [.53 .255 .3 .1];
str = 'Empirical NMSE';
ee = annotation('textbox',dim,'String',str);
ee.FontSize = 16;
ee.FontName = 'Times new roman';
ee.EdgeColor = 'none';

ff = annotation('arrow',[.483 .535],[.31 0.325]);
ff.LineWidth = 0.6;
ff.LineStyle = '-';
% export_fig(gcf,'FontMode','fixed','FontSize',10,'Color','nmse.pdf')
end
