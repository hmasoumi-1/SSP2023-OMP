function prob_plot(eta,Prob_exact_supp_eqnorm,Prob_exact_supp1,Prob_exact_supp2)
figure
box on
plot(eta,Prob_exact_supp_eqnorm,'r-s','LineWidth',1.8,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(eta,Prob_exact_supp1,'r-','LineWidth',1.8)
plot(eta,Prob_exact_supp2,'r--','LineWidth',1.8)
% xlim([0.0099,0.66]*sqrt(2*log(32)));
xlim([0.0099,0.66]);
grid on
set(gca, 'FontName', 'Times new roman','FontSize',18)
legend('$d_{\mathrm{max}}/ d_{\mathrm{min}} = 1$','$d_{\mathrm{max}}/ d_{\mathrm{min}} \approx 2.43$',...
    '$d_{\mathrm{max}}/ d_{\mathrm{min}} \approx 2.83$','Interpreter','latex','FontSize', 18)
xlabel('$\frac{\mathbf{x}_{\mathrm{min}}}{\sigma\sqrt{2\log(N)}}$','Interpreter','latex', 'FontSize', 22)
ylabel('Probability of exact support recovery','FontSize', 18)
% export_fig(gcf,'FontMode','fixed','FontSize',10,'Color','abc1.pdf')

end