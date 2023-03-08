%%
% First run parts of main_pfc.m and main_tci.m and save results
% in out_pfc and out_tci for comparison
clear all;close all;
Q_rand = [];
load('Q_rand_20pct.mat')
main_pfc;
main_tci;
%%

CEP_REF_plot = CEP_REF*ones(length(out.time), 1);
BIS_REF_plot = ContentrationToBIS(CEP_REF_plot./C50P, CER_REF./C50R, gamma, sigma, Emax);

% Divided by 60 due to conversion into minutes
figure;
subplot(3,1,1);
plot(out_pfc.time./60, BIS_REF_plot, 'k');
%title("Primerjava PFC in TCI", 'Interpreter', 'latex');
hold on;
plot(out_pfc.time./60, out_pfc.BIS_combined, 'r');
plot(out_pfc.time./60, out_tci.BIS_combined, 'b');
legend('ref','PFC', 'TCI', 'Interpreter', 'latex');
xlabel("$t [min]$", 'Interpreter', 'latex');
ylabel("$BIS [\%]$", 'Interpreter', 'latex');
ylim([0,Emax*1.05]);

subplot(3,1,2);
plot(out_pfc.time./60, CEP_REF_plot, 'k');
hold on;
plot(out_pfc.time./60, out_pfc.CePrs, 'r');
plot(out_pfc.time./60, out_tci.CePrs, 'b');
legend('ref','PFC', 'TCI', 'Interpreter', 'latex', 'Location','southeast');
xlabel("$t [min]$", 'Interpreter', 'latex');
ylabel("$x_e [mg/ml]$", 'Interpreter', 'latex');
ylim([0, max(max(out_tci.CePrs), max(out_pfc.CePrs))*1.05]);

subplot(3,1,3);
plot(out_pfc.time./60, out_pfc.u_reg_sat, 'r');
hold on;
plot(out_pfc.time./60, out_tci.u_reg_sat, 'b');
legend('PFC', 'TCI', 'Interpreter', 'latex');
xlabel("$t [min]$", 'Interpreter', 'latex');
ylabel("$u_{Prop} [mg/ml/min]$", 'Interpreter', 'latex');
ylim([-0.5, 41]);

%%
cumulative_I_PFC = cumtrapz(out_pfc.time./60,out_pfc.u_reg_sat);
dosage_PFC = cumulative_I_PFC(end);

cumulative_I_TCI = cumtrapz(out_tci.time./60,out_tci.u_reg_sat);
dosage_TCI = cumulative_I_TCI(end);

if dosage_PFC < dosage_TCI
    msg_winner = 'PFC';
    dosage_winner = dosage_PFC;
    dosage_looser = dosage_TCI;
else
    msg_winner = 'TCI';
    dosage_winner = dosage_TCI;
    dosage_looser = dosage_PFC;
end
pct_better = (dosage_looser - dosage_winner)/dosage_looser*100;

fprintf('\n');
fprintf('*******************START INFO*******************\n');
fprintf('Propofol dosage -- PFC: %.2f, TCI: %.2f\n', dosage_PFC, dosage_TCI);
fprintf('%.2f %% less dosage using %s.\n', pct_better, msg_winner);
fprintf('********************END INFO********************\n');
fprintf('\n');
