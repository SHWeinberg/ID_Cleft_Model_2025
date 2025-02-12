

%%
s = 'k-';
icleft = iEC(1:end-1);
iintra = setdiff(1:Nnodes-1,icleft);
phi_i = phi_mat(iintra,:);
phi_cleft = phi_mat(icleft,:);

% collect Vm / ionic currents
Vm = phi_mat(Iind(:,1),:) - phi_mat(Iind(:,2),:);

[~,ind] = sort(Iind(:,1));
Iind_cable = Iind(ind,:);
Vm_cable = Vm(ind,:);
tup = tup(ind,:); trepol = trepol(ind,:);

INa_all = I_all(p.iina:Ncurrents:end,:);
INa_all = INa_all(ind,:); % sort to match Vm order

% sort axial / pre-/post-junctional membrane indices to match Vm order
[ind_axial, ~] = ind2sub([length(ind) length(indices.ind_axial)], find(ind == indices.ind_axial));
[ind_disc_pre, ~] = ind2sub([length(ind) length(indices.ind_disc_pre)], find(ind == indices.ind_disc_pre));
[ind_disc_post, ~] = ind2sub([length(ind) length(indices.ind_disc_post)], find(ind == indices.ind_disc_post));
INa_axial = INa_all(ind_axial,:);
INa_disc_pre = INa_all(ind_disc_pre,:);
INa_disc_post = INa_all(ind_disc_post,:);
% Note there is only 1 post-junctional membrane for cell 1 and 
% 1 pre-junctional membrane for cell N for all values of Mdisc, 
% which is coupled directly to the bulk
% To exclude this junction, remove first / last membrane for post / pre

% indices for a specific junction (pre and post)
Njunc = 1;
ind_post = ind_disc_post((Njunc-1)*Mdisc+2:Njunc*Mdisc+1);
Vm_post = Vm_cable(ind_post,:); INa_post = INa_all(ind_post,:);
ind_pre = ind_disc_pre((Njunc-1)*Mdisc+1:Njunc*Mdisc);
Vm_pre = Vm_cable(ind_pre,:); INa_pre = INa_all(ind_pre,:);

Na_cleft_all = S_mat(iEC,:);
Na_cleft = Na_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);       
K_cleft_all = S_mat(iEC+Nnodes,:);
K_cleft = K_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);     
Ca_cleft_all = S_mat(iEC+2*Nnodes,:);
Ca_cleft = Ca_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:); 
A_cleft_all = S_mat(iEC+3*Nnodes,:);
A_cleft = A_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);




%             figure(100); subplot(2,1,1);F
%             plot(ts, Vm_post,'k',ts, Vm_post(ind_conn,:),'r--'); title('single GJ');
% 



tup_axial = tup(ind_axial(1:Nint:end),:);
i1 = round(.25*Ncell); i2 = round(.75*Ncell);
cv_est = 100*(i2-i1)*(p.L/1000)./(tup_axial(i2,:)-tup_axial(i1,:));  % mm/ms = m/s, 100*m/s = cm/s



phi_i = phi_mat(iintra,:);
phi_axial = phi_i(ind_axial,:); 
% 
%         
figure
plot(ts_save,Ca_cleft(:,1:end-1),'color',[1 0 0 0.5])
ylim([1,3])
%         
figure
plot(ts_save,Na_cleft(:,1:end-1),'color',[0 0 1 0.5])
ylim([120,150])



% 
% figure
% plot(ts_save,K_cleft(:,1:end-1),'color',[0 1 0 0.5])
% ylim([4.5,8])
% 
% figure
% plot(ts_save,A_cleft(:,1:end-1),'color',[0 0 0 0.5])
% ylim([140,150])
% 
% figure
% plot(ts_save,phi_axial(5,2:end))
% 
% figure;
% h = pcolor(ts, 1:length(iintra), phi_mat(iintra,:));
% set(gca,'ydir','reverse');
% h.LineStyle = 'none';
% set(gca,'ytick',1:length(iintra), 'yticklabel',iintra);
% colorbar
% caxis([-80,0])
% 
% subplot(2,2,3); plot(ts_save, phi_i,s); hold on;
% 
% subplot(2,2,2);
% h = pcolor(ts, 1:length(icleft), phi_mat(icleft,:));
% set(gca,'ydir','reverse');
% h.LineStyle = 'none';
% set(gca,'ytick',1:20:length(icleft), 'yticklabel',icleft(1:20:end));
% 
% subplot(2,2,4); plot(ts, phi_cleft, s); hold on;
% 
% figure
% plot(INa_all(ind_pre,:)'.*1e6)
% ylim([-140,20])
% %


%% analysis
  
list_higher_post = [];
list_higher_pre = [];
for i = 1:100
    
    if min(INa_pre(i,1:end)) > min(INa_post(i,1:end))
        list_higher_post = [list_higher_post, i];
    else
        list_higher_pre = [list_higher_pre, i];
    end
end

figure
plot3(FEM_data.partition_centers(list_higher_post, 1), ...
      FEM_data.partition_centers(list_higher_post, 2), ...
      FEM_data.partition_centers(list_higher_post, 3), '.','color','red', 'markersize',10)
hold on
plot3(FEM_data.partition_centers(list_higher_pre, 1), ...
      FEM_data.partition_centers(list_higher_pre, 2), ...
      FEM_data.partition_centers(list_higher_pre, 3), '.','color','blue', 'markersize',10)

  
  
figure
plot(ts_save,INa_pre(list_higher_post,1:end),'color',[0 0 1 0.5])
hold on
plot(ts_save,INa_post(list_higher_post,1:end),'color',[0 1 1 0.5])
xlim([400.8,401.6])

figure
plot(ts_save,Na_cleft(list_higher_post,1:end-1),'color',[0 0 1 0.5])
hold on
plot(ts_save,Na_cleft(list_higher_pre,1:end-1),'color',[0 1 1 0.5])

xlim([400.8,405.6])

