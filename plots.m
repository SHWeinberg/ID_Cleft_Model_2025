clear
addpath(genpath('/users/PAS1622/nickmoise/.matlab/down'))

%%
%get data from juncs, axial
load("/fs/scratch/PAS1622/nickmoise/ID_2025/gj_chan_loc_base_D1/" + ...
     "msh1_baseline_cycle1000_beats10_D1_gj_loc8_chan_loc8.mat")


Njunc = 25;

icleft = iEC(1:end-1);
iintra = setdiff(1:Nnodes-1,icleft);
phi_i = phi_save(iintra,:);
phi_cleft = phi_save(icleft,:);
% collect Vm / ionic currents
Vm = phi_save(Iind(:,1),:) - phi_save(Iind(:,2),:);

[~,ind] = sort(Iind(:,1));
Iind_cable = Iind(ind,:);
Vm_cable = Vm(ind,:);
tup = tup(ind,:); trepol = trepol(ind,:);

INa_all = I_save(p.iina:Ncurrents:end,:);

% sort axial / pre-/post-junctional membrane indices to match Vm order
[ind_axial, ~] = ind2sub([length(ind) length(indices.ind_axial)], find(ind == indices.ind_axial));
[ind_disc_pre, ~] = ind2sub([length(ind) length(indices.ind_disc_pre)], find(ind == indices.ind_disc_pre));
[ind_disc_post, ~] = ind2sub([length(ind) length(indices.ind_disc_post)], find(ind == indices.ind_disc_post));
INa_axial = INa_all(ind_axial,:);
INa_disc_pre = INa_all(ind_disc_pre,:);
INa_disc_post = INa_all(ind_disc_post,:);

ind_post = ind_disc_post((Njunc-1)*Mdisc+2:Njunc*Mdisc+1);
Vm_post = Vm_cable(ind_post,:); INa_post = INa_all(ind_post,:);
ind_pre = ind_disc_pre((Njunc-1)*Mdisc+1:Njunc*Mdisc);
Vm_pre = Vm_cable(ind_pre,:); INa_pre = INa_all(ind_pre,:);

Na_cleft_all = S_save(iEC,:);
Na_cleft = Na_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);       
K_cleft_all = S_save(iEC+Nnodes,:);
K_cleft = K_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);     
Ca_cleft_all = S_save(iEC+2*Nnodes,:);
Ca_cleft = Ca_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:); 
A_cleft_all = S_save(iEC+3*Nnodes,:);
A_cleft = A_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);





tup_axial = tup(ind_axial(1:Nint:end),:);
i1 = round(.25*Ncell); i2 = round(.75*Ncell);
cv_est = 100*(i2-i1)*(p.L/1000)./(tup_axial(i2,:)-tup_axial(i1,:));  % mm/ms = m/s, 100*m/s = cm/s



phi_i = phi_save(iintra,:);
phi_axial = phi_i(ind_axial,:); 





%
% x = ts_save;
% thresh_activation = -20;
% % local_cv = [];
% for i = 1:50
% 
%     data_iso_int = phi_axial(i,:);
% 
% 
% 
%     data_iso1 = data_iso_int(:,1:end-1);
%     data_iso2 = data_iso_int(:,2:end);
%     ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
%     y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
%     m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
%     activ_time(i) = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
% 
% 
% end
% 
% local_cv = (100./diff(activ_time))./10; %cm/s
% local_cv = local_cv(5:45); %cm/s

% figure
% plot(local_cv)
% % ylim([36,45])
% ylim([12,18])


%% CV resti data
data_dir = "/fs/scratch/PAS1622/nickmoise/ID_2025/" + ...
           "gj_chan_loc_base_D1/";



data_list = dir(data_dir);
data_list = natsortfiles(data_list);
data_list(1:2) = [];

scale_gj_loc_vec = [8 4 2 1 1/2 1/4 1/8];
scale_chan_loc_vec = [8 4 2 1 1/2 1/4 1/8];
CV_mat = zeros(7,7);
[X,Y] = meshgrid(scale_gj_loc_vec,scale_chan_loc_vec);


thresh_activation = 0;
cycle_vec = [200:5:300 350:50:1000];
for i_file = 1:length(data_list)
    
    filename = data_dir + "/" +  data_list(i_file).name;
    load(filename, 'phi_axial_all','ts')
    phi_axial_all(:,end) = phi_axial_all(:,end-1);
    
    ind_last_2beat = find(ts>(cycle_vec(end)*8-10));
    
    phi_axial_all = phi_axial_all(:,ind_last_2beat);
    
    x = ts(ind_last_2beat);
    cell_index = [10,40]; % cells from which CV is comp
    activ_time = zeros(2,20);
    for i_cell = 1:length(cell_index)
        data_iso_int = phi_axial_all(cell_index(i_cell),:);
        data_iso1 = data_iso_int(:,1:end-1);
        data_iso2 = data_iso_int(:,2:end);
        ind_th = find(data_iso1<thresh_activation & data_iso2> thresh_activation); % find indices below/above threshold
        y1 = data_iso1(ind_th); y2 = data_iso2(ind_th);
        m_slope = (y2-y1)./(x(ind_th+1)-x(ind_th));  % linear slope
        activ_time_cell = x(ind_th) - (data_iso_int(ind_th)-thresh_activation)./m_slope;
        activ_time(i_cell,1:length(activ_time_cell)) = activ_time_cell;
    end
    
    beats1 = activ_time(1,:) ~= 0;
    beats1(beats1==0) = [];
    beats2 = activ_time(2,:) ~= 0;
    beats2(beats2==0) = [];
    
    cell_dist = (cell_index(2) - cell_index(1))*100; %100 = cell length in um
    CV_cycles(i_file) = cell_dist./(activ_time(2,length(beats2)) - activ_time(1,length(beats1))) ./10;    
    
    load(filename, 'scale_gj_loc', 'scale_chan_loc')

    CV_mat(find(scale_gj_loc_vec == scale_gj_loc), ...
           find(scale_chan_loc_vec == scale_chan_loc)) = CV_cycles(i_file);


    % figure
    % imagesc(ts,1:50,phi_axial_all)
    % ax = gca;
    % ax.DataAspectRatio = [100 1 1];
    
end

% plot(ts,phi_axial_all(5,:))

% plot(CV_cycles)

figure
surf(log2(X),log2(Y),CV_mat)
%%
% 
% %%
% figure
plot(ts_save,Ca_cleft(:,1:end),'color',[1 0 0 0.5])
ylim([1,3])
% %         
% figure
% plot(ts_save,Na_cleft(:,1:end),'color',[0 0 1 0.5])
% ylim([120,150])
% 
% figure
% plot(ts_save,phi_axial')



% for i = 1:25
%     phi = phi_axial(i,:);
%     phi(phi>-70) = 1;
%     phi(phi<=-70) = 0;
%     [~,peak_time] = findpeaks(phi);
% 
%     peak_list(i) = ts_save(peak_time(end));
% end

% figure
% hold on
% for i = 1:25
%     plot(ts_save,phi_axial(i,1:end))
% end


CV = 0.01./(diff(peak_list)./1000);
%         
%         figure
%         plot(ts_save,K_cleft(:,1:end-1),'color',[0 1 0 0.5])
%         ylim([4.5,8])
%         
%         figure
%         plot(ts_save,A_cleft(:,1:end-1),'color',[0 0 0 0.5])
%         ylim([140,150])
% 
figure
plot(ts_save,phi_axial(5,1:end))
% 
%         figure;
%         h = pcolor(ts, 1:length(iintra), phi_save(iintra,:));
%         set(gca,'ydir','reverse');
%         h.LineStyle = 'none';
%         set(gca,'ytick',1:length(iintra), 'yticklabel',iintra);
%         colorbar
%         caxis([-80,0])
% 
% subplot(2,2,3); plot(ts_save, phi_i,s); hold on;
% 
%         subplot(2,2,2);
%         h = pcolor(ts, 1:length(icleft), phi_save(icleft,:));
%         set(gca,'ydir','reverse');
%         h.LineStyle = 'none';
%         set(gca,'ytick',1:20:length(icleft), 'yticklabel',icleft(1:20:end));
% 
%         subplot(2,2,4); plot(ts, phi_cleft, s); hold on;

%         figure
%         plot(INa_all(ind_pre,:)'.*1e6)
%         ylim([-140,20])
%%
