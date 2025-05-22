%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P_act_opt_multi_initial_cond.m ver.250522%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
clear global
close all

%Parameter settings
U_para=2;
cu=0.1;
cd=-1;
mu0=0.38;
mu2=-0.01;
sigma2=0.08;
nu=U_para+2i*cu;
gamma=1+1i*cd;

T=5;

addpath(genpath('DMSUITE'));
n=220;
b=0.239;

T_tilde=5;
p_tilde=3;
z_min_act=-17;
z_max_act=9;
case_num_initial_max=50000;

r_list=[10 5];
p_list=[3 8];
rp_num_max=size(r_list,2);
index_select_manual_list=zeros(max(p_list),rp_num_max);
index_select_manual_list(:,1)=[104; 98; 110; 0; 0; 0; 0; 0];
index_select_manual_list(:,2)=[104; 98; 110; 94; 114; 91; 117; 100];
case_num_act_max=50000;
run_num_max=5;
for rp_num=1:rp_num_max
    p=p_list(1,rp_num);
    for row_num=1:p
        if index_select_manual_list(row_num,rp_num)==0
            fprintf('Number of acuators manually placed does not correspond to p=%d (rp_num=%d)\n',p,rp_num)
            pause
        end
    end
end

sigma_min=10^(-10);

terminal_disp_min_ex=-4;
terminal_disp_min=10^terminal_disp_min_ex;
terminal_disp_max_ex=0;
terminal_disp_max=10^terminal_disp_max_ex;
input_disp_min_ex=-2;
input_disp_min=10^input_disp_min_ex;
input_disp_max_ex=6;
input_disp_max=10^input_disp_max_ex;
prob_dens_disp_min_ex=-6;
prob_dens_disp_min=10^prob_dens_disp_min_ex;
prob_dens_disp_max_hist_ex=4;
prob_dens_disp_max_hist=10^prob_dens_disp_max_hist_ex;
bin_num_max=200;

scrsz=get(0,'ScreenSize');
fig_size=[700 600];
font_name='Times New Roman';
font_size=30;
alpha=0.3;

%Discretization of the model
[z,D]=herdif(n,2,b);
D1=D(:,:,1);
D2=D(:,:,2);

M=diag(mu0-cu^2+mu2*(z.^2)/2);
A=-nu*D1+gamma*D2+M;
B=zeros(n,n);
for row_num=1:n
    for col_num=1:n
        B(row_num,col_num)=1/sqrt(2*pi*sigma2)*exp(-(z(row_num,1)-z(col_num,1))^2/(2*sigma2));%ガウス分布に従う
    end
end
C=zeros(n,n);
C(1,1)=sqrt((z(2,1)-z(1,1))/2);
for row_num=2:n-1
    C(row_num,row_num)=sqrt((z(row_num+1,1)-z(row_num-1,1))/2);
end
C(n,n)=sqrt((z(n,1)-z(n-1,1))/2);

%Preparation regarding bin
bin_edge_hist_terminal=10.^linspace(terminal_disp_min_ex,terminal_disp_max_ex,bin_num_max+1);
bin_edge_hist_input=10.^linspace(input_disp_min_ex,input_disp_max_ex,bin_num_max+1);
bin_edge_contour_terminal=10.^linspace(terminal_disp_min_ex,terminal_disp_max_ex,bin_num_max+1);
bin_edge_contour_input=10.^linspace(input_disp_min_ex,input_disp_max_ex,bin_num_max+1);
bin_width_hist_terminal=zeros(1,bin_num_max);
bin_width_hist_input=zeros(1,bin_num_max);
bin_width_contour_terminal=zeros(1,bin_num_max);
bin_width_contour_input=zeros(1,bin_num_max);
for bin_num=1:bin_num_max
    bin_width_hist_terminal(1,bin_num)=bin_edge_hist_terminal(1,bin_num+1)-bin_edge_hist_terminal(1,bin_num);
    bin_width_hist_input(1,bin_num)=bin_edge_hist_input(1,bin_num+1)-bin_edge_hist_input(1,bin_num);
    bin_width_contour_terminal(1,bin_num)=bin_edge_contour_terminal(1,bin_num+1)-bin_edge_contour_terminal(1,bin_num);
    bin_width_contour_input(1,bin_num)=bin_edge_contour_input(1,bin_num+1)-bin_edge_contour_input(1,bin_num);
end
bin_width_contour_terminal_mesh=zeros(bin_num_max,bin_num_max);
bin_width_contour_input_mesh=zeros(bin_num_max,bin_num_max);
for col_num=1:bin_num_max
    bin_width_contour_terminal_mesh(:,col_num)=bin_width_contour_terminal.';
end
for row_num=1:bin_num_max
    bin_width_contour_input_mesh(row_num,:)=bin_width_contour_input;
end
bin_area_contour=bin_width_contour_terminal_mesh.*bin_width_contour_input_mesh;

for run_num=1:run_num_max
    filename_run_num=sprintf('-%d',run_num);
    for rp_num=1:rp_num_max
        %Singular value decomposition
        r=r_list(1,rp_num);
        p=p_list(1,rp_num);
        index_select_manual=index_select_manual_list(1:p,rp_num);
        fprintf('run_num=%d,r=%d,p=%d\n',run_num,r,p);
        phi=C*expm(A*T)*B;
        [U,S,V]=svds(phi,r);
        W=V*S;
        filename_basic=sprintf('n%dr%dp%dT%d',n,r,p,T);
        filename_initial=sprintf('_pt%dTt%dz%.1f-%.1fci%d',p_tilde,T_tilde,z_min_act,z_max_act,case_num_initial_max);
        filename_bin=sprintf('_t%.1f-%.1fi%.1f-%.1fb%d',terminal_disp_min_ex,terminal_disp_max_ex,input_disp_min_ex,input_disp_max_ex,bin_num_max);
        
        %Computation of response without control
        filename=strcat('alpha_for_initial_state_',sprintf('n%dpt%dz%.1f-%.1fci%d',n,p_tilde,z_min_act,z_max_act,case_num_initial_max),filename_run_num,'.mat');
        if rp_num==1
            population=find((z>=z_min_act).*(z<=z_max_act));
            alpha_tilde_list=zeros(n,case_num_initial_max);
            for case_num_initial=1:case_num_initial_max
                index_tilde=randsample(population,p_tilde);
                kappa=randn([p_tilde 1])+1i*randn([p_tilde 1]);
                alpha_tilde_list(index_tilde,case_num_initial)=kappa;
            end
            save(filename,'alpha_tilde_list')
        else
            load(filename)
        end
        x0_nocontrol_list=expm(A*T_tilde)*B*alpha_tilde_list;
        yT_nocontrol_list=C*expm(A*T)*x0_nocontrol_list;
        
        %Computation of response with optimized actuators
        for type_obj=1:3
            %Actuator selection
            if type_obj==1
                fprintf('Computing response when actuators are placed by trace\n')
                obj=sum(abs(W).^2,2);
                [obj_select,index_select]=maxk(obj,p);
            elseif type_obj==2
                fprintf('Computing response when actuators are placed by determinant\n')
                if p<=r
                    [obj_select,index_select]=F_act_dg_under1(W,p);
                else
                    [obj_select_under,index_select_under]=F_act_dg_under1(W,r);
                    [obj_select,index_select]=F_act_dg_over1(W,p,obj_select_under,index_select_under);
                end
            else
                fprintf('Computing response when actuators are placed manually\n')
                index_select=index_select_manual;
                obj_select=zeros(p,1);
            end
            z_select=z(index_select,1);
        
            %Computation of probability density distribution
            [~,~,count_contour,~]=F_norm_count2(A,B,C,phi,T,index_select,x0_nocontrol_list,yT_nocontrol_list,case_num_initial_max,sigma_min,bin_edge_hist_terminal,bin_edge_hist_input,bin_edge_contour_terminal,bin_edge_contour_input,terminal_disp_min,terminal_disp_max,input_disp_min,input_disp_max);
            prob_dens_contour=count_contour./(case_num_initial_max*bin_area_contour);
            if type_obj==1
                filename=strcat('prob_dens_',filename_basic,filename_initial,filename_bin,'_tr',filename_run_num,'.mat');
            elseif type_obj==2
                filename=strcat('prob_dens_',filename_basic,filename_initial,filename_bin,'_det',filename_run_num,'.mat');
            else
                filename=strcat('prob_dens_',filename_basic,filename_initial,filename_bin,'_manual',filename_run_num,'.mat');
            end
            save(filename,'prob_dens_contour')
        end
        
        %Computation of response with randomly placed actuators
        fprintf('Computing response when actuators are placed randomly\n')
        index_select_list=zeros(p,case_num_act_max);
        for case_num_act=1:case_num_act_max
            index_select_list(:,case_num_act)=randsample(n,p);
        end
        filename=strcat('act_position_rand_',sprintf('n%dp%d_ca%d',n,p,case_num_act_max),filename_run_num,'.mat');
        save(filename,'index_select_list')

        count_hist_terminal_total=zeros(1,bin_num_max);
        count_hist_input_total=zeros(1,bin_num_max);
        count_contour_total=zeros(bin_num_max,bin_num_max);
        count_outside_total=zeros(3,2);
        for case_num_act=1:case_num_act_max
            fprintf('case_num_act=%d/%d\n',case_num_act,case_num_act_max)
            index_select=index_select_list(:,case_num_act);
            [~,~,count_contour,~]=F_norm_count2(A,B,C,phi,T,index_select,x0_nocontrol_list,yT_nocontrol_list,case_num_initial_max,sigma_min,bin_edge_hist_terminal,bin_edge_hist_input,bin_edge_contour_terminal,bin_edge_contour_input,terminal_disp_min,terminal_disp_max,input_disp_min,input_disp_max);
            count_contour_total=count_contour_total+count_contour;
        end
        case_num_max=case_num_initial_max*case_num_act_max;
        prob_dens_contour=count_contour_total./(case_num_max*bin_area_contour);
        filename=strcat('prob_dens_',filename_basic,filename_initial,filename_bin,'_rand',sprintf('_ca%d',case_num_act_max),filename_run_num,'.mat');
        save(filename,'prob_dens_contour')
    end
    
    %Showing probability density distribution
    for rp_num=1:rp_num_max
        close all
        r=r_list(1,rp_num);
        p=p_list(1,rp_num);
        filename_basic=sprintf('n%dr%dp%dT%d',n,r,p,T);
        filename_initial=sprintf('_pt%dTt%dz%.1f-%.1fci%d',p_tilde,T_tilde,z_min_act,z_max_act,case_num_initial_max);
        filename_bin=sprintf('_t%.1f-%.1fi%.1f-%.1fb%d',terminal_disp_min_ex,terminal_disp_max_ex,input_disp_min_ex,input_disp_max_ex,bin_num_max);
        for type_obj=1:4
            if type_obj==1
                filename=strcat('prob_dens_',filename_basic,filename_initial,filename_bin,'_tr',filename_run_num,'.mat');
                filename_type='_tr';
            elseif type_obj==2
                filename=strcat('prob_dens_',filename_basic,filename_initial,filename_bin,'_det',filename_run_num,'.mat');
                filename_type='_det';
            elseif type_obj==3
                filename=strcat('prob_dens_',filename_basic,filename_initial,filename_bin,'_manual',filename_run_num,'.mat');
                filename_type='_manual';
            else
                filename=strcat('prob_dens_',filename_basic,filename_initial,filename_bin,'_rand',sprintf('_ca%d',case_num_act_max),filename_run_num,'.mat');
                filename_type=strcat('_rand',sprintf('_ca%d',case_num_act_max));
            end
            load(filename)
            figure('Position',[300 0.1*scrsz(4) fig_size]);
            hold on
            histogram2('XBinEdges',bin_edge_contour_terminal,'YBinEdges',bin_edge_contour_input,'BinCounts',prob_dens_contour,'DisplayStyle','tile')
            colormap(gca,jet(500))
            colorbar
            xlim([terminal_disp_min terminal_disp_max])
            ylim([input_disp_min input_disp_max])
            xlabel('$\|y(T)\|_2/\|\hat{y}(T)\|_2$','Interpreter','latex');
            ylabel('$\|\alpha\|_2/\|x_0\|_2$','Interpreter','latex')
            ax=gca;
            ax.FontName=font_name;
            ax.FontSize=font_size;
            ax.XScale='log';
            ax.YScale='log';
            ax.ColorScale='log';
            ax.TickDir='out';
            box on
            grid on
            filename=strcat('terminal_input_',filename_basic,filename_initial,sprintf('_t%.1f-%.1fi%.1f-%.1fb%d',terminal_disp_min_ex,terminal_disp_max_ex,input_disp_min_ex,input_disp_max_ex,bin_num_max),filename_type,filename_run_num,'.pdf');
            exportgraphics(gcf,filename)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%