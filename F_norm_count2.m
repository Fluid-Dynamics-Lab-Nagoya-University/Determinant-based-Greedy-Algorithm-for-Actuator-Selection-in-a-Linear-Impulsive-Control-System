%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%F_norm_count2.m ver.241024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [count_hist_terminal,count_hist_input,count_contour,count_outside]=F_norm_count2(A,B,C,phi,T,index_select,x0_nocontrol_list,yT_nocontrol_list,case_num_initial_max,sigma_min,bin_edge_hist_terminal,bin_edge_hist_input,bin_edge_contour_terminal,bin_edge_contour_input,terminal_disp_min,terminal_disp_max,input_disp_min,input_disp_max)
    %Computation of the terminal output   
    index_select_sort=sort(index_select);
    phis=phi(:,index_select_sort.');
    [Ut,St,Vt]=svd(phis);
    [n,p]=size(phis);
    Stp=zeros(p,n);
    for l=1:p
        if St(l,l)>sigma_min
            Stp(l,l)=1/St(l,l);
        end
    end
    alphas_list=-Vt*Stp*Ut'*yT_nocontrol_list;
    
    Bs=B(:,index_select_sort);
    Bs_alphas_list=Bs*alphas_list;
    x0_control_list=x0_nocontrol_list+Bs_alphas_list;
    
    yT_control_list=C*expm(A*T)*x0_control_list;
    
    %Computation of the 2 norm of the terminal output 
    norm_x0_nocontrol=zeros(1,case_num_initial_max);
    norm_alphas=zeros(1,case_num_initial_max);
    norm_yT_nocontrol=zeros(1,case_num_initial_max);
    norm_yT_control=zeros(1,case_num_initial_max);
    for case_num_initial=1:case_num_initial_max
        norm_x0_nocontrol(1,case_num_initial)=norm(x0_nocontrol_list(:,case_num_initial));
        norm_alphas(1,case_num_initial)=norm(alphas_list(:,case_num_initial));
        norm_yT_nocontrol(1,case_num_initial)=norm(yT_nocontrol_list(:,case_num_initial));
        norm_yT_control(1,case_num_initial)=norm(yT_control_list(:,case_num_initial));
    end
    rate_terminal=norm_yT_control./norm_yT_nocontrol;
    rate_input=norm_alphas./norm_x0_nocontrol;
    
    count_hist_terminal=histcounts(rate_terminal,bin_edge_hist_terminal);%not used
    count_hist_input=histcounts(rate_input,bin_edge_hist_input);%not used
    count_contour=histcounts2(rate_terminal,rate_input,bin_edge_contour_terminal,bin_edge_contour_input);
    
    %Counting the number of cases outside the region (not used)
    count_outside=zeros(3,2);
    count_outside(1,1)=sum(rate_terminal<terminal_disp_min);
    count_outside(1,2)=sum(rate_terminal>terminal_disp_max);
    count_outside(2,1)=sum(rate_input<input_disp_min);
    count_outside(2,2)=sum(rate_input>input_disp_max);
    count_outside(3,1)=case_num_initial_max-sum((terminal_disp_min<=rate_terminal).*(rate_terminal<=terminal_disp_max).*(input_disp_min<=rate_input).*(rate_input<=input_disp_max));%コンターで領域外にプロットされるケース数（領域の外枠線上は領域内と判定）
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%