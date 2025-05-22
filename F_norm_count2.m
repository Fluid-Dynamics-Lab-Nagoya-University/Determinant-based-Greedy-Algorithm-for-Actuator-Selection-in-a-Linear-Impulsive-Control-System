%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%F_norm_count2.m ver.241024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [count_hist_terminal,count_hist_input,count_contour,count_outside]=F_norm_count2(A,B,C,phi,T,index_select,x0_nocontrol_list,yT_nocontrol_list,case_num_initial_max,sigma_min,bin_edge_hist_terminal,bin_edge_hist_input,bin_edge_contour_terminal,bin_edge_contour_input,terminal_disp_min,terminal_disp_max,input_disp_min,input_disp_max)
    %制御入力（alphas_list）の計算    
    %phisの特異値分解
    index_select_sort=sort(index_select);%列ベクトル
    phis=phi(:,index_select_sort.');%phi=C*expm(A*T)*B
    [Ut,St,Vt]=svd(phis);
    %Sigma_tilde_plus(Stp)の計算
    [n,p]=size(phis);
    Stp=zeros(p,n);
    for l=1:p
        if St(l,l)>sigma_min
            Stp(l,l)=1/St(l,l);
        end
    end
    %alphas_listの計算
    alphas_list=-Vt*Stp*Ut'*yT_nocontrol_list;%alphas（列ベクトル）がcase_num_initial_max個並んだ行列．

    %x0_control（制御入力を考慮した初期状態）の計算
    Bs=B(:,index_select_sort);
    Bs_alphas_list=Bs*alphas_list;%Bs*alphas（初期状態のうち，制御のためのインパルス入力によって生み出された部分．列ベクトル）がcase_num_max_initial個並んだ行列．
    x0_control_list=x0_nocontrol_list+Bs_alphas_list;%x0_control（制御入力を考慮した初期状態．列ベクトル）がcase_num_max_initial個並んだ行列．

    %yT_control（制御がある時のyの終端状態）の計算
    yT_control_list=C*expm(A*T)*x0_control_list;%yT_control（制御がある時のyの終端状態．列ベクトル）がcase_num_max_initial個並んだ行列．
    
    %2ノルム（割り算済み）の計算
    %各種ノルムを格納するための変数の準備（normという関数が存在するので変数名に注意）
    norm_x0_nocontrol=zeros(1,case_num_initial_max);%x0_nocontrol（制御入力を考慮していない初期状態）の2ノルム
    norm_alphas=zeros(1,case_num_initial_max);%alphas（制御入力の係数）の2ノルム
    norm_yT_nocontrol=zeros(1,case_num_initial_max);%yT_nocontrol（制御がない時のyの終端状態）の2ノルム
    norm_yT_control=zeros(1,case_num_initial_max);%yT_control（制御がある時のyの終端状態）の2ノルム
    %2ノルムの計算
    for case_num_initial=1:case_num_initial_max
        norm_x0_nocontrol(1,case_num_initial)=norm(x0_nocontrol_list(:,case_num_initial));
        norm_alphas(1,case_num_initial)=norm(alphas_list(:,case_num_initial));
        norm_yT_nocontrol(1,case_num_initial)=norm(yT_nocontrol_list(:,case_num_initial));
        norm_yT_control(1,case_num_initial)=norm(yT_control_list(:,case_num_initial));
    end
    %基準となる2ノルムで割り算
    rate_terminal=norm_yT_control./norm_yT_nocontrol;%yT_control（制御がある時のyの終端状態）の2ノルムを，yT_nocontrol（制御がない時のyの終端状態）の2ノルムで割り算したもの．1×case_num_initialの行ベクトル．
    rate_input=norm_alphas./norm_x0_nocontrol;%alphas（制御入力の係数）の2ノルムを，x0_nocontrol（制御入力を考慮しない初期状態）の2ノルムで割り算したもの．1×case_num_initialの行ベクトル．

    %2ノルム（割り算済み）の分布計算
    count_hist_terminal=histcounts(rate_terminal,bin_edge_hist_terminal);%終端状態の2ノルム（割り算済み）のヒストグラム領域内部のカウント分布．1×bin_numの行ベクトル．
    count_hist_input=histcounts(rate_input,bin_edge_hist_input);%制御入力の2ノルム（割り算済み）のヒストグラム領域内部のカウント分布．1×bin_numの行ベクトル．
    count_contour=histcounts2(rate_terminal,rate_input,bin_edge_contour_terminal,bin_edge_contour_input);%終端状態と制御入力の2ノルム（割り算済み）のコンター領域内部のカウント分布．bin_num×bin_numの行列．

    %領域外のケース数を計算
    count_outside=zeros(3,2);%(3,2)成分は使用しない
    count_outside(1,1)=sum(rate_terminal<terminal_disp_min);%終端状態に関するヒストグラムで，領域の左側にプロットされる（最小目盛より小さくなる）ケース数
    count_outside(1,2)=sum(rate_terminal>terminal_disp_max);%終端状態に関するヒストグラムで，領域の右側にプロットされる（最大目盛より大きくなる）ケース数
    count_outside(2,1)=sum(rate_input<input_disp_min);%制御入力に関するヒストグラムで，領域の左側にプロットされる（最小目盛より小さくなる）ケース数
    count_outside(2,2)=sum(rate_input>input_disp_max);%制御入力に関するヒストグラムで，領域の右側にプロットされる（最大目盛より大きくなる）ケース数
    count_outside(3,1)=case_num_initial_max-sum((terminal_disp_min<=rate_terminal).*(rate_terminal<=terminal_disp_max).*(input_disp_min<=rate_input).*(rate_input<=input_disp_max));%コンターで領域外にプロットされるケース数（領域の外枠線上は領域内と判定）
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%