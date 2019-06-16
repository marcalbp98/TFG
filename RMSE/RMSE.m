%Computing RMSE
iterations= 1000;
cd '/home/marshi/Desktop/Uni/BDBI/3r/TFG/model1'
for i= 1:iterations
    prm_1_kf_3;
    armse_raw_a(i)= rmse_raw_a;
    armse_raw_r(i)= rmse_raw_r;
    armse_kf_a(i)= rmse_kf_a;
    armse_kf_r(i)= rmse_kf_r;
end

Model1_Raw_Root_Mean_Square_Error_of_rate=mean(armse_raw_r)
Model1_Raw_Root_Mean_Square_Error_of_adaptation=mean(armse_raw_a)
Model1_kf_Root_Mean_Square_Error_of_rate=mean(armse_kf_r)
Model1_kf_Root_Mean_Square_Error_of_adaptation=mean(armse_kf_a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/home/marshi/Desktop/Uni/BDBI/3r/TFG/model2'
for i= 1:iterations
    prm_2_kf_3;
    brmse_raw_a(i)= rmse_raw_a;
    brmse_raw_r(i)= rmse_raw_r;
    brmse_raw_s(i)= rmse_raw_s;
    brmse_kf_s(i)= rmse_kf_s;
    brmse_kf_a(i)= rmse_kf_a;
    brmse_kf_r(i)= rmse_kf_r;
end
Model2_Raw_Root_Mean_Square_Error_of_rate=mean(brmse_raw_r)
Model2_Raw_Root_Mean_Square_Error_of_adaptation=mean(brmse_raw_a)
Model2_Raw_Root_Mean_Square_Error_of_depression=mean(brmse_raw_s)
Model2_kf_Root_Mean_Square_Error_of_rate=mean(brmse_kf_r)
Model2_kf_Root_Mean_Square_Error_of_adaptation=mean(brmse_kf_a)
Model2_kf_Root_Mean_Square_Error_of_depression=mean(brmse_kf_s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/home/marshi/Desktop/Uni/BDBI/3r/TFG/model3'
for i= 1:iterations
    prm_3_kf_3;
    crmse_raw_re(i)= rmse_raw_re;
    crmse_raw_ri(i)= rmse_raw_ri;
    crmse_raw_ae(i)= rmse_raw_ae;
    crmse_raw_ai(i)= rmse_raw_ai;
    crmse_raw_se(i)= rmse_raw_se;
    crmse_raw_si(i)= rmse_raw_si;
    crmse_kf_re(i)= rmse_kf_re;
    crmse_kf_ri(i)= rmse_kf_ri;
    crmse_kf_ae(i)= rmse_kf_ae;
    crmse_kf_ai(i)= rmse_kf_ai;
    crmse_kf_se(i)= rmse_kf_se;
    crmse_kf_si(i)= rmse_kf_si;
end
Model3_Raw_Root_Mean_Square_Error_of_rate_e=mean(crmse_raw_re)
Model3_Raw_Root_Mean_Square_Error_of_rate_i=mean(crmse_raw_ri)
Model3_Raw_Root_Mean_Square_Error_of_adaptation_e=mean(crmse_raw_ae)
Model3_Raw_Root_Mean_Square_Error_of_adaptation_i=mean(crmse_raw_ai)
Model3_Raw_Root_Mean_Square_Error_of_depression_e=mean(crmse_raw_se)
Model3_Raw_Root_Mean_Square_Error_of_depression_i=mean(crmse_raw_si)
Model3_kf_Root_Mean_Square_Error_of_rate_e=mean(crmse_kf_re)
Model3_kf_Root_Mean_Square_Error_of_rate_i=mean(crmse_kf_ri)
Model3_kf_Root_Mean_Square_Error_of_adaptation_e=mean(crmse_kf_ae)
Model3_kf_Root_Mean_Square_Error_of_adaptation_i=mean(crmse_kf_ai)
Model3_kf_Root_Mean_Square_Error_of_depression_e=mean(crmse_kf_se)
Model3_kf_Root_Mean_Square_Error_of_depression_i=mean(crmse_kf_si)
