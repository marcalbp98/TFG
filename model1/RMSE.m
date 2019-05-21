%Computing RMSE
iterations= 1000;
for k= 1:iterations
    kf_mse;
    rmse_a(k)= rmse_raw_a;
    rmse_r(k)= rmse_raw_r;
end
clearvars -except rmse_r rmse_a
Root_Mean_Square_Error_of_rate=mean(rmse_r)
Root_Mean_Square_Error_of_adaptation=mean(rmse_a)