function  RMSE_temp=RMSEtest(Prediction,measurement)
t1=sum((Prediction(:)-measurement(:)).^2,'omitnan');

if isnan(Prediction)
   Prediction(isnan(Prediction))=[];
end
t3=numel(Prediction);

if isnan(measurement)
   measurement(isnan(measurement))=[];
end
t4=numel(measurement);

t2 = min([t3 t4]);

RMSE_temp = sqrt( t1/ t2 );
end


