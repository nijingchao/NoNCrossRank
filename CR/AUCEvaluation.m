%% AUC value evaluation
function AUCEvaluation(RankRecord, ExpandSeeds, AllGeneID)

ROCn = zeros(6, length(RankRecord));
topn = zeros(length(RankRecord{1}), length(RankRecord));

for j = 1:length(ExpandSeeds)
    
    for k = 1:length(RankRecord{1})
        
        real_row = AllGeneID(RankRecord{j}(k)); % ID of gene at rank k
        
        if real_row == ExpandSeeds(j)
            topn(k, j) = 1;
        else
            topn(k, j) = 0;
        end
        
    end
    
    ROCn(1,j) = AUCValue(topn(:,j), 50);
    ROCn(2,j) = AUCValue(topn(:,j), 100);
    ROCn(3,j) = AUCValue(topn(:,j), 300);
    ROCn(4,j) = AUCValue(topn(:,j), 500);
    ROCn(5,j) = AUCValue(topn(:,j), 700);
    ROCn(6,j) = AUCValue(topn(:,j), 1000);
    
end

avg_ROCn = mean(ROCn,2);
disp(['AUC50:' num2str(avg_ROCn(1))]);
disp(['AUC100:' num2str(avg_ROCn(2))]);
disp(['AUC300:' num2str(avg_ROCn(3))]);
disp(['AUC500:' num2str(avg_ROCn(4))]);
disp(['AUC700:' num2str(avg_ROCn(5))]);
disp(['AUC1000:' num2str(avg_ROCn(6))]);

end