function [selected_primer, selected_position, selected_dg] = SelectedPrimer(primer_feature, num)
    % 将第三列转换为数值数组
    primer_dg = cell2mat(primer_feature(:,3));

    % 这将返回一个逻辑数组，其中每个元素的值表示对应的行是否满足条件
    rowsToKeep = primer_dg >= -12.5 & primer_dg <= -10.5;

    % 如果没有任何行满足条件，则返回空数组
    if ~any(rowsToKeep)
        selected_primer = {'AAAAAAAAAA'};
        selected_position = {0};
        selected_dg = {0};
        return;  % 提前结束函数执行
    end

    % 找出满足条件的行索引
    rowIndices = find(rowsToKeep);

    % 确保请求的数量不超过满足条件的行数
    num = min(num, length(rowIndices));

    % 从满足条件的行索引中随机选择
    chosen_indices = randsample(rowIndices, num);

    % 根据随机选择的索引提取对应的primer和position
    selected_primer = primer_feature(chosen_indices, 1);
    selected_position = primer_feature(chosen_indices, 2);
    selected_dg = primer_feature(chosen_indices, 3);
end
