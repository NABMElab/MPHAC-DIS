function [result] = calc_eff(arr)
result{1} = ToeholdDG(arr{1},'-T',63,'-S',0.22,'-DD');
result{2} = ToeholdDG(arr{6},'-T',63,'-S',0.22,'-DD');
result{3} = onecycle_extension(result{1});
result{4} = onecycle_extension(result{2});
result{5} = equaleff(result{3},result{4},20);


