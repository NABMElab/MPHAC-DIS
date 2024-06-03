function [message] = Remap(input)

number = input-0;
number= double(number);

% if sum(number > 255) > 0
%    error('Input converts to number greater than 255');
% end

message = '';

for i=1:length(number)
   
    temp = '';
    
    if fix(number(i)/64) == 0 
        temp = [temp, 'G'];
    elseif fix(number(i)/64) == 1
        temp = [temp, 'C'];
    elseif fix(number(i)/64) == 2
        temp = [temp, 'T'];
    elseif fix(number(i)/64) == 3
        temp = [temp, 'A'];
    end
        
    
    if mod(floor(number(i)/8),8) == 0
        temp = [temp, 'CA'];
    elseif mod(floor(number(i)/8),8) == 1
        temp = [temp, 'CT'];
    elseif mod(floor(number(i)/8),8) == 2
        temp = [temp, 'GA'];
    elseif mod(floor(number(i)/8),8) == 3
        temp = [temp, 'GT'];
    elseif mod(floor(number(i)/8),8) == 4
        temp = [temp, 'TC'];
    elseif mod(floor(number(i)/8),8) == 5
        temp = [temp, 'TG'];
    elseif mod(floor(number(i)/8),8) == 6
        temp = [temp, 'AC'];
    elseif mod(floor(number(i)/8),8) == 7
        temp = [temp, 'AG'];
    end
    
    
     if rem(number(i),8) == 0
        temp = [temp, 'CA'];
    elseif rem(number(i),8) == 1
        temp = [temp, 'CT'];
    elseif rem(number(i),8) == 2
        temp = [temp, 'GA'];
    elseif rem(number(i),8) == 3
        temp = [temp, 'GT'];
    elseif rem(number(i),8) == 4
        temp = [temp, 'TC'];
    elseif rem(number(i),8) == 5
        temp = [temp, 'TG'];
    elseif rem(number(i),8) == 6
        temp = [temp, 'AC'];
    elseif rem(number(i),8) == 7
        temp = [temp, 'AG'];
     end
    
    message = [message, temp];
end

