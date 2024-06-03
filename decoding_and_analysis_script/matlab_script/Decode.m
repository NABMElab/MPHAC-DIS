function [Mess_in_num] = Decode(Message, varargin)

ignoreerror = 0;

if ~isempty(varargin)
    ignoreerror = 1;
end

Mess_in_num = [];

for j=1:(length(Message)/5)
       
    Mess_chunk{j} = Message(5*j-4:5*j);
    Mess_chunk{j} = upper(Mess_chunk{j});
        
        if Mess_chunk{j}(1) == 'G' 
            Decode1(j) = 0;
        elseif Mess_chunk{j}(1) == 'C'
            Decode1(j) = 1;
        elseif Mess_chunk{j}(1) == 'T'
            Decode1(j) = 2;
        elseif Mess_chunk{j}(1) == 'A'
            Decode1(j) = 3;
        else
            if ignoreerror
                Decode1(j) = -1;
            else
                error('error in Decode1!');
                disp(j)
            end
        end
        
       if strcmp(Mess_chunk{j}(2:3), 'CA') 
            Decode2(j) = 0;
        elseif strcmp(Mess_chunk{j}(2:3), 'CT') 
            Decode2(j) = 1;
        elseif strcmp(Mess_chunk{j}(2:3), 'GA') 
            Decode2(j) = 2;
        elseif strcmp(Mess_chunk{j}(2:3), 'GT')  
            Decode2(j) = 3;
        elseif strcmp(Mess_chunk{j}(2:3), 'TC') 
            Decode2(j) = 4;
        elseif strcmp(Mess_chunk{j}(2:3), 'TG') 
            Decode2(j) = 5;
        elseif strcmp(Mess_chunk{j}(2:3), 'AC') 
            Decode2(j) = 6;
        elseif strcmp(Mess_chunk{j}(2:3), 'AG') 
            Decode2(j) = 7;
       else
           if ignoreerror
                Decode2(j) = -1;
            else
                error('error in Decode2!');
                disp(j)
            end
        end
        
        if strcmp(Mess_chunk{j}(4:5), 'CA') 
            Decode3(j) = 0;
        elseif strcmp(Mess_chunk{j}(4:5), 'CT') 
            Decode3(j) = 1;
        elseif strcmp(Mess_chunk{j}(4:5), 'GA') 
            Decode3(j) = 2;
        elseif strcmp(Mess_chunk{j}(4:5), 'GT')  
            Decode3(j) = 3;
        elseif strcmp(Mess_chunk{j}(4:5), 'TC') 
            Decode3(j) = 4;
        elseif strcmp(Mess_chunk{j}(4:5), 'TG') 
            Decode3(j) = 5;
        elseif strcmp(Mess_chunk{j}(4:5), 'AC') 
            Decode3(j) = 6;
        elseif strcmp(Mess_chunk{j}(4:5), 'AG') 
            Decode3(j) = 7;
        else
            if ignoreerror
                Decode3(j) = -1;
            else
                error('error in Decode3!');
                disp(j)
            end
        end
        
        if Decode1(j) ~= -1 &&  Decode2(j) ~= -1 && Decode3(j) ~= -1
            Decode(j) = Decode1(j)*8^2 + Decode2(j)*8 + Decode3(j);
        else
            Decode(j) = -1;
        end
        
%         Decode(j) = Decode1(j)*8^2 + Decode2(j)*8 + Decode3(j);
        Mess_in_num = [Mess_in_num Decode(j)];
     
end




        
        
            
    