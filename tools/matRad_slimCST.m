function [cst] = matRad_slimCST(cst)

j = 1;
while j <= length(cst)
    if contains(cst{j,2},'h_') || contains(cst{j,2},'h.') || contains(cst{j,2},'H_') || contains(cst{j,2},'H.')
        cst(j,:)=[];
        j = j-1;
    end
    if contains(cst{j,2},'sophagus')
        cst{j,2} = 'Esophagus';
    end
    
    j = j+1;  
end
end