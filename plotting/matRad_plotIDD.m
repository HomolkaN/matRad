function f = matRad_plotIDD(doseCube,doseCube2)
if ~exist('doseCube2','var')
    doseCube2 = [];
end
cubeSize = size(doseCube);
IDD = sum(reshape(doseCube,[cubeSize(1),cubeSize(2)*cubeSize(3)]),2);
f = figure;
plot([1:cubeSize(1)],IDD);

if ~isempty(doseCube2)
    cubeSize2 = size(doseCube2);
    if cubeSize2 ~= cubeSize
        error('Dose cube sizes not matching')
    end
    IDD2 = sum(reshape(doseCube2,[cubeSize(1),cubeSize(2)*cubeSize(3)]),2);
    hold on
    plot([1:cubeSize(1)],IDD2);
end
xlabel('Slice number in beam direction');
ylabel('IDD [Gy]');
end