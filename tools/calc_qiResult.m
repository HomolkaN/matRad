function qiResult = calc_qiResult(voi,varargin)

for v = 1:length(voi)
    voiName = regexprep(strip(varargin{1}(voi(v)).name), ' ', '_');

    qiResult.(voiName){1,2} = 'mean';
    qiResult.(voiName){1,3} = 'max';
    qiResult.(voiName){1,4} = 'min';
    qiResult.(voiName){1,5} = 'D_2';
    qiResult.(voiName){1,6} = 'D_5';
    qiResult.(voiName){1,7} = 'D_95';
    qiResult.(voiName){1,8} = 'D_98';
    qiResult.(voiName){1,9} = 'CI_2Gy';
    qiResult.(voiName){1,10} = 'HI_2Gy';


    qiResult.(voiName){1,1} = voiName;


    for i = 1:length(varargin)
        qiResult.(voiName){i+1,1} = inputname(i+1);
        qiResult.(voiName){i+1,2} = varargin{i}(voi(v)).mean;
        qiResult.(voiName){i+1,3} = varargin{i}(voi(v)).max;
        qiResult.(voiName){i+1,4} = varargin{i}(voi(v)).min;
        qiResult.(voiName){i+1,5} = varargin{i}(voi(v)).D_2;
        qiResult.(voiName){i+1,6} = varargin{i}(voi(v)).D_5;
        qiResult.(voiName){i+1,7} = varargin{i}(voi(v)).D_95;
        qiResult.(voiName){i+1,8} = varargin{i}(voi(v)).D_98;
        qiResult.(voiName){i+1,9} = varargin{i}(voi(v)).CI_2Gy;
        qiResult.(voiName){i+1,10} = varargin{i}(voi(v)).HI_2Gy;
    end
end



end
