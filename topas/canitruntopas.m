function minHistoriesForAllSpots = canitruntopas(w,runs)

if ~exist("runs","var")
    runs = 1;
end
minRelWeight = 0;
minParticlesBixel = round(max([minRelWeight*1e6*max(w(:)),1]));

minHistoriesForAllSpots = ceil(max(sum(w)./w*runs*abs(minParticlesBixel-0.5))/10000)*10000;

end

