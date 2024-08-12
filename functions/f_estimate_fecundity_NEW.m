%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [total_fecundity, mean_HT, var_HT] = f_estimate_fecundity_NEW(coral_cm2, HT, F_list, fecund_min_size, a, b)

% Select colony sizes with 100% gravid
I=find(coral_cm2>= fecund_min_size);

if isempty(I)==0
    adult_coral_sizes = coral_cm2(I);
    adult_relfitness = F_list(I);
    
    % Allometric relationship based on Hall and Hughes 1996:
    % Egg volume = exp(a+b*log(size))
    
    % all_egg_volumes = exp(a + b*log(adult_relfitness.*adult_coral_sizes)) ; % mm3 of eggs produced by each colony
    % Correction 29/09/2022: logs in Hall&Hughes were log10, not natural logs!!
    % Also, makes more sense to apply relative fitness on number of eggs rather than colony size (log...)
    all_egg_volumes = adult_relfitness.*(10.^(a + b*log10(adult_coral_sizes))) ; % mm3 of eggs produced by each colony
    
    total_fecundity = floor(sum(sum(all_egg_volumes))/0.1) ; %0.1 mm3 is the average volume of an egg
    HT_fecund = full(HT(I)); % need to convert sparse matrix HT
    mean_HT = sum(all_egg_volumes.*HT_fecund)/sum(all_egg_volumes);
    var_HT = var(HT_fecund,all_egg_volumes);
 
else % if no fecund corals, then no larvae, and no HT
    total_fecundity = 0;
    mean_HT = NaN;
    var_HT = NaN;
end
