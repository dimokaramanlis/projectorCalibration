

function curr_corrected = scalePhotodiodeToNanoAmps(curr, unit)

curr_corrected = zeros(size(curr));

for ii = 1:numel(curr_corrected)
    
    switch unit{ii}
        case {'nA'}
            scale = 1;
        case {'µA','uA'}
            scale = 1e3;
        case {'mA'}
            scale = 1e6;
        case {'sA'}
            scale = 1;
    end
    curr_corrected(ii) = curr(ii) * scale;
end
    
end


