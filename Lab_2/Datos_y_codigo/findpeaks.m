function [peaks, locations] = findpeaks(data, min_peak_distance)
    % Encuentra picos locales en un vector de datos
    % data: Vector de datos de entrada
    % min_peak_distance: Distancia mínima entre picos
    
    % Inicializa las salidas
    peaks = [];
    locations = [];
    
    % Número de puntos en los datos
    num_points = length(data);
    
    % Verifica que el tamaño mínimo de la distancia sea al menos 1
    if nargin < 2
        min_peak_distance = 1;
    end
    
    % Itera sobre los puntos del vector de datos
    for i = 2:num_points-1
        if data(i-1) < data(i) && data(i) > data(i+1)
            % Si el punto es un pico, agrega a las salidas
            peaks = [peaks; data(i)];
            locations = [locations; i];
        end
    end
    
    % Elimina picos demasiado cercanos
    if min_peak_distance > 1
        [locations, sort_idx] = sort(locations);
        peaks = peaks(sort_idx);
        valid_idx = true(length(locations), 1);
        
        for i = 2:length(locations)
            if locations(i) - locations(i-1) < min_peak_distance
                valid_idx(i) = false;
            end
        end
        
        locations = locations(valid_idx);
        peaks = peaks(valid_idx);
    end
end
