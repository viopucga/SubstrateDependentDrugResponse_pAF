function isLocalMax = find_local_max(vector)
    % Esta función encuentra los máximos locales de un vector
    % Entrada:
    %   vector: un vector numérico
    % Salida:
    %   isLocalMax: un array lógico que indica los máximos locales
    
    % Inicializar el array lógico con valores falsos
    isLocalMax = false(size(vector));
    
    % Iterar a través del vector para marcar los máximos locales
    for i = 2:length(vector)-1
        % Un máximo local ocurre si el valor actual es mayor que sus vecinos
        if vector(i) > vector(i-1) && vector(i) > vector(i+1)
            isLocalMax(i) = true;
        end
    end
end
