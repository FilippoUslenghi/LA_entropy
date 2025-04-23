%% Progetto punto 3)
function [entropy] = entropy_calculation(signalVertices_2)
    %Calcolare l'entropia di Shannon dei valori calcolati al punto 2).
    %Per prima cosa stimerei la distribuzione di probabilità semplicemente con un istogramma a prescindere della densità locale
    %dei punti. 
    %In seguito dovremo tenere conto della densità locale dei punti (ci possono essere parti dell'atrio in cui sono stati misurati molti segnali per 
    %cm^2 mentre in altre meno)

    %Distribuzione di probabilità dei valori (values), ovvero dei range
    %calcolati al punto 2
    %hold on
    
    %signalVertices_ = cell2mat(signalVertices_);

    %save('array_signalVertices',signalVertices_, 'v7.3');
    
    % std_signal = std(signalVertices_2, 'omitnan');
    figure()
    % signalVertices_2 = signalVertices_2 * 0.003;
    h = histogram(signalVertices_2,'Normalization','probability');
    %h.BinWidth = 100;
    h.BinWidth = 0.3;
    
    % numBins = h.NumBins;

    %knnh = histogram(knnsignalVertices_,'Normalization','probability');

    % Calculate the entropy

    p = h.Values;
    w = h.BinWidth;

    %if we want to get rid of the zeros, otherwise if we want to
    %consider the zeros too we need to set the 0 to 0.0...1 or the log -> Nan

    p = p(p ~= 0);  
    
    % g = - sum(p.*log(p));

    entropy = log(w) - sum(p.*log(p));
    disp(entropy)
    %entropy2 = - sum(p.*log(p));
    %disp(entropy2)



    %% Bin width

    % optBINS finds the optimal number of bins for a one-dimensional
    % data set using the posterior probability for the number of bins
    % This algorithm uses a brute-force search trying every possible
    % bin number in the given range. This can of course be improved.
    % Generalization to multidimensional data sets is straightforward.
    %
    % Usage:
    % optM = optBINS(data,maxM);
    % Where:
    % data is a (1,N) vector of data points
    % maxM is the maximum number of bins to consider
    %
    % Ref: K.H. Knuth. 2012. Optimal data-based binning for histograms
    % and histogram-based probability density models, Entropy.

    %uno= optBINS(signalVertices_, 114)
end
