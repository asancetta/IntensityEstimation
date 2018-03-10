function [Lambda, Xb, hT]  = get_integratedIntensity(X, ...
                                             time,...
                                             isJump,...
                                             betaWeib,...
                                             betaG)
                                           
Xb            = X*betaG;                     


durations     = weibullMeasure(time, isJump,betaWeib,0);
Lambda        = Xb.*durations;
                                           
indFirstGood  = find( ~isnan(Lambda), 1, 'first');                                          
Lambda(indFirstGood:end) = cumsum(Lambda(indFirstGood:end));

end


