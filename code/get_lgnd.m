lgnd = cell(1,2*(length(y0) + size(E,1)));
                for i = 1:length(y0)
                    if exist('visfac','var') && visfac(i) ~= 1
                        lgnd{i} = sprintf('%i $y_{\\textrm{ref},%i}$', visfac(i), i);
                        % lgnd2{i} = sprintf('%i $err_%i$', visfac(i), i);
                    else
                        lgnd{i} = sprintf('$y_{\\textrm{ref},%i}$', i);
                        % lgnd2{i} = sprintf('$err_%i$', i);
                    end
                    if exist('visfac','var') && visfac(i) ~= 1
                        lgnd{length(y0) + size(E,1) + i} = sprintf('%i $y_{\\textrm{num},%i}$', visfac(i), i);
                    else
                        lgnd{length(y0) + size(E,1) + i} = sprintf('$y_{\\textrm{num},%i}$', i);
                    end
                end
                if exist('visfac','var') && visfac(end) ~= 1
                    for i = 1:size(E,1)
                        lgnd{length(y0) + i} = [num2str(visfac(end)) '$\Sigma y_{\textrm{ref},i}$'];
                        lgnd{2*length(y0) + size(E,1) + i} = [num2str(visfac(end)) '$\Sigma y_{\textrm{num},i}$'];
                    end
                else
                    if all(all(E))
                        for i = 1:size(E,1)
                            lgnd{length(y0) + i} = '$\Sigma y_{\textrm{ref},i}$';
                            lgnd{2*length(y0) + size(E,1) + i} = '$\Sigma y_{\textrm{num},i}$';
                        end
                    else
                        for i = 1:size(E,1)
                            lgnd{length(y0) + i} = sprintf('\Sigma e_{%ii}y_i',i);
                            lgnd{2*length(y0) + size(E,1) + i} = sprintf('\Sigma e_{%ii}y_i num',i);
                        end
                    end
                end

                if strcmp(test,'linmod4x400')
                    lgnd{length(y0)+size(E,1)} = '$\mathbf{n}^T \mathbf{y}_{\textrm{ref}}$';
                    lgnd{2*(length(y0)+size(E,1))} = '$\mathbf{n}^T \mathbf{y}_{\textrm{num}}$';
                end

                if ~exist('dtfunc','var')
                    dtfunc = @(dt) dt;
                end
                if ~exist('loc','var')
                    loc = 'ne';
                end
lgnd2 = {'err'};