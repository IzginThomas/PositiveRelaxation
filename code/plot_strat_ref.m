                    % tstart = tstart/3600;
                    %tend = tend/3600 + 15;
                    % tend = tend/3600 + 3;
                    tex_temp = tex;%/3600;
                    %tex_temp = tex/3600 + 12;

                    %yex = (diag([9.906e1; 6.624e8; 5.326e11; 1.697e16; 4.000e16; 1.093e9])*yex')';
                    % yex = yex / scaleode;
                    subplot(3,2,5)
                    plot(tex_temp,yex(:,1),plotopt{:})
                    title('$O^{1D}$')
                    xlim([tstart tend])
                    xlabel('t',txtopt{:})

                    subplot(3,2,2)
                    plot(tex_temp,yex(:,2),plotopt{:})
                    title('O')
                    xlim([tstart tend])
                    xlabel('t',txtopt{:})

                    subplot(3,2,3)
                    plot(tex_temp,yex(:,3),plotopt{:})
                    title('$3\cdot O_3$')
                    xlim([tstart tend])
                    xlabel('t',txtopt{:})

                    subplot(3,2,1)
                    plot(tex_temp,yex(:,4) - 2*1.697e+16,plotopt{:})
                    title('$2\cdot(O_2 - 1.697\textrm{e+}16)$')
                    xlim([tstart tend])
                    xlabel('t',txtopt{:})

                    subplot(3,2,4)
                    plot(tex_temp,yex(:,5),plotopt{:})
                    title('NO')
                    xlim([tstart tend])
                    xlabel('t',txtopt{:})

                    subplot(3,2,6)
                    plot(tex_temp,yex(:,6),plotopt{:})
                    title('$2\cdot NO_2$')
                    xlim([tstart tend])
                    xlabel('t',txtopt{:})