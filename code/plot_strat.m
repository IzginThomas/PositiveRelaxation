linespec = 'o-';
%t = t/3600 + 12;
%t = t+12;
%yex = (diag([9.906e1; 6.624e8; 5.326e11; 1.697e16; 4.000e16; 1.093e9])*yex')';
tex_temp = tex;%/3600;
ttemp = t;%/3600;
% tex_temp = tex + 12;
% ttemp = t + 12;
subplot(3,2,5)
plot(tex_temp,yex(:,1),plotopt{:})
title('$O^{1D}$',txtopt{:})
ax = gca;
ax.ColorOrderIndex = 1;
grid on;
hold on
plot(ttemp,(Y(1,:))',linespec,plotopt{:})
hold off
xlim([tstart tend])
xlabel('t',txtopt{:})
set(gca,txtopt{:})

subplot(3,2,2)
plot(tex_temp,yex(:,2),plotopt{:})
title('O',txtopt{:})
ax = gca;
ax.ColorOrderIndex = 1;
grid on;
hold on
plot(ttemp,(Y(2,:))',linespec,plotopt{:})
hold off
xlim([tstart tend])
xlabel('t',txtopt{:})
set(gca,txtopt{:})

subplot(3,2,3)
plot(tex_temp,yex(:,3),plotopt{:})
title('$3\cdot O_3$',txtopt{:})
ax = gca;
ax.ColorOrderIndex = 1;
grid on;
hold on
plot(ttemp,(Y(3,:))',linespec,plotopt{:})
hold off
xlim([tstart tend])
xlabel('t',txtopt{:})
set(gca,txtopt{:})

subplot(3,2,1)
plot(tex_temp,yex(:,4) - 2*1.697e+16,plotopt{:})
title('$2\cdot(O_2 - 1.697\textrm{e+}16)$',txtopt{:})
ax = gca;
ax.ColorOrderIndex = 1;
grid on;
hold on
plot(ttemp,(Y(4,:))' - 2*1.697e+16,linespec,plotopt{:})
hold off
xlim([tstart tend])
xlabel('t',txtopt{:})
set(gca,txtopt{:})

subplot(3,2,4)
plot(tex_temp,yex(:,5),plotopt{:})
title('NO',txtopt{:})
ax = gca;
ax.ColorOrderIndex = 1;
grid on;
hold on
plot(ttemp,(Y(5,:))',linespec,plotopt{:})
hold off
xlim([tstart tend])
xlabel('t',txtopt{:})
set(gca,txtopt{:})

subplot(3,2,6)
plot(tex_temp,yex(:,6),plotopt{:})
title('$2\cdot NO_2$',txtopt{:})
ax = gca;
ax.ColorOrderIndex = 1;
grid on;
hold on
plot(ttemp,(Y(6,:))',linespec,plotopt{:})
hold off
xlim([tstart tend])
xlabel('t',txtopt{:})
sgtitle(str,txtopt{:})
set(gca,txtopt{:})
