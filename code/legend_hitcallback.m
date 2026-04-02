%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 10/26/2017                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function legend_hitcallback(src,evnt)

if strcmp(evnt.Peer.Visible,'on')
  evnt.Peer.Visible = 'off';
else
  evnt.Peer.Visible = 'on';
end

end