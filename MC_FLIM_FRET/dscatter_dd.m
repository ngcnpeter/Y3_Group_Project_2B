% dscatter_dd
% Produce phasor plots for Multi-channel FLIM-FRET analysis. (Only plot
% phasors for donor channel.)
% WeiYue Chen, Jan 2014
function [h_dd,h_da,h_a] = dscatter_dd(G_dd,S_dd,msize,cmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclude the +-Inf and NAN values so that the accumarray works.
G_dd(find(G_dd==-Inf))=-10000000000;
G_dd(find(G_dd==Inf|isnan(G_dd)==1))=10000000000;
S_dd(find(S_dd==-Inf))=-10000000000;
S_dd(find(S_dd==Inf|isnan(S_dd)==1))=10000000000;
% Define region for color variation for donor channel phasors
minx_dd = 0;
maxx_dd = 1;
miny_dd = 0;
maxy_dd = 1;

nbins_dd = [min(numel(unique(G_dd)),1000) ,min(numel(unique(S_dd)),1000) ];

edges1_dd = linspace(minx_dd, maxx_dd, nbins_dd(1)+1);
edges1_dd = [-Inf edges1_dd(2:end-1) Inf];
edges2_dd = linspace(miny_dd, maxy_dd, nbins_dd(2)+1);
edges2_dd = [-Inf edges2_dd(2:end-1) Inf];

[n_dd,p_dd] = size(G_dd);
bin_dd = zeros(n_dd,2);

[dum,bin_dd(:,2)] = histc(G_dd,edges1_dd);
[dum,bin_dd(:,1)] = histc(S_dd,edges2_dd);

H_dd = accumarray(bin_dd,1,nbins_dd([2 1])) ./ n_dd; 

H_dd = H_dd./max(H_dd(:));
ind_dd = sub2ind(size(H_dd),bin_dd(:,1),bin_dd(:,2));
col_dd = H_dd(ind_dd);
col_use_dd = min(64,round(63*(col_dd-min(col_dd))/(max(col_dd)-min(col_dd))+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_dd = scatter(G_dd,S_dd,msize,'filled','cdata',col_use_dd);
colormap(cmap)




