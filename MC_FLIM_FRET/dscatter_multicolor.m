% dscatter_multicolor
% Produce phasor plots for Multi-channel FLIM-FRET analysis. (Plot phasors
% for donor channel, FRET channel and acceptor channel)
% WeiYue Chen, Jan 2014
function [h_dd,h_da,h_a] = dscatter_multicolor(G_dd,S_dd,G_da,S_da,G_AM,S_AM,msize,cmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclude the +-Inf and NAN values so that the accumarray works.
G_dd(find(G_dd==-Inf))=-10000000000;
G_dd(find(G_dd==Inf|isnan(G_dd)==1))=10000000000;
S_dd(find(S_dd==-Inf))=-10000000000;
S_dd(find(S_dd==Inf|isnan(S_dd)==1))=10000000000;

G_da(find(G_da==-Inf))=-10000000000;
G_da(find(G_da==Inf|isnan(G_da)==1))=10000000000;
S_da(find(S_da==-Inf))=-10000000000;
S_da(find(S_da==Inf|isnan(S_da)==1))=10000000000;

G_AM(find(G_AM==-Inf))=-10000000000;
G_AM(find(G_AM==Inf|isnan(G_AM)==1))=10000000000;
S_AM(find(S_AM==-Inf))=-10000000000;
S_AM(find(S_AM==Inf|isnan(S_AM)==1))=10000000000;

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

%%%%%%%%%%%%%%%%%%%%%%
% Define region for color variation for donor channel phasors
minx_da = 0;
maxx_da = 1;
miny_da = 0;
maxy_da = 1;


nbins_da = [min(numel(unique(G_da)),1000) ,min(numel(unique(S_da)),1000) ];

edges1_da = linspace(minx_da, maxx_da, nbins_da(1)+1);
edges1_da = [-Inf edges1_da(2:end-1) Inf];
edges2_da = linspace(miny_da, maxy_da, nbins_da(2)+1);
edges2_da = [-Inf edges2_da(2:end-1) Inf];

[n_da,p_da] = size(G_da);
bin_da = zeros(n_da,2);

[dum,bin_da(:,2)] = histc(G_da,edges1_da);
[dum,bin_da(:,1)] = histc(S_da,edges2_da);
% **bin = bin+1;
H_da = accumarray(bin_da,1,nbins_da([2 1])) ./ n_da; 

H_da = H_da./max(H_da(:));
ind_da = sub2ind(size(H_da),bin_da(:,1),bin_da(:,2));
col_da = H_da(ind_da);
col_use_da = min(64,round(63*(col_da-min(col_da))/(max(col_da)-min(col_da))+1))+64;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define region for color variation for donor channel phasors
minx_a = 0;
maxx_a = 1;
miny_a = 0;
maxy_a = 1;


nbins_a = [min(numel(unique(G_AM)),1000) ,min(numel(unique(S_AM)),1000) ];

edges1_a = linspace(minx_a, maxx_a, nbins_a(1)+1);
edges1_a = [-Inf edges1_a(2:end-1) Inf];
edges2_a = linspace(miny_a, maxy_a, nbins_a(2)+1);
edges2_a = [-Inf edges2_a(2:end-1) Inf];

[n_a,p_a] = size(G_AM);
bin_a = zeros(n_a,2);

[dum,bin_a(:,2)] = histc(G_AM,edges1_a);
[dum,bin_a(:,1)] = histc(S_AM,edges2_a);

H_a = accumarray(bin_a,1,nbins_a([2 1])) ./ n_a; 

H_a = H_a./max(H_a(:));
ind_a = sub2ind(size(H_a),bin_a(:,1),bin_a(:,2));
col_a = H_a(ind_a);
col_use_a = min(64,round(63*(col_a-min(col_a))/(max(col_a)-min(col_a))+1))+128;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_da = scatter(G_da,S_da,msize,'filled','cdata',col_use_da);
colormap(cmap)
hold on
h_dd = scatter(G_dd,S_dd,msize,'filled','cdata',col_use_dd);
colormap(cmap)
h_a = scatter(G_AM,S_AM,2.5,'filled','cdata',col_use_a);



