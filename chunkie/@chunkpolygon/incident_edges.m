
function [ichleft,ichright] = incident_edges(chnkr)
%INCIDENT_EDGES for a chunkpolygon, find the incident edges for 
% each vertex. Left means those chunks which precede the vertex
% 
%
% Syntax:
%
% [ichleft,ichright] = incident_edges(chnkr);
%
% Input:
%
% chnkr - a chunkpolygon object
%
% Output:
%
% ichleft - cell array, ichleft{i} is an array of the chunks which 
%               are to the left of the given vertex in the chunk ordering
% ichright - cell array, ichright{i} is an array of the chunks which 
%               are to the right of the given vertex in the chunk ordering
%
%

ichleft={};
ichright={};

[inds,adjs,info] = sortinfo(chnkr);

if (info.ier ~= 0)
    warning('error in sort when finding incident edges, doing nothing');
    return
end

nvert=chnkr.nvert;

ichleft = cell(nvert,1);
ichright = cell(nvert,1);
for i = 1:nvert
    ichleft{i} = [];
    ichright{i} = [];
end

ncomp = info.ncomp;
nchs = info.nchs;

istart=1;
for i = 1:ncomp
    iend = istart+nchs(i)-1;
    ichs = inds(istart:iend); ichs = ichs(:);
    leftvert = -adjs(1,istart);
    rightvert = -adjs(2,iend);
    if (rightvert > 0)
        ichleft{rightvert} = [ichleft{rightvert}; ichs];
    end
    if (leftvert > 0)
        ichright{leftvert} = [ichright{leftvert}; ichs];
    end
    istart = istart + nchs(i);
end

    
end