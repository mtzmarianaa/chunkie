classdef chunkpolygon < chunker
%CHUNKPOLYGON subclass of CHUNKER corresponding to a polygonal domain
% In addition to the usual chunker representation of the curve,
% this stores vertex and edge information. 
%

% author: Travis Askham (askhamwhat@gmail.com)


    properties(Access=public)
        igraddatarows
    end
    
    methods
        function obj = chunkpolygon(p)
            if nargin < 1
                p = chunkerpref();
            else
                p = chunkerpref(p);
            end
            obj@chunker(p);
            obj.igraddatarows = [];
        end
                    
        
        function obj = cleardata(obj)
            datakeep = obj.datastor(obj.igraddatarows,:);
            obj.datastor = datakeep;
            obj.igraddatarows = 1:size(datakeep,1);
        end
        
        [ileft,iright] = incident_edges(obj);
        grad = vertgrad(obj);
        grad = vertdsdtgrad(obj);
        
    end
    methods(Static)
        obj = chunkerpoly(verts,cparams,pref,edgevals)
    end
end
