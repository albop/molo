classdef CompeconInterpolator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        d = 0;
        n_x = 0;
        cdef = {};
        nodes = {};
        orders = [];
        grid = [];
        values = [];
        mdim_values = [];
        coeff = [];
        B = [];
        smin = [];
        smax = [];
        interp_type = [];
    end
    
    methods
        function LI = CompeconInterpolator(smin, smax, orders, interp_type)
            d = length(orders);

            if nargin == 3
                interp_type = 'lin';
            end

            cdef=fundefn(interp_type, orders, smin, smax);
            nodes = funnode(cdef);
            grid = gridmake(nodes);
            
            LI.cdef = cdef;
            LI.grid = grid;
            LI.smin = smin;
            LI.smax = smax;
            LI.d = d;
            LI.nodes = nodes;
            LI.orders = orders;
            LI.interp_type = interp_type;
        end

        function set_values(LI, values)
            [LI.coeff, LI.B]=funfitxy(LI.cdef, LI.grid, values);
            LI.n_x = size(values,2);
            LI.values = values;
            %LI.mdim_values = reshape(values, [ LI.n_x, LI.orders]);
        end
        
        function [interp_values] = eval(LI, s)
            interp_values = funeval(LI.coeff, LI.cdef, s);
        end
        
        
        
    end
    
end

