classdef MultilinearInterpolator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        d = 0;
        n_x = 0;

        grid = [];

        orders = [];
        values = [];
        mdim_values = [];
        coeff = [];
        B = [];
        smin = [];
        smax = [];
        interp_type = [];
        mvalues = [];
    end
    
    methods
        function self = MultilinearInterpolator(smin, smax, orders)
            
            d = length(orders);

            if d ~= 4
                error('multilinear interpolation only implemented for dimension 4');
            end

            grid = mlinspace(smin,smax,orders);
            
            self.grid = grid;
            self.smin = smin;
            self.smax = smax;
            self.d = d;
            self.orders = orders;
            self.interp_type = 'lin';
        end

        function set_values(self, values)
            self.values = values;
            self.n_x = size(values,2);
        end
        
        function [interp_values] = eval(self, s)
            interp_values = multilinear_interpolation(self.smin,self.smax,self.orders,self.values,s);
        end
        
        
        
    end
    
end

