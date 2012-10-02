function [z] = multilinear_interpolation(smin,smax,orders,x,y)

	% smin : d x 1 : lower bounds
	% smax : d x 1 : upper bounds
	% orders : d x 1 : number of points in each dimension
	% x : Ns x d : values on the grid
	% y : N x d : points where to interpolate

	% z: interpolated values


	d = length(orders);
	N = size(y,1);
    nx = size(x,2);
	qq = zeros(N,d);
	mm = zeros(N,d);

    
    ss = zeros(size(y));
    for i = 1:d
        ss(:,i) = (y(:,i)-smin(i))/(smax(i)-smin(i));
    end

	for i=1:d
		s = ss(:,i);
		n = orders(i);
		delta = 1/(n-1);
		r = s/delta;
		q = floor(r);
		q = (q<0) .* 0 + (q>= n-2) .* (n-2) + (0<=q).*(q<n-2).*q;
		m = r-q;
		mm(:,i) = m;
		qq(:,i) = q;
	end
	qq = qq+1; % matlab uses indices starting at 1

	[b,g] = strange_construction( x, qq, orders);

%    z0 = zeros(N,nx);
%    for i=1:nx  % TOD: get rid of this loop
%        z0(:,i) = b(:,i) + recursive_evaluation_old(d, {}, mm, g, i);
%    end

    z = zeros(N,nx);
    z = b + recursive_evaluation(d, {}, mm, g);
    
end

function [e] = recursive_evaluation(d,ind,mm,g)
    if length(ind) == d
        S = struct;
        S.type = '()';
        S.subs = [':', ':', ind];
%        e = subsref(g, S);
        e = g(:,:,ind{:});
    else
        j = length(ind);
        ind1 = [ind, 1];
        ind2 = [ind, 2];
        nx = size(g,2);
        tmm = repmat( mm(:,j+1), 1, nx );
        a =  recursive_evaluation(d,ind1,mm,g);
        b = recursive_evaluation(d,ind2,mm,g);
        e = (1-tmm) .* a + tmm .* b;
    end
end

 
function [e] = recursive_evaluation_old(d,ind,mm,g,i)
    if length(ind) == d
        S = struct;
        S.type = '()';
        S.subs = [':', i, ind];
        e = subsref(g, S);
    else
        j = length(ind);
        ind1 = [ind, 1];
        ind2 = [ind, 2];
        e = (1-mm(:,j+1)) .* recursive_evaluation_old(d,ind1,mm,g,i) ...
            + mm(:,j+1) .* recursive_evaluation_old(d,ind2,mm,g,i);
    end
end


function [b,g] = strange_construction(a, q, dims)

	% a : nx x (d1 * ... * dk)
	% q : N x k
	% y : N
	% g : 2**k  ( nx * 2 * ... * 2 )
	
    
	N = size(q,1);
    nx = size(a,2); 

    d = length(dims);
	k = length(dims); % size(q,2);

	q = q-1; % let start indices at 0 ;-)
	
	cdims = cumprod(dims);

	lin_q = q(:,1);
	for i = 2:k
		lin_q = lin_q + q(:,i)*cdims(i-1);
	end;

    cart_prod = cartesianProduct ( repmat( {[0 1]}, 1,d) );

	lin_cp = cart_prod(:,1);
	for i = 2:k
		lin_cp = lin_cp + cart_prod(:,i)*cdims(i-1);
	end;

	lin_q = lin_q + 1;
	
	b = a(lin_q,:);

	g = zeros( N, nx, size(cart_prod,1));
	for i = 1:size(cart_prod,1)
		g(:,:,i) = a(lin_q+lin_cp(i),:) - b;
    end

    ndims = [[ N, nx ], 2*ones(1,d)];
	%g = reshape(g,N,nx,2,2,2,2);
    g = reshape(g,ndims);

	return

end


function result = cartesianProduct(sets)
    % taken from stackoverflow 
    c = cell(1, numel(sets));
    [c{:}] = ndgrid( sets{:} );
    result = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
end
