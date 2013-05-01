   
function [depth, median, inf, sup]=wfbplotImages(data,weight,depth,method)
% Produces functional boxplots or enhanced functional boxplots of the given functional data. 
% It can also be used to carry out functional data ordering based on band depth. 
%
%Arguments
%
%	data: a n-by-p-by-q-by... functional data matrix where n is the number of images, and p-by-q-by... is the demension of the images
%         alternatively, a functional data object or functional parameter
%         object can also be specified. 
%	x:	the x coordinates of curves. Defaults to 1:p where p is the number of x coordinates.
%	method: the method to be used to compute band depth. Can be one of "BD2", "BD3", "MBD" or "Both" 
%	        with a default of "MBD".
%	depth:	a vector giving band depths of curves. If missing, band depth computation is conducted.
%	show: logical. If TRUE (the default) then a functional boxplot is produced. If not, band depth
%	      and outliers are returned.
%	prob: a vector giving the probabilities of central regions in a decreasing order, then an enhanced 
%	      functional boxplot is produced. Defaults to be 0.5 and a functional boxplot is plotted.
%	color:	a vector giving the colors of central regions from light to dark for an enhanced functional 
%	       boxplot. Defaults to be magenta for a functional boxplot.
%	outliercol: color of outlying curves. Defaults to be red.
%	barcol: color of bars in a functional boxplot. Defaults to be blue.
%	fullout: logical for plotting outlying curves. If FALSE (the default) then only the part outside the
%	         box is plotted. If TRUE, complete outling curves are plotted.
%	factor: the constant factor to inflate the middle box and determine fences for outliers. Defaults to 
%	        be 1.5 as in a classical boxplot.
%
%Details
%
%	For functional data, the band depth (BD) or modifed band depth (MBD) allows for ordering a sample of 
%	curves from the center outwards and, thus, introduces a measure to define functional quantiles and 
%	the centrality or outlyingness of an observation. A smaller rank is associated with a more central 
%	position with respect to the sample curves. "BD2" uses two curves to determine a band and "BD3" uses 
%	three curves. BD usually provides many ties (curves have the same depth values), but MBD does not. 
%	The method "Both" uses BD2 first and then uses MBD to break ties.
%
%Value
%
%	depth:	band depths of given curves.
%	outpoint: column indices of detected outliers.
%
%Author(s)
%
%	Ying Sun: sunwards@stat.tamu.edu
%
%	Marc G. Genton: genton@stat.tamu.edu
%
%References
%
%	Sun, Y. and Genton, M. G. (2011), "Functional Boxplots," Journal of Computational and Graphical Statistics,
%	to appear.
%
%	Lopez-Pintado, S. and Romo, J. (2009), "On the concept of depth for functional data," Journal of the American
%	Statistical Association, 104, 718-734.
%

% %%default values
% 
% factor=1.5; 
% fullout='False'; 
% barcol='b'; 
% outliercol='r'; 
% color='m'; 
% prob=0.5; 
% show='True';
% method='MBD'; 
% depth=[]; 

% Examples
%
% 	clear all;
%
% 	ncasem = 39;
% 	ncasef = 54;
% 	nage   = 31;
% 
% 	fid = fopen('hgtm.dat','rt');
% 	hgtmmat = reshape(fscanf(fid,'%f'),[nage,ncasem]);
% 
% 	fid = fopen('hgtf.dat','rt');
% 	hgtfmat = reshape(fscanf(fid,'%f'),[nage,ncasef]);
% 
% 	age = [ 1:0.25:2, 3:8, 8.5:0.5:18 ]';
%
%%fbplot of boys' height
%	fbplot(hgtmmat,age,depth,method,show,prob,color,outliercol,barcol,fullout,factor)
%	xlim([0.5,18.5])
%	ylim([60,200])
%	xlabel('Age (Years)')
%	ylabel('Height (cm)')
%	title('Boys')
%
%%fbplot of girls' height
%	fbplot(hgtfmat,age,depth,method,show,prob,color,outliercol,barcol,fullout,factor)
%	xlim([0.5,18.5])
%	ylim([60,200])
%	xlabel('Age (Years)')
%	ylabel('Height (cm)')
%	title('Girls')
 
function combinat=combinat(n,p)
if n<p 
combinat=0;
else
   combinat=nchoosek(n,p);
end
end

function resultados=estaEntre(v,matrizDatos)
% The input in this function is the data matrix and the pair of indexes and the output is a vector of dimension n, where the coordinate k takes the value 1 or 0 whether the observation k is inside the band or not.
dim = length( size(matrizDatos) ) - 1;
Z = matrizDatos;
resultados = zeros(1, size(matrizDatos, 1));
if dim == 2
	[inf, sup] = band( Z(v, :, :) );
	[n, p, q] = size(matrizDatos);
	for iI = 1:n
		tmp = reshape(Z(iI, :, :), p, q);
		resultados(iI) = ( sum( sum( and(tmp <= sup, tmp >= inf) ) ) == (p*q) );
	end
elseif dim == 3
	[inf, sup] = band( Z(v, :, :, :) );
	[n, p, q, r] = size(matrizDatos);
	for iI = 1:n
		tmp = reshape(Z(iI, :, :, :), p, q, r);
		resultados(iI) = ( sum( sum( sum( and( tmp <= sup, tmp >= inf ) ) ) ) == (p*q*r) );
	end
else
	error('The size of the image is not correct: %d', dim );
end
end

function [contg]=BD2(matrizDatos, weight)
% With this function we calculate the band depth of every observation in the matrix of data: matrizDatos. The rows are the observations and the columns are the variables.
tic
n=size(matrizDatos, 1);
cont=zeros(1,n); % The depths are initialized to 0.
sum_weight = 0;
%for i=1:((n+2)/3)
for i=1:(n-1)
   for j=(i+1):(n) % We choose pairs of indexes in all the possible ways.
      cont=cont + weight(i) * weight(j) * estaEntre([i j],matrizDatos); % With this subfuntion we check whether the observations are inside the band defined by observation i and j.  
      sum_weight = sum_weight + weight(i) * weight(j);               
   end
end	
   contg=cont/combinat(n,2)/sum_weight;
toc
end

function [contg]=BD3(matrizDatos, weight)
% With this function we calculate the band depth with J=3. The imput is the data matrix, where the rows are the observations and the columns the variables. 
tic
n=size(matrizDatos, 1); 
cont=zeros(1,n);      % Initialize the depth of each observation to zero.
                      % Select three observations from the sample in all the possible ways. 
sum_weight = 0;
for i=1:(n-2)
   for j=(i+1):(n-1)
      for k=(j+1):n
         cont = cont + weight(i) * weight(j) * weight(k) * estaEntre([i j k],matrizDatos);   % In this subfunction we check which observations from the sample is inside the band delimeted by observations i,j and k.           
         sum_weight = sum_weight + weight(i) * weight(j) * weight(k);
      end
   end
end	
contg=cont/combinat(n,3)/sum_weight;
toc
end

function resultado=a(v,matrizDatos) 
dim = length( size(matrizDatos) ) - 1;
Z = matrizDatos;
resultado = zeros(1, size(matrizDatos, 1));
if dim == 2
	[inf, sup] = band( Z(v, :, :) );
	[n, p, q] = size(matrizDatos);
	for iI = 1:n
		tmp = reshape(Z(iI, :, :), p, q);
		resultado(iI) = sum( sum( and(tmp <= sup, tmp >= inf) ) ) / (p*q);
	end
elseif dim == 3
	[inf, sup] = band( Z(v, :, :, :) );
	[n, p, q, r] = size(matrizDatos);
	for iI = 1:n
		tmp = reshape(Z(iI, :, :, :), p, q, r);
		resultado(iI) = sum( sum( sum( and( tmp <= sup, tmp >= inf ) ) ) ) / (p*q*r);
	end
else
	error('The size of the image is not correct: %d', dim );
end
end

function [contg]=MBD(matrizDatos, weight)
% This function calculates the generalized band depth of a set of data
tic
n=size(matrizDatos, 1); % size of the data matrix
cont=zeros(1,n); % The initial value of the depth of the observations from the matrix is zero
sum_weight = 0;
for i=1:(n-1)
    for j=(i+1):(n) % consider all possible pairs of functions
        
        cont = cont + weight(i) * weight(j) * a([i j],matrizDatos);  % In cont we save the generalized band depth of each observation (function)
        sum_weight = sum_weight + weight(i) * weight(j);
    end
end
contg=cont;
contg=cont/combinat(n,2)/sum_weight;
toc
end

function [bandInf, bandSup] = band( data )
% This function computes the band of the given data
if length(size(data)) == 3
	bandInf = ones( size(data, 2), size(data, 3) );
	bandSup = zeros( size(data, 2), size(data, 3) );	
elseif length(size(data)) == 4
	bandInf = ones( size(data, 2), size(data, 3), size(data, 4) );
	bandSup = zeros( size(data, 2), size(data, 3), size(data, 4) );
else
	error('The size of the data is not correct: %d\n', size(data) );
end
for iI = 1:size(data, 1)
	if length(size(data)) == 3
		bandInf = bitand( bandInf, reshape(data(iI, :, :), size(data, 2), size(data,3)) );
		bandSup = bitor(  bandSup, reshape(data(iI, :, :), size(data, 2), size(data,3)) );
    elseif length(size(data)) == 4
		bandInf = bitand( bandInf, reshape(data(iI, :, :, :), size(data, 2), size(data, 3), size(data, 4)) );
		bandSup = bitor(  bandSup, reshape(data(iI, :, :, :), size(data, 2), size(data, 3), size(data, 4)) );
	else
		error('The size of the image is not correct: %d\n', size(data)-1 );
	end
end
end

%%default values

if nargin<4, 	method='MBD'; end
if nargin<3, 	depth=[]; end
if nargin<2,	weight=ones(size(data,1), 1)/size(data,1); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute band depth for images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(depth)
   if strcmp(method,'BD2')
      depth=BD2(data, weight);
   elseif strcmp(method,'BD3')
      depth=BD3(data, weight);
   elseif strcmp(method,'MBD')
      depth=MBD(data, weight);
   elseif strcmp(method,'Both')
      depth=round(BD2(data, weight)*10000)+MBD(data, weight);  
   end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort band depth and find the median, the outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[dp_s,index]=sort(depth,'descend');
weight_sum = 0;
index_50_flag = 0;
index_993_flag = 0;
for iI = 1:length(index)
	weight_sum = weight_sum + weight( index( iI ) );
	if weight_sum >= 0.993 && index_993_flag == 0
		index_993 = iI;
		index_993_flag = 1;
	elseif weight_sum >= 0.5 && index_50_flag == 0
		index_50 = iI;
		index_50_flag = 1;
	end	
end
m = index_50;
centralRegionId = index(1:m);
dim = length(size(data))-1;
if dim == 2
    median = reshape(data(index(1), :, :), size(data,2), size(data,3));
	center=data(index(1:m), :, :);
	out=data(index((m+1):end), :, :);
elseif dim == 3
    median = reshape(data(index(1), :, :), size(data,2), size(data,3));
	center=data(index(1:m), :, :, :);
    out=data(index((m+1):end), :, :, :);
end
[inf, sup] = band(center);
depth
index

% find the outerliers

%{    
% compute the 5% and 95% whiskers
    
     if prob(pp)==0.5 %check outliers
	dist=factor*(sup-inf);
	upper=sup+dist;
	lower=inf-dist;
	%upper = sup_95;
 	%lower = inf_95;
	%outlier column
	outly=sum(or(data<lower'*ones(1,n),data>upper'*ones(1,n)));

	% make the curves outside of 0.993 weights also be outlier
	outly(index(index_993:end)) = 1;

	outpoint=find(outly);
	out=data(:,outpoint);
	weight_out = weight(outpoint);
        good=data;
        good(:,outpoint)=[];
	maxcurve=max(good,[],2)';
	mincurve=min(good,[],2)';
	%dataBoxplot(1, :) = mincurve;
	%dataBoxplot(2, :) = inf;
	%dataBoxplot(3, :) = data(:, index(1));
	%dataBoxplot(4, :) = sup;
	%dataBoxplot(5, :) = maxcurve;
	if sum(outly)>0
		if show 
			%plot(x,out, '--r', 'LineWidth', 2);
			hold all;
			weight_out = weight_out ./ sum( weight_out);
			for iOut = 1:length(weight_out)
				plot( x, out(:, iOut), 'Color', [1.0-weight_out(iOut) 1.0-weight_out(iOut) 1.0-weight_out(iOut)], 'LineStyle', '--' );
			end
			%xflip = [x(1 : end - 1)' fliplr(x')];
			%for iOut = 1:length(weight_out)
			%	yflip = [out(1:end-1,iOut)' fliplr(out(:,iOut)')];
			%	patch(xflip, yflip, 'r', 'EdgeAlpha', 0.1, 'FaceColor', 'none');
			%end
		end
	end
	barval=(x(1)+x(tp))/2;
        loc=find(sort([x;barval])==barval);
	bar=loc(1);

	if show
            hold all;
            [xinv,xindex]=sort(x,'descend');
            xx=[x;xinv];
            supinv=sup(xindex);
            yy=[inf,supinv];
            h=fill(xx,yy,color(pp,:));
            if prob(pp)==0.5
                set(h,'edgecolor',barcol,'LineWidth',2);
            else
                set(h,'edgecolor','none');
            end
        end

	if show
            hold all;
            line([x(bar) x(bar)],[maxcurve(bar) sup(bar)],'Color',barcol,'LineWidth',2);
            hold all;
            line([x(bar) x(bar)],[mincurve(bar) inf(bar)],'Color',barcol,'LineWidth',2);
        end
            
%             if show 
%                 hold all;
%                 plot(x,maxcurve,'Color','b','LineWidth',2);
%                 plot(x,mincurve,'Color','b','LineWidth',2);
%                 if fullout
%                     if sum(outly)>0 
%                         hold all;
%                         plot(x,out,'Color',outliercol, 'LineStyle', '--', 'LineWidth', 2);
%                     end
%                 end
%             end
            
        end
        
%         if show
%             hold all;
%             plot(x,data(:,index(1)),'Color','k', 'LineStyle', '-', 'LineWidth',2);
%         end
		
% 		if show 
% 			hold all;
%             [xinv,xindex]=sort(x,'descend');
%             xx=[x;xinv];
%     		supinv=sup(xindex);
%             yy=[inf,supinv];
% 			h=fill(xx,yy,color(pp,:));
% 			if prob(pp)==0.5
%                 set(h,'edgecolor',barcol,'LineWidth',2);
%             else 
%                 set(h,'edgecolor','none');
% 			end
% 		end
%  		if show 
%  			hold all;
%  			plot(x,data(:,index(1)),'Color','k','LineWidth',2);
%  			plot(x,maxcurve,'Color','b','LineWidth',2);
%  			plot(x,mincurve,'Color','b','LineWidth',2);
%  			if fullout
%  				if sum(outly)>0 
%  					hold all;
%  					plot(x,out,'Color',outliercol);
%  				end
%  			end
%  		end

 		if fullout
 			if sum(outly)>0 
 				%hold all;
 				%plot(x,out,'Color',outliercol);
 			end
 		end
    end
    
    if show 
        hold all;
 		plot(x,data(:,index(1)),'Color','k','LineWidth',2);
 		plot(x,maxcurve,'Color','b','LineWidth',2);
 		plot(x,mincurve,'Color','b','LineWidth',2);
    end
    hold off;
    depth;
    outpoint;
    medianCurveId = index(1);

    % find the curve's position in the boxplot
    sum_weight = 0;
    for iI = 1:length(weight)
	curve_down = min( data( :, index(1:iI) ), [], 2 );
	if sum( sampleData < curve_down )
        	sum_weight = sum_weight + weight( index(iI) );
	else
	        break;
	end
    end
    posPercentile = 1 - sum_weight;
    if iI >= length( weight )
    	curve_down = min( data, [], 2 );
	if sum( sampleData < curve_down )
		[value, pos] = min( sampleData - curve_down );
		posPercentile = value / curve_down(pos);	
	end
    end
%}
end
