function [X_scaled,X_normalized]=lhsdesign_modified(n,min_ranges_p,max_ranges_p,maxiter)
%lhsdesign_modified is a modification of the Matlab Statistics function lhsdesign.
%It might be a good idea to jump straight to the example to see what does
%this function do.
%The following is the description of lhsdesign from Mathworks documentation
% X = lhsdesign(n,p) returns an n-by-p matrix, X, containing a latin hypercube sample of n values on each of p variables.
%For each column of X, the n values are randomly distributed with one from each interval (0,1/n), (1/n,2/n), ..., (1-1/n,1), and they are randomly permuted.

%lhsdesign_modified provides a latin hypercube sample of n values of
%each of p variables but unlike lhsdesign, the variables can range between
%any minimum and maximum number specified by the user, where as lhsdesign
%only provide data between 0 and 1 which might not be very helpful in many
%practical problems where the range is not bound to 0 and 1
%
%Inputs: 
%       n: number of radomly generated data points
%       min_ranges_p: [1xp] or [px1] vector that contains p values that correspond to the minimum value of each variable
%       max_ranges_p: [1xp] or [px1] vector that contains p values that correspond to the maximum value of each variable
%       maxiter: maximum number of iterations to perform in an attempt to improve the design (default=5)
%Outputs
%       X_scaled: [nxp] matrix of randomly generated variables within the
%       min/max range that the user specified
%       X_normalized: [nxp] matrix of randomly generated variables within the
%       0/1 range 
%
%Example Usage: 
%       [X_scaled,X_normalized]=lhsdesign_modified(100,[-50 100 ],[20  300],10);
%       figure
%       subplot(2,1,1),plot(X_scaled(:,1),X_scaled(:,2),'*')
%       title('Random Variables')
%       xlabel('X1')
%       ylabel('X2')
%       grid on
%       subplot(2,1,2),plot(X_normalized(:,1),X_normalized(:,2),'r*')
%       title('Normalized Random Variables')
%       xlabel('Normalized X1')
%       ylabel('Normalized X2')
%       grid on


p=length(min_ranges_p);
[M,N]=size(min_ranges_p);
if M<N
    min_ranges_p=min_ranges_p';
end
    
[M,N]=size(max_ranges_p);
if M<N
    max_ranges_p=max_ranges_p';
end

slope=max_ranges_p-min_ranges_p;
offset=min_ranges_p;

SLOPE=ones(n,p);
OFFSET=ones(n,p);

for i=1:p
    SLOPE(:,i)=ones(n,1).*slope(i);
    OFFSET(:,i)=ones(n,1).*offset(i);
end
X_normalized = lhs_maximin(n,p,maxiter);

X_scaled=SLOPE.*X_normalized+OFFSET;

% ---------------------
function X = lhs_maximin(n,p,maxiter)
%LHS_MAXIMIN Generate a latin hypercube sample using the max-mi distance criterion.
%   X=LHSDESIGN(N,P,MAXITER) generates a latin hypercube sample X containing N
%   values on each of P variables with MAXITER maximum number of iterations to perform in an
%   attempt to improve the design (default=5). For each column, the N values are
%   randomly distributed with one from each interval (0,1/N), (1/N,2/N),
%   ..., (1-1/N,1), and they are randomly permuted.
%
%   Latin hypercube designs are useful when you need a sample that is
%   random but that is guaranteed to be relatively uniformly distributed
%   over each dimension.

if isnan(maxiter) || isempty(maxiter) || maxiter<0; maxiter = 5; end
% Start with a plain lhs sample over a grid
X = getsample(n,p);

% Create designs, save best one
bestscore = score(X);
for j=2:maxiter
    x = getsample(n,p);
    newscore = score(x);
    if newscore > bestscore
        X = x;
        bestscore = newscore;
    end
end

% ---------------------
function x = getsample(n,p)
x = rand(n,p);
for i=1:p
    x(:,i) = rank(x(:,i));
end
x = x - rand(size(x));
x = x / n;

% ---------------------
function s = score(x)
% compute score function, larger is better

if size(x,1)<2
    s = 0;       % score is meaningless with just one point
    return
end
% Maximize the minimum point-to-point difference
[~,dist] = knnsearch_1(x,x,2);
s = min(dist(:,1));

% -----------------------
function r=rank(x)

% Similar to tiedrank, but no adjustment for ties here
[~, rowidx] = sort(x);
r(rowidx) = 1:length(x);
r = r(:);

% -----------------------
function [idx,D]=knnsearch_1(varargin)
% KNNSEARCH   Linear k-nearest neighbor (KNN) search
% IDX = knnsearch(Q,R,K) searches the reference data set R (n x d array
% representing n points in a d-dimensional space) to find the k-nearest
% neighbors of each query point represented by each row of Q (m x d array).
% The results are stored in the (m x K) index array, IDX. 
%
% IDX = knnsearch(Q,R) takes the default value K=1.
%
% IDX = knnsearch(Q) or IDX = knnsearch(Q,[],K) does the search for R = Q.
%
% Rationality
% Linear KNN search is the simplest appraoch of KNN. The search is based on
% calculation of all distances. Therefore, it is normally believed only
% suitable for small data sets. However, other advanced approaches, such as
% kd-tree and delaunary become inefficient when d is large comparing to the
% number of data points. On the other hand, the linear search in MATLAB is
% relatively insensitive to d due to the vectorization. In  this code, the 
% efficiency of linear search is further improved by using the JIT
% acceleration of MATLAB. Numerical example shows that its performance is
% comparable with kd-tree algorithm in mex.
%
% By Yi Cao at Cranfield University on 25 March 2008

% Check inputs
[Q,R,K,fident] = parseinputs(varargin{:});

[N,M] = size(Q);
L=size(R,1);
idx = zeros(N,K);
D = idx;

if K==1
    % Loop for each query point
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [D(k),idx(k)]=min(d);
    end
else
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [s,t]=sort(d);
        idx(k,:)=t(1:K);
        D(k,:)=s(1:K);
    end
end
if nargout>1
    D=sqrt(D);
end

function [Q,R,K,fident] = parseinputs(varargin)
% Check input and output
error(nargchk(1,3,nargin));

Q=varargin{1};

if nargin<2
    R=Q;
    fident = true;
else
    fident = false;
    R=varargin{2};
end

if isempty(R)
    fident = true;
    R=Q;
end

if ~fident
    fident = isequal(Q,R);
end

if nargin<3
    K=1;
else
    K=varargin{3};
end
