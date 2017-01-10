% min ||B.*LR'-S||+lambda*(||L||+||R'||)+||H*L*R'||+||L*R'*T||
% L:n*r R t*r
function [L,R,v]=estics(S,B,r,lambda,MaxIter)
	if nargin < 4, lambda = 1e-3;     end
	if nargin < 5, MaxIter = 50;      end

	[n,t]=size(S);

	v=1e50;
	Lopt=rand(n,r);
	
	Il=diag(ones(r,1));
	Ir=diag(ones(t,1));
	Zt=zeros(t);
	Zn=zeros(n);
	Ot=ones(t);
	On=ones(n);
	for i=1:MaxIter
		Ropt=myInverse(B,Lopt,lambda,S)';
		Lopt=myInverse(B',Ropt,lambda,S')';
		%Ropt=myInverseFeature(B,Lopt,lambda,S,H)';
		%Lopt=myInverseFeature(B',Ropt,lambda,S',T')';
		X=Lopt*Ropt';
        vopt1=ev(B.*X-S);
        vopt2=lambda*(ev(Lopt)+ev(Ropt'));
        % vopt3=ev(H*X);
        % vopt4=ev(X*T);
        vopt=vopt1+vopt2;
        % vopt=vopt1+vopt2+vopt3+vopt4;
		%vopt=ev(B.*X-S)+lambda*(ev(Lopt)+ev(Ropt'))+ev(H*X)+ev(X*T);
        %disp([vopt1 vopt2 vopt3 vopt4]);
		if vopt<v
			v=vopt;
			L=Lopt;
			R=Ropt;
		end
	end

function y=ev(X)
	Xt=X.*X;
	y=sqrt(sum(Xt(:)));

% best fit [B.*(LY); sqrt(lambda)*Y]=[S;0]
% B(n*t) S(n*t) L(n*r) Y(r*t)
% A(m*k), B(m*n), M(m*n), Y(k*n), minimizing ||B.*(LY-S)||
function Y = myInverse(B,L,lambda,S)
	[n,t]=size(B);
	r=size(L,2);
	Y=zeros(r,t);
	Ir=diag(ones(r,1));
	Zr=zeros(r,1);
	for i=1:t
		Bi=diag(B(:,i));
		P=[Bi*L; sqrt(lambda)*Ir];
		Q=[S(:,i); Zr];
		Y(:,i)=(P'*P)\(P'*Q);
    end
    
% best fit [B.*(LY); sqrt(lambda)*Y; FLY]=[S;0;0]
function Y = myInverseFeature(B,L,lambda,S,F)
	[n,t]=size(B);
	r=size(L,2);
    fx=size(F,1);
	Y=zeros(r,t);
	Ir=diag(ones(r,1));
	Zr=zeros(r,1);
    Zf=zeros(fx,1);
	for i=1:t
		Bi=diag(B(:,i));
		P=[Bi*L; sqrt(lambda)*Ir; F*L];
		Q=[S(:,i); Zr; Zf];
		Y(:,i)=(P'*P)\(P'*Q);
    end
    
% fit [LY; sqrt(lambda)*Y]=[S;0], standard liner squares
function Y = myInverseSimple(A,B,M)
  [m,n]=size(B);
  k = size(A,2);
 Y=zeros(k,n);
  for i=1:n
	Mi=diag(M(:,i));
	Y(:,i) = (A'*Mi*Mi*A)\(A'*Mi*B(:,i));
  end

function Y = liliInverse(A,B,M)

  [m,n] = size(B);
  k = size(A,2);

  % this is much more efficient than unique(M','rows') and
  % almost never fails
  [Mu,I,J] = unique(rand(1,m)*M);
  
  % initialize to 0
  Y = zeros(k,n);
  for i = 1:length(I)
    ind = find(M(:,I(i)));
    if (~isempty(ind))
      Ai = A(ind,:);
      cols = find(J == i);
      Bi = B(ind,cols);
      Y(:,cols) = (Ai'*Ai)\(Ai'*Bi);
    end
  end
