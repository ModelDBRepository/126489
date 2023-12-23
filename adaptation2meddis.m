% 
% Derive synapse model parameters from adaptation characteristics for Meddis model
%

%%%%%%%%% Synapse Characteristics %%%%%%%%%%%%
% Defined in Zhang and Carney (2005)

fiber_types = [0.1 10 60];	% specified by spontaneous rate

Synapse.Asp = 50;
Synapse.PTS = 1+9*(Synapse.Asp/(9+Synapse.Asp));

Synapse.Ass = 350;
Synapse.tauR = 0.002;
Synapse.tauST = 0.060;
Synapse.Ar_over_Ast = 6.0;


for i = 1:length(fiber_types)
	Synapse.Asp = fiber_types(i);
	Synapse.PTS = 1+9*(Synapse.Asp/(9+Synapse.Asp)); % This will change as Asp changes?

	Asp         = Synapse.Asp;
	Ass         = Synapse.Ass;
	PTS         = Synapse.PTS;
	tauR        = Synapse.tauR;
	tauST       = Synapse.tauST;
	Ar_over_Ast = Synapse.Ar_over_Ast;

	Aon         = PTS*Ass;
	Ar          = (Aon-Ass)*(Ar_over_Ast)/(1.+Ar_over_Ast);
	Ast         = Aon-Ass-Ar;


	%/* Following equations is from X. Zhang and L.H. Carney (JASA, 2005) */

	Sr      = 1/tauR+1/tauST;
	Sr2     = Ar/tauR+Ast/tauST;
	Pr      = 1/(tauR*tauST);

	k2      = Sr2/(Aon-Asp);
	k1      = Asp/Aon*k2;

	beta    = (Ass-Asp)*k1*k2/(Asp*k2-Ass*k1);
	a1      = (beta+k2)*beta;
	a2      = -(Sr-k2)*(beta+k2);
	a3      = Pr;

	rv      = roots([a1 a2 a3]);
	z       = rv(2); % pick the positive one, if not maybe we should rethink it
	if((z>=1) | (z<=0))
	  error('Error : u must be 0-1');
	end;

	u = 1-z;
	y = beta*z;
	x = Sr-k2-y;
	m = Asp*(y+k1*z)/(y*k1);

	Meddis.x    = x;
	Meddis.y    = y;
	Meddis.u    = u;
	Meddis.m    = m;
	Meddis.k1   = k1;
	Meddis.k2   = k2;
	Meddis.l    = 2580; % This value is fixed;
	Meddis.r    = Meddis.l*u/(1-u);

Synapse
Meddis

end;
