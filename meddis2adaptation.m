% 
% Derive adaptation characteristcs from parameter values in meddis model
%
% After determine that our simplied version of meddis model can fit the original
% version very well, we come to get the parameters of simplified version from
% Synapse Adaptation Characteristics, And then explore how the change of k will
% affect the time constant of the adaptation (like in Westerman_001.m)
%
% derive the synapse adaptation characteristics from
% parameters in Meddis paper
%
% Waveform 1 : simulation
% Waveform 2 : Analytical results
% Waveform 3 : Analytical Results based on the simplified value 
%              (assume l+r -> inf, and rate = k*q
%
% Waveform 4 : Parameter Estimated from the desired PST characteristics
%
%

% From Sumnar et al., 2002
% fiber_types is specified by different M in meddis model
fiber_types = [10 13 8];

% We specify the spontatneous rate of these fibers to determine k value at rest
fiber_types_Asp = [60 10 0.1];

% We specify the sustained rate to determine k value at high stimulus level
Meddis2002.Ass = 350;  % specify the sustained rate for these fibers (used to determine k value at high level

Meddis2002.y = 10;
Meddis2002.x = 66.3;
Meddis2002.l = 2580;
Meddis2002.r = 6580;
Meddis2002.u = Meddis2002.r/(Meddis2002.r+Meddis2002.l);

for i = 1:length(fiber_types)

    model = Meddis2002;
    model.M = fiber_types(i);
    model.Asp = fiber_types_Asp(i);
    
    y = model.y;
    x = model.x;
    u = model.u;
    M = model.M;

    Asp = model.Asp;
    Ass = model.Ass;
    
    % Based on Eq. 22 (Zhang and Carney, 2005), maximum response rate of model outptu when k goes to infinite
    Rmax = y*M/(1-u);
    if (Rmax < (Ass+10))
        Ass = Rmax-10;
        model.Ass = Ass;
    end;
    % Use another Ass
    Ass = Rmax-10;
    model.Ass = Ass;
    
    if Ass < Asp
        Asp = Ass;
        model.Asp = Asp;
    end;
    
    % k value is derived based Eq. 24 in Zhang and Carney, 2005
    % R = kyM/[y+k(1-u)];
    k1 = y*Asp/(y*M-(1-u)*Asp);
    k2 = y*Ass/(y*M-(1-u)*Ass);

    model.k1 = k1;
    model.k2 = k2;
    
    % adaptation components 
    % -1/tR, -1/tST are roots of of denominator of Q(s) in Eq. 17
    % Ar = k2*Q(s)(s+1/tR), s = -1/tR
    % Ast = k2*Q(s)(s+1/tST), s = -1/tST

    % polynomial function of denominator of Q(s) in Eq. 17
    a = [1 x+y+k2 x*(y+k2*(1-u))];
    rv = roots(a);
    
    tauValues = -1./rv;
    if tauValues(1) < tauValues(2)
        tauST = tauValues(2);
        tauR = tauValues(1);
    else
        tauST = tauValues(1);
        tauR = tauValues(2);
    end;
    
    rootR = -1/tauR;
    rootST = -1/tauST;
    
    q0_ = Asp/k1;
    w0_ = Asp*u/x;      % Eq. 15
    % based on Eq. 17
    qAr = ((rootR*q0_+y*M)*(rootR+x) + w0_*x*rootR)/(rootR-rootST)/rootR;
    Ar = k2*qAr;
    
    qAst = ((rootST*q0_+y*M)*(rootST+x) + w0_*x*rootST)/(rootST-rootR)/rootST;
    Ast = k2*qAst
    
    adaptation.tauST_infinitek = 1/(x-u*x);       
    adaptation.tauST = tauST;
    adaptation.tauR = tauR;
    adaptation.Asp = Asp;
    adaptation.Ass = Ass;
    adaptation.Ar = Ar;
    adaptation.Ast = Ast;

    adaptation
    model
end;
