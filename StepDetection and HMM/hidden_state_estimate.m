function output = hidden_state_estimate(steps,chpt,L,stsz)
% transition states assignment 1=initial dummy state to assure all possible
% initial states / 2=bound to bound / 3=bound to unbound / 4=unbound to bound / 5=bound to
% bound (reduced stepsize)/ 6=unbound to unbound (enhanced stepsize)

xemi=1:1:100; % possible emission values (step sizes)
w=1; % width of step distribution for each transition

% transition matrix between states
TRANS=[0 .33 .33 .33 0 0;
    0 .495 .495 0 .01 0;
    0 0 0 .99 0 .01;
    0 .4995 .4995 0 .001 0;
    0 .4995 .4995 0 .001 0;
    0 0 0 .99 0 .01];

% propability of emission values (stepsizes), only discrete values possible
EMIS=[0.*xemi;
    1/sqrt(2*pi*2^2).*exp(-(xemi-round((xemi)./(2*stsz)).*(2*stsz)).^2./(2*w^2));
    1/sqrt(2*pi*2^2).*exp(-(xemi-stsz-round((xemi-stsz)./(2*stsz)).*(2*stsz)).^2./(2*w^2));
    1/sqrt(2*pi*2^2).*exp(-(xemi-stsz-round((xemi-stsz)./(2*stsz)).*(2*stsz)).^2./(2*w^2));
    1/sqrt(2*pi*2^2).*exp(-(xemi-stsz-round((xemi-stsz)./(2*stsz)).*(2*stsz)).^2./(2*w^2));
    1/sqrt(2*pi*2^2).*exp(-(xemi-round((xemi)./(2*stsz)).*(2*stsz)).^2./(2*w^2))];

% translate transition states into 1HB and 2HB
statei=[0 2 2 1 2 1];
statef=[0 2 1 2 2 1];


% perform HMM estimate
HMMstate=hmmviterbi(steps,TRANS,EMIS);

% generate state array of length L changing states at chpt
eststate=zeros(L,1);
for k=1:length(chpt)
    switch k
        case 1
            eststate(1:chpt(k)-1)=statei(HMMstate(k));
            eststate(chpt(k):chpt(k+1))=max(statef(HMMstate(k)),statei(HMMstate(k+1)));
        case length(chpt)
            eststate(chpt(k):end)=statef(HMMstate(k));
        otherwise
            eststate(chpt(k):chpt(k+1))=max(statef(HMMstate(k)),statei(HMMstate(k+1)));
    end
end

% generate output
output.HMMest=eststate;
output.transitions=HMMstate;

end