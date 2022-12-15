function [SF,chpt,estsigma]=automated_step_detection(trace,varargin)
%input arguments
%first: position trace with steps
%second: width of moving median filter (this removes short back and forth
%stepping
%third: step size filter (in a single run all steps below the threshold are
%removed and the step function is recalculated, note that this can again
%generate new steps with size below the threshold which are not removed)

%align trace dimension for evaluation
sztr=size(trace);
if sztr(1)~=1
    trace=trace';
end


%default values for filter values if no user input is definded
filter_w=1;
min_step_h=0;

%variable user input
if ~isempty(varargin)
    switch length(varargin)
        case 1
            filter_w=varargin{1};
        case 2
            filter_w=varargin{1};
            min_step_h=varargin{2};
        otherwise
            error(['To many input arguments'])
    end
end

%estimated trace std
estsigma=std(trace(2:end)-trace(1:end-1))./sqrt(2);

%perform step fit
[TF,SF]=ischange(trace,'Threshold',12*estsigma^2);

%correct for short for- and backsteps
cSF=SF;
if filter_w>1
    cSF=movmedian(SF,filter_w);
    cSF(1:filter_w)=SF(1:filter_w);
    cSF(end-filter_w+1:end)=SF(end-filter_w+1:end);
end

%remove steps smaller min_step_h
ccSF=cSF;
if min_step_h>0
    level=unique(cSF,'stable');
    steps=abs(level(2:end)-level(1:end-1));
    if any(steps<min_step_h)
        l=find(steps<min_step_h);
        for m=l
            newlevel=mean(ccSF(ccSF==level(m) | ccSF==level(m+1)));
            ccSF(ccSF==level(m) | ccSF==level(m+1))=newlevel;
            level(m:m+1)=newlevel;
        end
    end
end

%generate output and adjust to input dimensions
if sztr(1)~=1
    SF=ccSF';
else
    SF=ccSF;
end
chpt=find((SF(2:end)-SF(1:end-1))~=0)+1;
end
