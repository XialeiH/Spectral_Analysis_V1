function [l6Clamped,l6Unclamped] = l6ns_l6_dynamic(state,context,sourceScale)
% Dynamic L6 input generated from one state, including the canonical clamp.

if nargin<3 || isempty(sourceScale); sourceScale=[1 1]; end
state=state(:);
n=numel(state)/3;
s=state(1:n); c=state(n+(1:n));
fixed=context.FixedPointStruct;
s=fixed.S+sourceScale(1)*(s-fixed.S);
c=fixed.C+sourceScale(2)*(c-fixed.C);
if context.Isaturation
    eRaw=(1-context.CWeight)*s+context.CWeight*c;
    eBase=L6Convert(eRaw,context.EKpUse);
    adjustment=eBase./eRaw;
    adjustment(~isfinite(adjustment))=1;
    s=s.*adjustment;
    c=c.*adjustment;
end
e=(1-context.CWeight)*s+context.CWeight*c;
field=reshape(e,context.MapSize);
padded=padarray(field,[1 1],'circular');
convolved=conv2(padded,context.L6Kernel,'same');
l6Unclamped=L6Convert(convolved(2:end-1,2:end-1),context.L6Parameters);
l6Unclamped=l6Unclamped(:)/3;
l6Clamped=min(max(l6Unclamped,1),40);
end
