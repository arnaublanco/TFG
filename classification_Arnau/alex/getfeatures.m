function [features, selected]=getfeatures(featmat, grouping, groupi)
%Usage: [ selfeatmat selected]= getfeatures(featmat, grouping, groupi)
% This _function returns a matrix containg only the features group
%       indicated by groupi. Selected are the corresponing row indexes
%     featmat: all _set feature matrix
%     grouping: grouping matrix
%     groupi: selected group indexes
%
%     Alex Perera, April 2002


if isscalar(grouping) && grouping==0
    selected=groupi;
else
    selected= [];
    for ig=groupi,
        inous_indexos = grouping(ig,:) ~= 0 ;
        selected = [selected grouping(ig,inous_indexos) ];
    end
end

features = featmat(:,selected);