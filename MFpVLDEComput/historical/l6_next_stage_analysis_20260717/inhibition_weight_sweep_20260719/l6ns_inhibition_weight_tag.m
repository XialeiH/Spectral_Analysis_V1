function tag = l6ns_inhibition_weight_tag(value)
% File-safe signed weight tag.

if value < 0
    prefix = 'm';
else
    prefix = '';
end
tag = [prefix strrep(sprintf('%.2f',abs(value)),'.','p')];
end
