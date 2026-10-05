function gke_checkcode(codeRoot)
files=dir(fullfile(codeRoot,'*.m'));
failed=false;
for index=1:numel(files)
    issues=checkcode(fullfile(files(index).folder,files(index).name),'-id');
    if ~isempty(issues)
        fprintf('FILE %s\n',files(index).name);
        disp(struct2table(issues));
        failed=true;
    end
end
assert(~failed,'MATLAB checkcode issues found.');
fprintf('MATLAB checkcode passed for %d files.\n',numel(files));
end
