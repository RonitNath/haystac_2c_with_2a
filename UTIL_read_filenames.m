function filenames = getFilenames(directory)
    try
        if exist(directory, "dir") ~= 7
            error("Directory does not exist");
        end

        files = dir(directory);

        filenames = {};

        for i = 1:length(files)
            if ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')
                filenames = [filenames; files(i).name];
            end
        end
    catch exception
        disp(['Error: ', exception.message]);

        filenames = {};
    end
end