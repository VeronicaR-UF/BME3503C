function best_slice = SliceDecider(fileroot, num_slices, FourD, png_dir)

cd(png_dir)
file_name_list = [];

path_root = 'png\';

if FourD == 'y'
    Area_Mat = zeros(1, num_slices*2);
    for ii = 1:1:num_slices*2
        if ii <= num_slices
            if ii < 10
                pathway_to_png = strcat(path_root, fileroot, "_t001_z00", int2str(ii), '.png');
            elseif ii <100
                pathway_to_png = strcat(path_root, fileroot, "_t001_z0", int2str(ii), '.png');
            else
                pathway_to_png = strcat(path_root, fileroot, "_t001_z", int2str(ii), '.png');
            end
        else
            if ii-num_slices < 10
                pathway_to_png = strcat(path_root, fileroot, "_t003_z00", int2str(ii-num_slices), '.png');
            elseif ii-num_slices <100
                pathway_to_png = strcat(path_root, fileroot, "_t003_z0", int2str(ii-num_slices), '.png');
            else
                pathway_to_png = strcat(path_root, fileroot, "_t003_z", int2str(ii-num_slices), '.png');
            end
        end
        Brain_PNG = im2double(imread(pathway_to_png));
        % Bad_Seg = imadjust(Brain_PNG, [0.3 0.7], []);
        Bad_Seg = imbinarize(Brain_PNG, 0.5);
        Area_Mat(1, ii) = bwarea(Bad_Seg);
        file_name_list = [file_name_list pathway_to_png];
    end
    [~, index_ma] = max(Area_Mat);
    best_slice = file_name_list(index_ma);
   
else
    Area_Mat = zeros(1, num_slices);
    for ii = 1:1:num_slices
        if ii <= num_slices
            if ii < 10
                pathway_to_png = strcat(path_root, fileroot, "_z00", int2str(ii), '.png');
            elseif ii <100
                pathway_to_png = strcat(path_root, fileroot, "_z0", int2str(ii), '.png');
            else
                pathway_to_png = strcat(path_root, fileroot, "_z", int2str(ii), '.png');
            end
        end
        Brain_PNG = im2double(imread(pathway_to_png));
        Bad_Seg = imbinarize(Brain_PNG, 0.5);
        Area_Mat(1, ii) = bwarea(Bad_Seg);
        file_name_list = [file_name_list pathway_to_png];
    end
    [~, index_ma] = max(Area_Mat);
    best_slice = file_name_list(index_ma);
end