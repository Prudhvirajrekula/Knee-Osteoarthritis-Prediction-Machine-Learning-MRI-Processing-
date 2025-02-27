function analyze_bone_masks_and_prepare_features()
    % Set paths
    folder_path = "/Users/ariperez/Desktop/Ari\'s\ Stuff/Pace/Final\ Project/preds"; % Change to your dataset folder path
    label_file  = "/Users/ariperez/Desktop/Ari\'s\ Stuff/Pace/Final\ Project/KL grade for baseline year and 2-year follow up.csv"; % Path to labels file
    output_features_csv = "/Users/ariperez/Desktop/Ari\'s\ Stuff/Pace/Final\ Project/new_distance_features.csv"; % Output feature file
    
    % Get all files in the folder
    files = dir(fullfile(folder_path, '*_pred.png'));
    
    % Extract patient IDs and slice numbers
    patient_slice_map = [];
    for i = 1:length(files)
        parts = regexp(files(i).name, '(\d+)_(\d+)_pred', 'tokens');
        if ~isempty(parts)
            patient_id = str2double(parts{1}{1});
            slice_number = str2double(parts{1}{2});
            patient_slice_map = [patient_slice_map; patient_id, slice_number]; %#ok<AGROW>
        end
    end
    
    % Get unique patient IDs
    patient_ids = unique(patient_slice_map(:, 1));
    % patient_ids = [9002116] % Used if we want to see only specific patients
    % startingSlice = 42;
    % endingSlice = 117;
    startingSlice = 27;
    imageSet = 111;
    gapStart = 0;
    gapEnd = 0;
    gap = gapEnd - gapStart + 1;
    interval = 1;
    num_slices = (imageSet + 1 - gap)/interval; % Fixed length for all feature vectors (27–137 slices). Removing slices. Every nth slice
    
    % Initialize features matrix
    features = NaN(length(patient_ids), num_slices);  % Initialize with NaNs to handle missing slices
    patient_indices = zeros(length(patient_ids), 1);
    
    % Process each patient
    for p = 1:length(patient_ids)
    % for p = 1:1
        patient_id = patient_ids(p);
        patient_indices(p) = patient_id;

        endingSlice = startingSlice + imageSet - 1;
        
        % Get slice numbers for the patient
        slice_numbers = patient_slice_map(patient_slice_map(:, 1) == patient_id, 2);
        valid_slices = slice_numbers(slice_numbers >= startingSlice & slice_numbers <= endingSlice);

        % Remove slices in gap
        valid_slices = valid_slices(~ismember(valid_slices, gapStart:gapEnd));
        
        % Initialize distances for this patient
        distances = NaN(1, num_slices);

        slice_interval = 1;
        slice_idx = 1;

        % Process each valid slice
        for slice_number = valid_slices'
            % Used when we want to skip slices
            if interval > 1 && mod(slice_interval, interval) != 0
                slice_interval = slice_interval + 1;
                continue; 
            end

            slice_interval = 1;
            file_name = sprintf('%d_%03d_pred.png', patient_id, slice_number);
            file_path = fullfile(folder_path, file_name);
            
            if exist(file_path, 'file') ~= 2
                continue;
            end

            % Use this to diagnose issues
            % fmt=['Processing: %s \n'];
            % fprintf(fmt, file_path);

            slice_image = imread(file_path);

            median_min_distance = findMedianOfSmallestDistances(slice_image, 3);
            distance_between_centroids = measure_bone_distance(slice_image);

            % Use this to diagnose issues
            % fmt=['min distance: %d \n'];
            % fprintf(fmt, min_distance);
            % Use this to diagnose issues
            % fmt=['centroid distance: %d \n'];
            % fprintf(fmt, distance_between_centroids);

            % Use this to diagnose issues
            % fmt=['Distance ratio: %d \n'];
            % fprintf(fmt, distance_ratio);

            if ~isnan(median_min_distance) && ~isnan(distance_between_centroids)
                distance_ratio = median_min_distance/distance_between_centroids * 100; % This number is pretty small so multiply by 100 for more digits
                distances(slice_idx) = distance_ratio;
                slice_idx = slice_idx + 1;
            else
              fmt=['%s is NAN \n'];
              fprintf(fmt, file_path);
            end
        end
        
        % Interpolate NaN values to ensure a complete feature vector
        nan_indices = isnan(distances);

        % Use this to diagnose issues
        % fmt=['Patient %d has %d distances =' repmat(' %1.0f ',1,numel(distances)) '\n'];
        % fprintf(fmt, patient_id, length(distances), distances);

        if all(nan_indices)
            fprintf('Patient %d has no valid slices.\n', patient_id);
            continue;
        end

        distances(nan_indices) = interp1(find(~nan_indices), distances(~nan_indices), find(nan_indices), 'linear', 'extrap');

        % Use for diagnosing issues
        % fprintf('Feature length:  %d.\n', length(features(p, :)));
        % fprintf('Distance length:  %d.\n', length(distances));

        features(p, :) = distances;


        % Used when we are generating images
        % pause
        %
        % close all
    end
    
    % Load KL grade labels
    labels = read_labels(label_file, patient_indices);
    
    % Remove rows with NaN in features or labels
    valid_rows = ~any(isnan(features), 2) & ~isnan(labels);
    features = features(valid_rows, :);
    patient_indices = patient_indices(valid_rows);
    labels = labels(valid_rows);
    
    % Combine features and labels
    dataset = [patient_indices, features, labels];
    
    % Save features and labels to a CSV file
    fid = fopen(output_features_csv, 'w');

    % Writeheader with appropriate column names
    slice_headers = arrayfun(@(x) sprintf('Slice%d', x), 1:num_slices, 'UniformOutput', false);
    headers = [{'PatientID'}, slice_headers, {'KLGrade'}];

    % Manually join the headers using sprintf
    header_line = sprintf('%s', headers{1});
    for k = 2:length(headers)
        header_line = sprintf('%s,%s', header_line, headers{k});
    end
    fprintf(fid, '%s\n', header_line);

    % Write data: Use a loop to write each patient's data
    for i = 1:size(dataset, 1)
        % Extract the row data for the current patient
        row_data = dataset(i, :);

        % Format the row properly: PatientID, Slice features, KLGrade
        fprintf(fid, '%d', row_data(1)); % Write PatientID first

        % Write slice data: Slice27 to Slice137
        for j = 2:length(row_data)-1
            fprintf(fid, ',%.2f', row_data(j)); % Write each slice value
        end

        % Write KLGrade in the last column without an extra comma at the end
        fprintf(fid, ',%d\n', row_data(end));
    end

    fclose(fid);
    fprintf('Feature dataset saved to %s\n', output_features_csv);

end

function labels = read_labels(label_file, patient_ids)
    % Open the CSV file
    fid = fopen(label_file, 'r');
    if fid == -1
        error('Could not open the file: %s', label_file);
    end
    
    % Read the header line to identify columns
    header_line = fgetl(fid);
    headers = regexp(header_line, ',', 'split'); % Use regexp to split by commas
    
    % Find the indices for the 'cases' and 'V00XRKL' columns
    case_idx = find(strcmp(headers, 'cases'));
    kl_grade_idx = find(strcmp(headers, 'V00XRKL'));
    
    if isempty(case_idx) || isempty(kl_grade_idx)
        error('The required columns "cases" and "V00XRKL" are missing from the file.');
    end
    
    % Read the rest of the data
    data = textscan(fid, repmat('%s', 1, numel(headers)), 'Delimiter', ',', 'HeaderLines', 0);
    fclose(fid);
    
    % Extract the cases and KL grades
    cases = str2double(data{case_idx});
    kl_grades = str2double(data{kl_grade_idx});
    
    % Match KL grades to patient IDs
    labels = NaN(size(patient_ids));
    for i = 1:length(patient_ids)
        match_idx = find(cases == patient_ids(i), 1);
        if ~isempty(match_idx)
            labels(i) = kl_grades(match_idx);
        end
    end
end

function distance = measure_bone_distance(slice_image)
    % Threshold the image to separate the bones
    binary_image = slice_image > 0; % Assuming non-zero pixels represent bones
    
    % Label connected regions
    labeled_image = bwlabel(binary_image);
    stats = regionprops(labeled_image, 'Centroid', 'BoundingBox');
    
    % Ensure there are at least two regions for femur and tibia
    if numel(stats) < 2
        distance = NaN; % Unable to calculate
        return;
    end
    
    % Sort regions by area to identify femur and tibia
    [~, sorted_idx] = sort(cellfun(@(x) prod(x(3:4)), {stats.BoundingBox}), 'descend');
    femur = stats(sorted_idx(1)).Centroid;
    tibia = stats(sorted_idx(2)).Centroid;

    % Calculate the distance between the centroids
    distance = norm(femur - tibia);

    % Use for diagnosing issues
    % distance
    return

    figure;

    imshow(labeled_image, []);
    hold on

    %%Plot Bounding Box
    for n=1:size(stats,1)
        rectangle('Position',stats(n).BoundingBox,'EdgeColor','g','LineWidth',2)
        s = regionprops(labeled_image,'centroid');
        centroids = cat(1,s.Centroid);

        hold on
        plot(centroids(:,1),centroids(:,2), 'b*')
    end

    plot([femur(1) tibia(1)], [femur(2) tibia(2)]);
    hold off
    title('Centroid plot');
    pause

end

function x = findMedianOfSmallestDistances(image, min_distance_count)

  % Finds the boundary objects in the image 
  image_no_border = imclearborder(image);
  [B,L8] = bwboundaries(image_no_border,8);
  objectAreas = [];

  for object = 1:length(B);
    boundaryObject = B{object};
    boundaryx = boundaryObject(:, 2);
    boundaryy = boundaryObject(:, 1);
    objectArea = polyarea(boundaryx, boundaryy);
    objectAreas(end+1) = objectArea;
  end

  if length(objectAreas) < 2
      x = NaN; % Unable to calculate
      return;
  end

  [mx,biggestIndex] = max(objectAreas);
  objectAreas(biggestIndex) = -1;
  [mx,secBiggestIndex] = max(objectAreas);

  boundary1 = B{biggestIndex};
  boundary2 = B{secBiggestIndex};

  % Below we take the 2 boundary objects and find the distances between all of the pixels in the 2 objects using pdist2
  distances = pdist2(boundary1, boundary2);

  % Find the minimum distance
  minDistance = min(distances(:));

  minDistances = sort(unique(distances(:)))(1:min_distance_count);
  x = median(minDistances); % Get the median of the distances
  return % we return here normally because we only want the data. We comment out this line when we want to see the image output

  % Save the coordinates of minimum distance lines
  [rows, cols] = find(distances == minDistance);

  figure;

  % subplot(1,3,1);
  % imshow(image, []);
  % title('Original Image');
  %
  % % Blots the boundary objects in red
  % subplot(1,3,2);
  % imshow(image, [])
  % hold on
  % plot(boundary1(:,2), boundary1(:,1), 'green', 'LineWidth', 1)
  % plot(boundary2(:,2), boundary2(:,1), 'red', 'LineWidth', 1)
  % hold off
  % title('L8 Boundary Image');
  %
  % subplot(1,3,3);
  imshow(image, [])
  hold on 

  plot(boundary1(:,2), boundary1(:,1), 'green', 'LineWidth', 1)
  plot(boundary2(:,2), boundary2(:,1), 'green', 'LineWidth', 1)

  hold on 

  boundary1x = boundary1(:, 2);
  boundary1y = boundary1(:, 1);
  boundary2x = boundary2(:, 2);
  boundary2y = boundary2(:, 1);

  % Finds all of the shortest distances and plots them in blue
  for k = 1:length(rows)
    index1 = rows(k);
    index2 = cols(k);
    xPlot = [boundary1x(index1), boundary2x(index2)];
    yPlot = [boundary1y(index1), boundary2y(index2)];
    plot(xPlot, yPlot, 'blue', 'LineWidth', 2);
  end   
  hold off
  title('Minimum Distances between Objects');
  pause
end
