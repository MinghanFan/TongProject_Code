function varargout = su_lc_tools(action, varargin)
%SU_LC_TOOLS Shared functions for finding and exporting Su et al. LC sections.
% Keep this file in the same folder as find_Su_LC_slices.m and
% export_LC_three_sets.m.  Run those two scripts, not this helper.
%
% The XYZ -> QuickNII conversion and plane-mesh intersection are the same
% ones that produced the previously checked CC1A_163_A5_1 Set 3 overlay.

switch lower(char(action))
    case 'loadmesh'
        [varargout{1},varargout{2}] = loadQuickMesh(varargin{:});
    case 'loadrecords'
        varargout{1} = loadRecords(varargin{:});
    case 'findrecord'
        varargout{1} = findRecord(varargin{:});
    case 'applyrecord'
        varargout{1} = applyRecord(varargin{:});
    case 'slicemask'
        [varargout{1},varargout{2}] = sectionMask(varargin{:});
    case 'saveimages'
        saveImages(varargin{:});
    otherwise
        error('Unknown su_lc_tools action: %s', string(action));
end
end

function records = loadRecords(jsonPath)
    assert(isfile(jsonPath),'Cannot find VisuAlign JSON: %s',jsonPath);
    % Retain the path-separator handling used in the tested scripts.
    J = jsondecode(strrep(fileread(jsonPath),'\','/'));
    assert(isfield(J,'slices') && ~isempty(J.slices), ...
        'No slices found in JSON: %s',jsonPath);
    records = J.slices;
end

function record = findRecord(records,sliceName)
    filenames = string({records.filename});
    hit = find(contains(filenames,string(sliceName)));
    assert(numel(hit)==1, ...
        'Expected one JSON record for %s, found %d.',string(sliceName),numel(hit));
    record = records(hit);
end

function s = applyRecord(s,record)
    a = double(record.anchoring(:)');
    assert(numel(a)==9 && all(isfinite(a)), ...
        'QuickNII record requires nine finite anchoring parameters.');
    s.o = a(1:3); s.u = a(4:6); s.v = a(7:9);
    if isfield(record,'markers') && ~isempty(record.markers)
        s.markers = double(record.markers);
    else
        s.markers = zeros(0,4);
    end
    assert(size(s.markers,2)==4,'VisuAlign markers must have four columns.');
end

function [Q,F] = loadQuickMesh(meshPath,atlasSize,voxelSizeUm,midlineUm)
    assert(isfile(meshPath),'Cannot find Su core mesh: %s',meshPath);
    assert(numel(atlasSize)==3,'Atlas size must be a three-element vector.');
    if exist('readSurfaceMesh','file')
        mesh = readSurfaceMesh(meshPath);
        V = double(mesh.Vertices);
        F = double(mesh.Faces);
    else
        [V,F] = readSimpleOBJ(meshPath);
    end
    assert(size(V,2)==3 && size(F,2)==3,'Expected a triangular OBJ mesh.');
    assert(all(isfinite(V(:))) && all(F(:)>=1) && all(F(:)<=size(V,1)), ...
        'OBJ has invalid vertices or face indices.');

    % Su meshes were generated after pooling left/right hemisphere cells
    % into the half-space x<=5700 um. Mirror the core across that midline.
    Vr = V;
    Vr(:,1) = 2*midlineUm - Vr(:,1);
    F = [F; F+size(V,1)];
    V = [V; Vr];

    % VALIDATED FOR THIS PROJECT'S ATLAS PREPROCESSING:
    % A = flip(flip(permute(annotationVolume,[3 1 2]),2),3).
    % Slice.getObliqueSlice samples round(o+a*u+b*v)+[0 1 1].
    % Convert Su physical x/y/z (um) to that continuous QuickNII frame.
    Q = [1+V(:,1)/voxelSizeUm, ...
         atlasSize(2)-1-V(:,2)/voxelSizeUm, ...
         atlasSize(3)-1-V(:,3)/voxelSizeUm];
end

function [mask,nSegments] = sectionMask(Q,F,o,u,v,W,H)
    segments = intersectMeshWithPlane(Q,F,o,u,v,W,H);
    nSegments = size(segments,1);
    mask = rasterizeMeshSection(segments,H,W);
end

function seg = intersectMeshWithPlane(Q,F,o,u,v,W,H)
    % Exact triangle/plane crossings in the unwarped QuickNII plane.
    % Output each segment as [x1 y1 x2 y2] in thumbnail pixel space.
    o = double(o(:)'); u = double(u(:)'); v = double(v(:)');
    n = cross(u,v);
    assert(norm(n)>0,'Degenerate QuickNII anchoring vectors.');
    n = n/norm(n);
    d = (Q-o)*n(:);
    td = reshape(d(F(:)),size(F));
    chosen = min(td,[],2)<=0 & max(td,[],2)>0;
    Fc = F(chosen,:);
    if isempty(Fc); seg = zeros(0,4); return; end

    entries = zeros(0,3);
    basis = [u(:),v(:)];
    for e = 1:3
        a = Fc(:,e);
        b = Fc(:,mod(e,3)+1);
        da = d(a); db = d(b);
        crossEdge = (da<=0 & db>0) | (db<=0 & da>0);
        if ~any(crossEdge); continue; end
        triIdx = find(crossEdge);
        a = a(crossEdge); b = b(crossEdge);
        t = da(crossEdge)./(da(crossEdge)-db(crossEdge));
        pts = Q(a,:) + t.*(Q(b,:)-Q(a,:));
        coeff = basis\(pts-o).';
        xy = [coeff(1,:)'*W+1,coeff(2,:)'*H+1];
        entries = [entries;[triIdx,xy]]; %#ok<AGROW>
    end
    if isempty(entries); seg = zeros(0,4); return; end
    entries = sortrows(entries,1);
    numPerTriangle = accumarray(entries(:,1),1,[size(Fc,1),1]);
    keep = numPerTriangle(entries(:,1))==2;
    entries = entries(keep,:);
    if isempty(entries); seg = zeros(0,4); return; end
    assert(mod(size(entries,1),2)==0,'Unpaired mesh intersection entries.');
    seg = reshape(entries(:,2:3).',4,[]).';
    nonzero = hypot(seg(:,1)-seg(:,3),seg(:,2)-seg(:,4))>1e-8;
    seg = seg(nonzero,:);
end

function mask = rasterizeMeshSection(seg,H,W)
    % Even-odd fill of intersections with horizontal rows. Multiple LC
    % components and disconnected contours are handled by crossing parity.
    mask = false(H,W);
    if isempty(seg); return; end
    low = max(1,ceil(min(min(seg(:,[2 4])))));
    high = min(H,floor(max(max(seg(:,[2 4])))));
    for r = low:high
        y0 = seg(:,2); y1 = seg(:,4);
        crosses = (y0<=r & y1>r) | (y1<=r & y0>r);
        if ~any(crosses); continue; end
        s = seg(crosses,:);
        xs = s(:,1)+(r-s(:,2)).*(s(:,3)-s(:,1))./(s(:,4)-s(:,2));
        xs = sort(xs);
        if mod(numel(xs),2)~=0
            warning('Odd mesh crossings in row %d; skipping this row.',r);
            continue;
        end
        for k = 1:2:numel(xs)
            first = max(1,ceil(xs(k)));
            last = min(W,floor(xs(k+1)));
            if last>=first
                mask(r,first:last) = true;
            end
        end
    end
end

function saveImages(histology,alignedAtlas,lcID,outDir,suMask)
    % Create histology.png, atlas.png, overlay.png using the same styles
    % as the tested original/manual and Su LC comparison scripts.
    if nargin<5; suMask=[]; end
    if ~isfolder(outDir); mkdir(outDir); end
    rgb=histology;
    if size(rgb,3)==1; rgb=repmat(rgb,[1 1 3]); end
    if ~isa(rgb,'uint8'); rgb=im2uint8(rgb); end
    [H,W,~] = size(rgb);
    assert(isequal(size(alignedAtlas),[H W]),'Aligned atlas size mismatch.');
    if ~isempty(suMask)
        assert(isequal(size(suMask),[H W]),'Su mask size mismatch.');
    end
    imwrite(rgb,fullfile(outDir,'histology.png'));

    borders = Slice.annotationBorders(Slice.rescaleAnnotation(alignedAtlas));
    oldLC = alignedAtlas==lcID;
    oldEdge = bwperim(oldLC,8);
    atlasView = uint8(25*ones(H,W,3));
    overlay = rgb;

    if isempty(suMask)
        % Set 1 / Set 2: same style as compare_LC_original_vs_manual.m.
        atlasView = setPixels(atlasView,oldLC,[110 30 30]);
        atlasView = setPixels(atlasView,borders,[215 215 215]);
        atlasView = setPixels(atlasView,oldEdge,[255 55 55]);
        overlay = blendPixels(overlay,borders,[255 255 255],0.50);
        overlay = blendPixels(overlay,oldEdge,[255 30 30],1.00);
    else
        % Set 3: same style as create_Set3_Su_LC.m.
        suMask = logical(suMask);
        suEdge = bwperim(suMask,8);
        atlasView = blendPixels(atlasView,suMask,[0 205 220],0.27);
        atlasView = blendPixels(atlasView,borders,[220 220 220],1.00);
        atlasView = blendPixels(atlasView,oldEdge,[255 55 55],1.00);
        atlasView = blendPixels(atlasView,suEdge,[0 255 255],1.00);
        overlay = blendPixels(overlay,suMask,[0 205 220],0.25);
        overlay = blendPixels(overlay,borders,[255 255 255],0.50);
        overlay = blendPixels(overlay,oldEdge,[255 45 45],1.00);
        overlay = blendPixels(overlay,suEdge,[0 255 255],1.00);
    end
    imwrite(atlasView,fullfile(outDir,'atlas.png'));
    imwrite(overlay,fullfile(outDir,'overlay.png'));
end

function img = setPixels(img,mask,color)
    for c = 1:3
        ch = img(:,:,c);
        ch(mask) = uint8(color(c));
        img(:,:,c) = ch;
    end
end

function img = blendPixels(img,mask,color,alpha)
    img = double(img);
    for c = 1:3
        ch = img(:,:,c);
        ch(mask) = (1-alpha)*ch(mask)+alpha*color(c);
        img(:,:,c) = ch;
    end
    img = uint8(round(img));
end

function [V,F] = readSimpleOBJ(meshPath)
    fid=fopen(meshPath,'r');
    assert(fid~=-1,'Unable to open OBJ: %s',meshPath);
    c=textscan(fid,'%s','Delimiter','\n','Whitespace','');
    fclose(fid);
    lines=c{1};
    vlines=lines(startsWith(lines,'v '));
    flines=lines(startsWith(lines,'f '));
    V=zeros(numel(vlines),3);
    F=zeros(numel(flines),3);
    for k=1:numel(vlines)
        vals=sscanf(vlines{k}(3:end),'%f');
        assert(numel(vals)>=3,'Invalid OBJ vertex.');
        V(k,:)=vals(1:3)';
    end
    for k=1:numel(flines)
        parts=regexp(strtrim(flines{k}(3:end)),'\s+','split');
        assert(numel(parts)==3,'OBJ contains a nontriangular face.');
        for j=1:3
            tok=regexp(parts{j},'^\d+','match','once');
            assert(~isempty(tok),'Invalid/negative OBJ face index.');
            F(k,j)=str2double(tok);
        end
    end
end
