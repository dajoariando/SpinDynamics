function [Xq,Yq,Uq] = get_interp_field_plane_data (file, xgridsize, ygridsize)
    %% load the file 
    dat = import_field_plane_COMSOL(file);

    %% reorganize the data
    x = dat(:,1);
    y = dat(:,2);
    u = dat(:,3);

    %% perform scattered interpolant
    F = scatteredInterpolant(x,y,u,'linear','linear');

    %% make interpolated grids
    xgrid = linspace(min(x),max(x),xgridsize);
    ygrid = linspace(min(y),max(y),ygridsize);
    [Xq,Yq] = meshgrid (xgrid,ygrid);

    %% calculate data from grid
    Uq = F(Xq,Yq);

    %% plot
    %figure;
    %mesh (Xq,Yq,Uq);

end