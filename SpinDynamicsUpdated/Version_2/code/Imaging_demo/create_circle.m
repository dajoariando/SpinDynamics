function circlePixels = create_circle (xlen, ylen, radius)

% create circle image
    [columnsInImage, rowsInImage] = meshgrid(1:xlen, 1:ylen);
    % Next create the circle in the image.
    centerX = xlen/2;
    centerY = ylen/2;

    circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
    %image(circlePixels) ;
    %colormap([0 0 0; 1 1 1]);
    %title('Binary image of a circle');
    
end