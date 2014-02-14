function array = getPixelArray(javaimage)

w = javaimage.getWidth();
h = javaimage.getHeight();

array = zeros(w, h);
for j = 0:w-1
    for k = 0:h-1
        %array(j+1,k+1) = javaimage.getRGB(j, k);
        array(j+1,k+1) = javaimage.getPixel(j,k,[]);
    end
end

end