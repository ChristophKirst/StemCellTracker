function mij = ijput(image)
%
% mij = ijput(image)
%
% description: 
%     creates a new image in ImageJ instance
%
% input:
%     image    image to pass to ImageJ
%
% output:
%     mij      instance of ImageJ to which the image was passed
% 
% todo: take careof image types
%
% See also: ij get


format = imformat(image);

if strcmp(format, 'hw')



		else if (object instanceof float[][][]) {
			float[][][] fs = (float[][][]) object;
			int height = fs.length;
			int width = fs[0].length;
			int stackSize = fs[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				FloatProcessor floatprocessor = new FloatProcessor(width, height);
				float[] fp = (float[]) floatprocessor.getPixels();
				i = 0;
				int h = 0;
				while (h < height) {
					int w = 0;
					while (w < width) {
						fp[i] = fs[h][w][sz];
						w++;
						i++;
					}
					i = ++h * width;
				}
				if (sz == 0)
					floatprocessor.resetMinAndMax();
				imagestack.addSlice("", floatprocessor);
			}
			imp=new ImagePlus(title, imagestack);