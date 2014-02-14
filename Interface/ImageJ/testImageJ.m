
initialize()

%% directly read a stack

stack  = imread_hdr('./Test/Images/Stack/StackIsoCropM');

% 
% figure(1)
% clf
% implot3d(stack)
% 
% figure(2)
% hist(log2(double(stack(:))+eps), 255)

%% imagej

universe = ijplot3d(stack);

%% some usefull commands

universe.contains('Figure')
%universe.resetView()
%universe.adjustView()
str = universe.allContentsString();
%methodsview(universe)


cont = universe.getContent('Figure');
methodsview(cont)

col = javax.vecmath.Color3f(1.0, 0, 0)
cont.setColor(col)


%% others

vol = cont.getContent()
methodsview(vol)


rend = vol.getRenderer()
methodsview(rend)


%% other stuff

bg = rend.getVolumeNode()
bg.numChildren()

sw = bg.getChild(0)
sw.numChildren()

gr = sw.getChild(4)  % set of surfaces
gr.numChildren()

bgs = gr.getChild(10)
bgs.numChildren()

shape3d = bgs.getChild(0)

app = shape3d.getAppearance()
geo = shape3d.getGeometry()

tex = app.getTexture()


%% change properties

cont = universe.getContent('Figure');
vol = cont.getContent();
bg = vol.getChild(0);
sw = bg.getChild(0);
og = sw.currentChild();
nplanes = og.numChildren;


for id = 0:nplanes-1

   sh = og.getChild(id).getChild(0);
   ap = sh.getAppearance();
   %tex = ap.getTextureUnitState(0).getTexture
   %texatrr = ap.getTextureUnitState(0).getTextureAttributes
   
   % set alpha value
   re = ap.getRenderingAttributes;
   re.setAlphaTestValue(0.3);
   
   % set color shading 
   %co = ap.getColoringAttributes;
   %co.setCapability(co.ALLOW_SHADE_MODEL_WRITE)
   %co.setShadeModel(co.SHADE_GOURAUD)
   
   % set transparency
   tr = ap.getTransparencyAttributes;
   tr.setTransparency(0.);
   
   
end
   



%% get alpha values of texture

cont = universe.getContent('Figure');
vol = cont.getContent();
bg = vol.getChild(0);
sw = bg.getChild(0);
og = sw.currentChild();
id = 50;

sh = og.getChild(id).getChild(0);
ap = sh.getAppearance();
tex = ap.getTextureUnitState(0).getTexture
texatrr = ap.getTextureUnitState(0).getTextureAttributes
img = tex.getImage(0).getImage
dat = img.getData;

al = img.getAlphaRaster


