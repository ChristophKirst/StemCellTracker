function compile

%% compiles code for ImageJ interface
jclasspath =  '/home/ckirst/programs/Fiji/jars';

system(['export CLASSPATH=' jclasspath]);
system(['javac -cp ' jclasspath filesep 'ij-1.48q.jar ' jclasspath   /home/ckirst/programs/Fiji/jarsMImageJ.java')

end


