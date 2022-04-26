
import cv2, sys
import numpy as np

class VideoGenerator:

    def __init__(self,VideoName,Particles):
        self.Name       = VideoName
        self.Particles  = Particles
        self.Width      = 2000
        self.Height     = 2000
        self.Channels   = 3
        self.framerate  = 30
        self.Padding    = 0.05
    
    def GetRanges(self): 
        xMin = self.Particles[0].Tra[0][0]
        xMax = self.Particles[0].Tra[0][0]
        yMin = self.Particles[0].Tra[1][0]
        yMax = self.Particles[0].Tra[1][0]
        for iFrame in range(len(self.Particles[0].Tra[0])):
            for particle in self.Particles:
                xMin = min(particle.Tra[0][iFrame],xMin)
                xMax = max(particle.Tra[0][iFrame],xMax)
                yMin = min(particle.Tra[1][iFrame],yMin)
                yMax = max(particle.Tra[1][iFrame],yMax)
        return [xMin,xMax] , [yMin,yMax] 

    def Generate(self):
        
        xRange , yRange = self.GetRanges()

        xRange[0] -= self.Padding*(xRange[1]-xRange[0])
        xRange[1] += self.Padding*(xRange[1]-xRange[0])
        yRange[0] -= self.Padding*(yRange[1]-yRange[0])
        yRange[1] += self.Padding*(yRange[1]-yRange[0])

        def GetPixel(X,Y):
            x = int((self.Width *float(X - xRange[0]))/(xRange[1]-xRange[0])//1)
            y = int((self.Height*float(Y - yRange[0]))/(yRange[1]-yRange[0])//1)
            return x,y


        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video  = cv2.VideoWriter(self.Name+".mp4", fourcc, float(self.framerate), (self.Width, self.Height))
    
        num_frames = len(self.Particles[0].Tra[0])

        print("Generating Video for:",self.Name)        
        for frame_count in range(num_frames):

            sys.stdout.write('\r')
            sys.stdout.write("%d%% : [%-1s] " % ( (100/num_frames)*frame_count,'='*int((float(50*frame_count)/num_frames)//1)))
            sys.stdout.flush()
            
            frame = np.zeros((self.Height,self.Width,self.Channels), dtype=np.uint8)
            for particle in self.Particles:
                x,y = GetPixel(particle.Tra[0][frame_count],particle.Tra[1][frame_count])
                for xC in range(10):
                    for yC in range(10):
                        for rbg in range(self.Channels):
                            frame[self.Height-(y+yC-5)][(x+xC-5)][rbg] = particle.Col[rbg]
            video.write(frame)
        print()
 
        video.release()