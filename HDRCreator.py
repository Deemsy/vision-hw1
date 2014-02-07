import sys
import getopt
import cv2
import numpy as np
import math 
class HDRCreator:
    def __init__(self):
        
        //images
        self.images=[]
        
        //shutter speed
        self.speeds=[]

        //hdr image
        self.hdr=None
    def loadimages(self,filename):
        file=open(filename)
        for line in file:
            tokens=line.strip().split()
            if len(tokens)<2:
                continue
            filename=tokens[0]
            exposure=float(tokens[1])
            image=cv2.imread(filename)
            self.images.append(image)
            self.speeds.append(math.log(exposure))
    def createHDR(self):
        channels=[]
        for channel in range(3):

            #create pixel x num images array to solve linear system
            numPixels=70
            A=np.zeros((numPixels,len(self.images)),np.uint32)
            B=self.images
            
            #get sampled pixel location values in all images
            for imageNum in range(len(self.images)):
                img=cv2.split(self.images[imageNum])[channel] #get the correct channel
                for pixel in range(numPixels):
                    A[pixel,imageNum]=img[pixel*10,256]

            #solve the system
            g=self.solveSystem(A,self.speeds,1)
        
            #calculate the radiance map
            #weighted average of intensities of each image
            rows=self.images[0].shape[0]
            cols=self.images[0].shape[1]
            radianceMap=np.zeros((rows,cols),np.uint32)
            normalization=np.zeros((rows,cols),np.uint32)
            vecWeight=np.vectorize(self.weight)
            for imageNum in range(len(self.images)):
                img=cv2.split(self.images[imageNum])[channel] #get the correct channel
                print imageNum,"/",len(self.images)
                normalization=normalization+vecWeight(img)
                temp=vecWeight(img)
                speed=self.speeds[imageNum]
                for r in range(rows):
                    for c in range(cols):
                        temp[r,c]*=(g[img[r,c]]-speed)
                radianceMap=radianceMap+temp
            np.divide(radianceMap,normalization) 
            channels.append(radianceMap)
        #combine channels
        combined=cv2.merge([channels[0],channels[1],channels[2]])

        #set hdr image
        self.hdr=combined

        # Hardcode in file saving for now
        #
        #save file with opencv
        #HDR saving isn't working for me, but JPG,TIFF,EXR do work
        #cv2.imwrite('output.hdr',combined)
        
        image=combined
        #HDR saving code found at https://gist.github.com/edouardp/3089602
        f = open("output.hdr", "wb")
        f.write("#?RADIANCE\n# Made with Python & Numpy\nFORMAT=32-bit_rle_rgbe\n\n")
        f.write("-Y {0} +X {1}\n".format(image.shape[0], image.shape[1]))
 
        brightest = np.maximum(np.maximum(image[...,0], image[...,1]), image[...,2])
        mantissa = np.zeros_like(brightest)
        exponent = np.zeros_like(brightest)
        scaled_mantissa = np.zeros_like(brightest)
        np.frexp(brightest, mantissa, exponent)
        for i in range(scaled_mantissa.shape[0]):
            for j in range(scaled_mantissa.shape[1]):
                if brightest[i,j]==0:
                    scaled_mantissa[i,j]=0
                else:
                    scaled_mantissa[i,j] = mantissa[i,j] * 256.0 / brightest[i,j]
        rgbe = np.zeros((image.shape[0], image.shape[1], 4), dtype=np.uint8)
        rgbe[...,0:3] = np.around(image[...,0:3] * scaled_mantissa[...,None])
        rgbe[...,3] = np.around(exponent + 128)
  
        rgbe.flatten().tofile(f)
        f.close()



    def weight(self,w):
        min=1
        max=255
        if w<= (min+max)/2:
            weighted= w-min
        else:
            weighted=max-w
        return weighted
    def solveSystem(self,Z,B,lamb):
        n=256
        rowsZ=Z.shape[0]
        colsZ=Z.shape[1]
        rowsA=rowsZ*colsZ+n+1
        colsA=n+rowsZ
        A=np.zeros((rowsA,colsA),np.uint32)
        b=np.zeros((rowsA,1),np.uint32)
        
        k=1
        for i in range(0,rowsZ):
            for j in range(0,colsZ):
                currentWeight=self.weight(Z[i,j]+1)
                A[k,Z[i,j]+1] = currentWeight
                A[k,n+i]=-currentWeight
                b[k,0]=currentWeight*B[j]
        #Set middle val = 1
        A[k,129]=1
        k+=1

        #Smoothness
        for i in range(0,n-2):
            wplus=self.weight(i+1)
            A[k,i]=lamb*wplus
            A[k,i+1]=-2*lamb*wplus
            A[k,i+2]=lamb*wplus
        sol=np.linalg.lstsq(A,b)[0]  #least squares for overdetermined system
        g=sol[1:n]
        lE=sol[n+1:]
        print len(g)
        print len(lE)
        return sol

def main():
    if len(sys.argv)<2:
        sys.exit("usage: HDRCreator imagelist")
    hd=HDRCreator()
    hd.loadimages(sys.argv[1])
    hd.createHDR()
    


if __name__ == "__main__":
    main()
