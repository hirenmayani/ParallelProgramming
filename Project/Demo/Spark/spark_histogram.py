#!/usr/bin/python

from PIL import Image
import pyspark
import numpy as np
import time

fpath = "BLU.BMP"
start_time = time.time()

def load_image(infilename) :
    img = Image.open( infilename )
    img.load()
    data = np.asarray( img, dtype="int32" )
    return data

def save_image( npdata, outfilename ) :
    img = Image.fromarray( np.asarray( np.clip(npdata,0,255), dtype="uint8"), "L" )
    img.save( outfilename )

img = load_image(fpath)

output = [0 for i in range(256*3)]

sc = pyspark.SparkContext.getOrCreate()
rddData = sc.parallelize(img.flatten()).map(lambda x: (x, 1)).reduceByKey(lambda x, y: x + y).collect()

print("--- %s seconds ---" % (time.time() - start_time))


import sys
import math
import numpy



# to covert ipixel to bit
def pixel_to_bit(pixel,threshold=128):
    """Convert the pixel value to a bit."""
    if pixel > threshold:
        return 1
    else:
        return 0

# to covert image to pixel array
def convert(im,center,inner_radius,outer_radius,fill=255,mark=False):
    """
    Convert the pix data into an integer array of 256 values.
    Each integer value is the LED pattern along a radius.
    
    Arguments:
    im :  An Image; currently the pixels in the imag must be a single
          value (and not, say, an RGB triple)
    center : An (x,y) tuple of the center of the sampling region
    inner_radius : The inner radius of the sampling region
    outer_radius : The outer radius of the sampling region
    fill : The value to use for samples that are outside the bounds of the image.
    mark : If True, the input image im is modified by inverting the value of
           each sampled point in the sampling region.
           
    Return value:
    A list of 256 lists of 32 binary values.  The 32 binary values are the
    samples for each "spoke".
    """
    pix = im.load()
    result = []
    for theta in numpy.linspace(0.0,2*numpy.pi,num=256,endpoint=False):
        spokedata = []
        for r in numpy.linspace(inner_radius,outer_radius,num=32):
            i = int(center[0] + r*math.cos(theta))
            j = int(center[1] + r*math.sin(theta))
            if i >= im.size[0] or i < 0 or j >= im.size[1] or j < 0:
                pixel = fill
            else:
                pixel = pix[i,j]
                if not mark is False:
                    pix[i,j] = 255-pix[i,j]
            bit = pixel_to_bit(pixel)
            spokedata.append(bit)
        result.append(spokedata)
    return result



       

#
# Main script begins here.
#
# First check for keyword arguments on the command line.
#
def client():
    kwargs = {'cx':None, 'cy':None, 'ri':None, 'ro':None, 'fill':None, 'mark':None }

    if len(sys.argv) > 2:
        for kwa in sys.argv[2:]:
            if not '=' in kwa:
                exit(-1)
            (w,v) = kwa.split('=')
            if w in kwargs:
                kwargs[w]=int(v)
            else:
                exit(-1)

    #
    # Open the image file and get the basic information about it.
    #
    try:
        im = Image.open(sys.argv[1])
    except:
        # Eventually this should give more useful information (e.g. file does not
        # exist, or not an image file, or ...
        print( "Unable to open %s" % sys.argv[1])
        exit(-1)

    width,height = im.size
    if im.mode == "RGB":
        exit(0)

    #
    # Set the default options for any that were not given on the command line.
    # (Do this after opening the file, because some of the default options depend
    # on the image data.)
    #
    if kwargs['cx'] is None:
        kwargs['cx'] = width/2
    if kwargs['cy'] is None:
        kwargs['cy'] = height/2
    if kwargs['ri'] is None:
        kwargs['ri'] = 2
    if kwargs['ro'] is None:
        kwargs['ro'] = width/2
    if kwargs['fill'] is None:
        pix = im.load()
        kwargs['fill'] = pix[0,0]
    if kwargs['mark'] is None:
        kwargs['mark'] = 0

    #
    # Do it!
    #
    samples = convert(im,(kwargs['cx'],kwargs['cy']),kwargs['ri'],kwargs['ro'],
                    fill=kwargs['fill'],mark=kwargs['mark'])

    #
    # Write the data to a new file.
    # This part needs to be rewritten to put the output file into the desired form.
    # For now, just write something readable to "radial_samples.dat"
    #
    f = open("radial_samples.dat",'w')
    for d in samples:
        for b in d:
            f.write("%d " % b)
        f.write('\n')
    f.close()

    #
    # If mark was set, the sampled pixels have been inverted.
    # Save the image to a new file, so we can look at the sample pattern.
    #
    if kwargs['mark'] != 0:
        im.save("tmp.png")