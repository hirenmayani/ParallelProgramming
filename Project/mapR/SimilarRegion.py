
import pyspark as ps
from scipy import linalg
from pprint import pprint
import io
import zipfile
import numpy as np
import hashlib
from tifffile import TiffFile
from pyspark.mllib.linalg import Vectors
from pyspark.mllib.linalg.distributed import RowMatrix

def getOrthoTif(zfBytes):
 #given a zipfile as bytes (i.e. from reading from a binary file),
# return a np array of rgbx values for each pixel
 bytesio = io.BytesIO(zfBytes)
 zfiles = zipfile.ZipFile(bytesio, "r")
 #find tif:
 for fn in zfiles.namelist():
  if fn[-4:] == '.tif':#found it, turn into array:
   tif = TiffFile(io.BytesIO(zfiles.open(fn).read()))
   return tif.asarray()

def splitter(x):
    step = x[1].shape[0]/500

    harr = np.vsplit(x[1], step)
    varr = []
    for arr in harr:
        varr+=np.hsplit(arr, step)

    return (x[0],varr)


def mapper(x):
    pixels = []
    for i in range(len(x[1])):
        pixels.append((x[0]+"-"+str(i),x[1][i]))

    return pixels

def filterPixels(pixel):
    return (pixel[0] == '3677454_2025195.zip-0' or pixel[0] == '3677454_2025195.zip-1' or pixel[0] == '3677454_2025195.zip-18' or pixel[0] == '3677454_2025195.zip-19')

def calculateIntensity(x):
    subImage = x[1]*1.0
    rbg_values = subImage[:,:,:3]
    rgb_values_mean = rbg_values.mean(2)
    rgb_values_mean = rgb_values_mean[:,:,None]
    infrared = subImage[:,:,3:]
    infrared = infrared/100
    #np.concatenate((rgb_values_mean, infrared))

    intensities = rgb_values_mean*infrared
    intensities = intensities.reshape((500,500))

    return (x[0],intensities)

def reduceResolutionSplit(x):
    step = 50

    harr = np.vsplit(x[1], step)

    varr = []
    for arr in harr:
        varr += np.hsplit(arr, step)
    #print(len(varr))
    varr = np.array(varr)
    varr = varr.mean((1,2))
    varr=varr.reshape((50,50))

    return (x[0], varr)

def reduceResolutionSplitBonus(x):
    step = 100

    harr = np.vsplit(x[1], step)

    varr = []
    for arr in harr:
        varr += np.hsplit(arr, step)
    #print(len(varr))
    varr = np.array(varr)
    varr = varr.mean((1,2))
    #print(varr.shape)
    varr=varr.reshape((100,100))

    return (x[0], varr)

def quantify(x):
    npArr = x[1]
    npArr[np.logical_and(npArr >= -1, npArr <= 1)] = 0
    npArr[npArr < -1] = -1
    npArr[npArr > 1] = 1
    return(x[0],npArr)

def orderedConcat(x,y):
    if x[0] == "r":
        return np.concatenate((x[1],y[1]))
    else:
        return np.concatenate((y[1], x[1]))


def generateSignature(x):

    size = int(len(x[1])/128)*128
    resizedArray = x[1][:size]
    chunkSize = int(len(x[1])/128)


    bitString = []
    #
    # print(x[1])
    for i in range(128):
        #print("chunk:"+str(i*chunkSize)+"-"+str((i+1)*chunkSize))
        chunk = resizedArray[i*chunkSize : (i+1)*chunkSize]
        digest = hashlib.md5(chunk).hexdigest()
        bitString.append(format(ord(digest[0]),"b")[:6])
    #print(len(bitString))
    return(x[0],bitString)


def divideIntoBands(x):

    rowsPerBand = 4
    bands = 128/rowsPerBand
    #banded = np.zeros((int(bands),rowsPerBand))
    #banded = np.empty((int(bands), rowsPerBand))
    #print(banded)
    banded = []
    for i in range(int(bands)):
        banded.append(x[1][i*rowsPerBand : (i+1)*rowsPerBand])
    print(banded)
    return (x[0],np.array(banded))
def hashBands(x):

    zip,bands = x[0],x[1]
    hashes = []
    for i in range(len(bands)):
        hashes.append((hashlib.md5(bands[i]).hexdigest()+"band"+str(i),zip))
    return hashes


def addToList(x,y):
    if not isinstance(x, list):
        x = [x]
    if not isinstance(y, list):
        y = [y]

    return x + y

def classify(x):
    candidates = x[1]
    #print(candidates)
    if not isinstance(candidates,list):
        candidates = [candidates]
    find = ['3677454_2025195.zip-0', '3677454_2025195.zip-1','3677454_2025195.zip-18','3677454_2025195.zip-19']
    match = []
    for zip in find:
        if zip in candidates:
            for candidate in candidates:
                if zip != candidate:
                    match.append((zip + "--" + candidate, 1))

    return match


def countTop20Candidates(zippedCandidates):
    count = []
    for zip,candidates in zippedCandidates.items():
        for candidate in candidates:
            #print(candidate)
            if zip!=candidate:
                count.append((zip+"--"+candidate,1))
    return count

def arrange(obj):
    pair,count = obj

    if count == 0:
        return

    zips = pair.split("--")
    zip = zips[0]
    candy = zips[1]
    return (zip,{candy:count})

def addToDict(x,y):
    if x is None and y is None:
        return
    if x is None:
        return y
    if y is None:
        return x
    for k,v in x.items():
        try:
            y[k] += v
        except:
            y[k] = v
    return y

def PCA(x):
    zip,img_diffs = x
    mu, std = np.mean(img_diffs, axis=0), np.std(img_diffs, axis=0)
    img_diffs_zs = (img_diffs - mu) / std
    #print(img_diffs_zs.shape)
    # run singular value decomposition.
    # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.linalg.svd.html
    # Read more about svd options from above link
    U, s, Vh = linalg.svd(img_diffs_zs, full_matrices=1)
    # Now reduce the number of dimensions by simply taking the first low_dim_p dimensions of columns (aka singular values).
    low_dim_p = 10
    img_diffs_zs_lowdim = U[:, 0:low_dim_p]
    return zip,img_diffs_zs_lowdim


#np.set_printoptions(threshold=np.nan)
sc = ps.SparkContext()


path = "/Users/krishnasharma/Downloads/a2_small_sample"
path1 = "hdfs:/data/small_sample"
path2 = "hdfs:/data/large_sample"

#Load files from server

tfrdd = sc.binaryFiles(path2)#="a2_small_sample")
tfrdd = tfrdd.map(lambda x:(x[0].split("/")[-1],x[1]))
tfrdd = tfrdd.map(lambda x:(x[0],getOrthoTif(x[1])))

#Split row-wise and column-wise to get 500x500 submimages
tfrdd = tfrdd.map(splitter)

#Obtain individual zips of the form - ('3677502_2035200.zip-5',np array (500x500))
tfrdd = tfrdd.flatMap(mapper)

#Print rgb values for particular files
pixels = tfrdd.filter(filterPixels)
pixels = pixels.map(lambda x:(x[0],np.array(x[1][0][0])))
pixels = pixels.collect()
for pixel in pixels:
    print(pixel)

#Calculate intensity values for each pixel
tfrdd = tfrdd.map(calculateIntensity)


#Reducing Resolution by taking mean over 10x10 subimages
tfrddOriginal = tfrdd.map(reduceResolutionSplit)

#Taking difference to reduce effect of shadow row-wise
row_diff = tfrddOriginal.map(lambda x:(x[0],np.diff(x[1])))

#Quantifying difference to -1,0,1
row_diff = row_diff.map(quantify)

#Taking difference to reduce effect of shadow col-wise
col_diff = tfrddOriginal.map(lambda x:(x[0],np.diff(x[1],axis=0)))

#Quantifying difference to -1,0,1
col_diff = col_diff.map(quantify)

#Mapping rows and cols to ensure then get concatenated in reduce in same order for all features
row_diff_rdd = row_diff.map(lambda x:(x[0],("r",x[1].flatten())))
col_diff_rdd = col_diff.map(lambda x:(x[0],("c",x[1].flatten())))

#Merge 2 feature vectors
merged = row_diff_rdd.union(col_diff_rdd)


#Concatinating row-wise and col-wise features
features = merged.reduceByKey(orderedConcat).persist()

#Printing Features for particular files
print()
printFeatures = features.filter(lambda pixel:pixel[0] == '3677454_2025195.zip-1' or pixel[0] == '3677454_2025195.zip-18')
printFeatures = printFeatures.collect()
#print(printFeatures)
for feature in printFeatures:
    print(feature)
print()


#Create 128bit signature
signatures = features.map(generateSignature)

#Divide signatures into bands for LSH
divided = signatures.map(divideIntoBands)

#Hash Each Band of an image to their respective bucket
hashedBands = divided.flatMap(hashBands)

#Find all images which hash to the same bucket in any of the bands
similar = hashedBands.reduceByKey(addToList)

#classify candidates in 4 zips -
similarMap = similar.flatMap(classify)

##Adding counts of 2 images matched in multiple bands
similarMap = similarMap.reduceByKey(lambda x,y: x+y)

#Rearranging in the form (zip:(candidate,count))
similarMap = similarMap.map(arrange)

#Adding multiple dictionaries
similarMap = similarMap.reduceByKey(addToDict)

#
mapped = similarMap.collect()
import operator
sorted_dicts = {}
for zip,candidates in mapped:
    sorted_dicts[zip] = sorted(candidates.items(), key=operator.itemgetter(1),reverse=True)
for zip,candidates in sorted_dicts.items():
    #print("-----------------------------------")
    #print(len(candidates))
    processC = candidates[:20]
    finalC = []
    for c in processC:
        finalC.append(c[0])
    sorted_dicts[zip] = finalC
    if zip == '3677454_2025195.zip-1' or zip == '3677454_2025195.zip-18':
        print(zip)
        pprint(sorted_dicts[zip])
        print()

#pprint(sorted_dicts)

allFeatures = features.collect()
featureDict = {}
for zip,feature in allFeatures:
    featureDict[zip] = feature

lowD = []

for zip,topCandidates in sorted_dicts.items():
    featureVector = []
    featureVector.append(featureDict[zip])
    for candidate in topCandidates:
        featureVector.append(featureDict[candidate])
    lowD.append((zip,np.array(featureVector)))

svd = sc.parallelize(lowD)
svd = svd.map(PCA)
zippedVectors = svd.collect()
#print(zippedVectors)
distDict = {}
for zip,vectors in zippedVectors:
    reference = vectors[0]
    dists = []
    for vector in vectors[1:]:
        dists.append(np.linalg.norm(reference - vector))
    distDict[zip] = dists

distDictSorted = {}
for zip,distances in distDict.items():
    topCandidates = sorted_dicts[zip]
    candy = {}
    #print(len(distances))
    #print(len(topCandidates))
    for i in range(len(topCandidates)):
        candy[topCandidates[i]] = distances[i]
    distDictSorted[zip] = sorted(candy.items(), key=operator.itemgetter(1))
    if zip == '3677454_2025195.zip-1' or zip == '3677454_2025195.zip-18':
        print(zip)
        pprint(distDictSorted[zip])
        print()


#pprint(distDictSorted)



#pprint(distDict)


###################################################################################################################################################
###################################################################################################################################################
################################################      BONUS BEGINS FROM HERE     #################################################################
###################################################################################################################################################
###################################################################################################################################################

print("BONUS BEIGNS HERE")
#Reducing Resolution by taking mean over 5x5 subimages
tfrddOriginal = tfrdd.map(reduceResolutionSplitBonus)

#Taking difference to reduce effect of shadow row-wise
row_diff = tfrddOriginal.map(lambda x:(x[0],np.diff(x[1])))

#Quantifying difference to -1,0,1
row_diff = row_diff.map(quantify)

#Taking difference to reduce effect of shadow col-wise
col_diff = tfrddOriginal.map(lambda x:(x[0],np.diff(x[1],axis=0)))

#Quantifying difference to -1,0,1
col_diff = col_diff.map(quantify)

#Mapping rows and cols to ensure then get concatenated in reduce in same order for all features
row_diff_rdd = row_diff.map(lambda x:(x[0],("r",x[1].flatten())))
col_diff_rdd = col_diff.map(lambda x:(x[0],("c",x[1].flatten())))

#Merge 2 feature vectors
merged = row_diff_rdd.union(col_diff_rdd)


#Concatinating row-wise and col-wise features
features = merged.reduceByKey(orderedConcat).persist()

#Create 128bit signature
signatures = features.map(generateSignature)

#Divide signatures into bands for LSH
divided = signatures.map(divideIntoBands)

#Hash Each Band of an image to their respective bucket
hashedBands = divided.flatMap(hashBands)

#Find all images which hash to the same bucket in any of the bands
similar = hashedBands.reduceByKey(addToList)

#classify candidates in 4 zips -
similarMap = similar.flatMap(classify)

similarMap = similarMap.reduceByKey(lambda x,y: x+y)
similarMap = similarMap.map(arrange)
similarMap = similarMap.reduceByKey(addToDict)
mapped = similarMap.collect()
import operator
sorted_dicts = {}
for zip,candidates in mapped:
    sorted_dicts[zip] = sorted(candidates.items(), key=operator.itemgetter(1),reverse=True)

for zip,candidates in sorted_dicts.items():
    #print("-----------------------------------")
    #print(len(candidates))
    processC = candidates[:20]
    finalC = []
    for c in processC:
        finalC.append(c[0])
    sorted_dicts[zip] = finalC

#pprint(sorted_dicts)
#pprint(sorted_dicts)

allFeatures = features.collect()
featureDict = {}
for zip,feature in allFeatures:
    featureDict[zip] = feature
lowD = []
for zip,topCandidates in sorted_dicts.items():
    featureVector = []
    featureVector.append(featureDict[zip])
    for candidate in topCandidates:
        featureVector.append(featureDict[candidate])
    lowD.append((zip,np.array(featureVector)))

distDict = {}
for i in range(len(lowD)):
    zip = lowD[i][0]
    rows = sc.parallelize(lowD[i][1])
    mat = RowMatrix(rows)

    svd = mat.computeSVD(10, computeU=True)
    U = svd.U  # The U factor is a RowMatrix.
    distances = U.rows.collect()
    reference = distances[0]
    dists = []
    for vector in distances[1:]:
        dists.append(np.linalg.norm(reference - vector))
    distDict[zip] = dists
distDictSorted = {}
for zip,distances in distDict.items():
    topCandidates = sorted_dicts[zip]
    candy = {}
    #print(len(distances))
    #print(len(topCandidates))
    for i in range(len(topCandidates)):
        candy[topCandidates[i]] = distances[i]
    distDictSorted[zip] = sorted(candy.items(), key=operator.itemgetter(1))
    if zip == '3677454_2025195.zip-1' or zip == '3677454_2025195.zip-18':
        print(zip)
        pprint(distDictSorted[zip])
        print()

#pprint(distDictSorted)



