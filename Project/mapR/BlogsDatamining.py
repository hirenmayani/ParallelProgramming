import pyspark
import re
from bs4 import BeautifulSoup
from pprint import pprint
from os import listdir
import string


def tagFinder(content):
    print(content)
    soup = BeautifulSoup(content[1], "lxml")
    date_post = soup.find_all(["date", "post"])
    dp_tups = []
    j = 0
    for i in range(int(len(date_post) / 2)):
        date = date_post[j].get_text(strip=True)[3:]
        post = date_post[j + 1].get_text(strip=True)
        dp_tups.append((date, post))
        j += 2
    return dp_tups

def findIndutries(dp_tup):
    if not dp_tup == "":
        date = "-".join(dp_tup[0].split(","))
    words = dp_tup[1].split(" ")
    indL = []
    translator = str.maketrans("","",string.punctuation)

    for word in words:
        x = word.lower().translate(translator)
        if x in broadcastVar.value:
            indL.append((date,x))
    return indL

def addLists(x,y):
    if not isinstance(x, list):
        x = [x]
    if not isinstance(y, list):
        y = [y]

    return x + y

""" Begin here """

sc = pyspark.SparkContext.getOrCreate()
path = "/Users/admin/Downloads/blogs"


""" Get all possible industry names: """
industryNames = []
for file in listdir(path):
    industryNames.append(file)


""" Create an rdd of all the filenames  """
industryNamesRDD = sc.parallelize(industryNames)


""" Use transformations until you are left with only a set of possible industries   """
industryNamesRDD = industryNamesRDD.map(lambda x: x.split('.')[3].lower())
industryNamesRDD = industryNamesRDD.distinct()
industryNamesRDD.unpersist()


""" Use an action to export the rdd to a set and make this a spark broadcast variable   """
industryNames = industryNamesRDD.collect()
#print(industryNames)
broadcastVar = sc.broadcast(industryNames)


""" Create an rdd for the contents of all files [i.e. sc.wholeTextFiles(file1,file2,...) ]  """
dataResRDD = sc.wholeTextFiles(path)


"""
Use transformations to search all posts across all blogs for mentions of industries, and record the frequency 
each industry was mentioned by month and year.
"""
RDD = dataResRDD.flatMap(tagFinder).flatMap(findIndutries).map(lambda x:((x[1],x[0]),1)).reduceByKey(lambda x,y: x+y)\
    .map(lambda x:(x[0][0],(x[0][1],x[1]))).reduceByKey(addLists).sortByKey(True)


""" Use an action to print the recorded frequencies in this format: """
pprint(RDD.collect())





