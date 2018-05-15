import pyspark
from random import random


# WordCount implementation below
def wordCountImplementation(data):
    sc = pyspark.SparkContext.getOrCreate()
    RDDRead = sc.parallelize(data)
    reducedRDD = RDDRead.flatMap(lambda x: x[1].lower().split()).map(lambda x: (x, 1)).reduceByKey(lambda x, y: x + y)
    print(reducedRDD.collect())


# Set Difference Implementation below
def setDifferenceImplementation(data):
    sc = pyspark.SparkContext.getOrCreate()
    RDDRead = sc.parallelize(data)
    m = RDDRead.flatMap(lambda x: [(i, x[0]) for i in x[1]]).reduceByKey(lambda x, y: x + y).map(
        lambda x: (x[1], x[0])).reduceByKey(func).map(lambda x: x if x[0] == 'R' else None)
    print(m.collect()[1])


def func(x, y):
    if x == None or y == None:
        return None

    if not isinstance(x, list):
        x = [x]
    if not isinstance(y, list):
        y = [y]

    return x + y


data = [(1, "The horse raced past the barn fell"),
        (2, "The complex houses married and single soldiers and their families"),
        (3, "There is nothing either good or bad, but thinking makes it so"),
        (4, "I burn, I pine, I perish"),
        (5, "Come what come may, time and the hour runs through the roughest day"),
        (6, "Be a yardstick of quality."),
        (7, "A horse is the projection of peoples' dreams about themselves - strong, powerful, beautiful"),
        (8,
         "I believe that at the end of the century the use of words and general educated opinion will have altered so much that one will be able to speak of machines thinking without expecting to be contradicted."),
        (9, "The car raced past the finish line just in time."),
        (10, "Car engines purred and the tires burned.")]
data1 = [('R', ['apple', 'orange', 'pear', 'blueberry']),
         ('S', ['pear', 'orange', 'strawberry', 'fig', 'tangerine'])]
data2 = [('R', [x for x in range(50) if random() > 0.5]),
         ('S', [x for x in range(50) if random() > 0.75])]

""" Put input here """
wordCountImplementation(data)
setDifferenceImplementation(data1)
setDifferenceImplementation(data2)
