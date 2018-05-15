import datetime

def td_total_seconds(td):
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6)

prev = datetime.datetime.now()
flarge="/home/hmayani/ParallelProgramming/Project/mapR/word_count_datafiles/enwiki-latest-abstract.xml"
fmid="/home/hmayani/ParallelProgramming/Project/mapR/word_count_datafiles/word_50MB.txt"
text_file = sc.textFile(fmid)
counts = text_file.flatMap(lambda line: line.split(" ")).map(lambda word: (word, 1)).reduceByKey(lambda a, b: a + b)
print(len(counts.collect()))
now =datetime.datetime.now()
td_total_seconds(now -prev)
