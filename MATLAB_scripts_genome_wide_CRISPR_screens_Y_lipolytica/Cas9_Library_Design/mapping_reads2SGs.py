import sys
import csv
input_file = sys.argv[1]
output_file = sys.argv[2]

data = []
with open("SGs.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    next(reader, None)
    for row in reader:
        s = row[0].strip()
        cols = s.split(',')
        data.append(cols)


hmap_count = {}


for row in data:
    seq = row[2]
    if seq in hmap_count:
        hmap_count[seq] += 1
    else:
        hmap_count[seq] = 1
    

hmap = {}
for row in data:
    seq = row[2]
    if hmap_count[seq] == 1:
        hmap[seq] = 0

n_line = 0
for line in open(input_file):
    line = line.strip()
    n_line += 1
    if n_line % 4 == 1 or n_line % 4 == 3 or n_line % 4 == 0:
        continue
    read = set([])
    for i in range(0, len(line)-20+1):
        seq = line[i:i+20]
        read.add(seq)
    for seq in read:
        if seq in hmap:
            hmap[seq] += 1
    if n_line % 100000 == 0:
        print n_line

fo = file(output_file, 'w')
fo.write("Number,Gene,Sequence,Count\n")
for row in data:
    sg_id = row[0]
    target = row[1]
    seq = row[2]
    if seq in hmap:
        count = hmap[seq]
    else:
        count = "N/A"
    fo.write(sg_id + "," + target + "," + seq+"," + str(count) + "\n")






