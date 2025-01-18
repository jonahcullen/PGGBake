import fileinput

# from cactus/doc/mc-paper/hprc/resources/compute_mapping_stats.py

# parse the SAM record and count how many reads
# in each profile (mapping quality x perfect alignement x aligment score)
#    mapping quality: number. '-1' for unmapped reads
#    perfect: boolean to specify if the reads aligned perfectly
#    alignment score: value of the AS tag

# to tally the number of reads/records in each profile
records = {}

for line in fileinput.input():
    if line[0] == "@":
        continue
    line = line.rstrip().split('\t')
    ## skip if secondary or supplementary alignment
    if int(int(line[1])/256)%2 == 1 or int(int(line[1])/2048)%2 == 1:
        continue
       #print('HERE ',line[0])
    ## handle unmapped reads
    if line[5] == "*":
        line[4] = '-1'
    ## extract MD and AS tags
    as_tag = False
    md_tag = False
    for ii in range(11, len(line)):
        tag = line[ii].split(':')
        if tag[0] == 'MD':
            md_tag = line[ii]
        elif tag[0] == 'AS':
            as_tag = tag[2]
        if as_tag and md_tag:
            break
    ## check if perfectly aligned
    perfect = False
    if line[5] == str(len(line[9])) + "M":
        if md_tag and md_tag == "MD:Z:" + str(len(line[9])):
            perfect = True
    ## increment records for this profile/rid
   #print('okay ', line[4])
    rid = '{}\t{}\t{}'.format(line[4], perfect, as_tag)
    if rid not in records:
        records[rid] = 1
    else:
        records[rid] += 1

# print the counts for each rid/profiles
for rec in records.keys():
    print('{}\t{}'.format(records[rec], rec))
