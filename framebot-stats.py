# program to extract statistics from the framebot file

framebot_filename = input('enter non_chimeric_framebot.txt file: ')
framebot = []
#non_chimeric_framebot_1-13.txt
with open(framebot_filename,'r') as infile:
    for line in infile:
        y = line.split('\t')
        framebot.append(y)

#pull out stats line from each entry
seq_stats = []
for i in framebot:
    if "STATS" in i:
        seq_stats.append(i[1:8])

#print(seq_stats)

headers = ['Target','Query','NuclLen','AlignLen','%Identity','Score','Frameshifts','Reversed']

dat = [headers] + seq_stats
'''
output_data = input("enter output file name: ")
with open(output_data, "w") as out_file:
    out_file.writelines('\t'.join(i) + '\n' for i in dat)
'''

