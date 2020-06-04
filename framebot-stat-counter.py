# count of framebot output based on identity%

stat_filename = input('enter framebot-stat file: ')
id_perc = []

with open(stat_filename,'r') as infile:
    next(infile) # skip header
    for line in infile:
        y = line.split('\t')
        id_perc.append(eval(y[4])) #append ID% column, eval to convert str to number (float)

#print(id_perc)


# count based on percent identity

count30 = 0
count40 = 0
count50 = 0
count60 = 0
count70 = 0
count80 = 0
count90 = 0
count100 = 0

for i in id_perc:
    if i > 30 and i < 40:
        count30 +=1
    if i > 40 and i < 50:
        count40 +=1
    if i > 50 and i < 60:
        count50 +=1
    if i >60 and i < 70:
        count60 +=1
    if i >70 and i < 80:
        count70 +=1
    if i >80 and i < 90:
        count80 +=1
    if i >90 and i < 100:
        count90 +=1
    if i == 100:
        count100 += 1

print(">30%", count30)
print(">40%",count40)
print(">50%",count50)
print(">60%",count60)
print(">70%",count70)
print(">80%",count80)
print(">90%",count90)
print("100%",count100)



