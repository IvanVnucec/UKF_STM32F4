f = open("Acc_2021-02-21_231529.txt", "r")

acc = []
for line in f.readlines()[1:]:
    aa = line.split()[1:]       
    bb = line.split()[1:]
    for i, a in enumerate(bb):
        bb[i] = f"{a}f"

    acc.append(bb)

acc = acc[0::8]

dataTowrite = []
dataTowrite.append(f"const float accData[{len(acc)}][3] = {{\n")
for a in acc[:-1]:
    dataTowrite.append(f"\t{{{a[0]}, {a[1]}, {a[2]}}},\n")

dataTowrite.append(f"\t{{{acc[-1][0]}, {acc[-1][1]}, {acc[-1][2]}}}\n")
dataTowrite.append("};\n")

f.close()
f = open("IMUdata.c", 'w')
f.writelines(dataTowrite)
f.close()


f = open("Gyr_2021-02-21_231529.txt", "r")
acc = []
for line in f.readlines()[1:]:
    aa = line.split()[1:]       
    bb = line.split()[1:]
    for i, a in enumerate(bb):
        bb[i] = f"{a}f"

    acc.append(bb)

acc = acc[0::8]

dataTowrite = []
dataTowrite.append(f"const float gyroData[{len(acc)}][3] = {{\n")
for a in acc[:-1]:
    dataTowrite.append(f"\t{{{a[0]}, {a[1]}, {a[2]}}},\n")

dataTowrite.append(f"\t{{{acc[-1][0]}, {acc[-1][1]}, {acc[-1][2]}}}\n")
dataTowrite.append("};\n")

f.close()
f = open("IMUdata.c", 'a')
f.writelines(dataTowrite)
f.close()



f = open("Mag_2021-02-21_231529.txt", "r")

acc = []
for line in f.readlines()[1:]:
    aa = line.split()[1:]       
    bb = line.split()[1:]
    for i, a in enumerate(bb):
        bb[i] = f"{a}f"

    acc.append(bb)

acc = acc[:len(dataTowrite)-2]
dataTowrite = []
dataTowrite.append(f"\nconst float magData[{len(acc)}][3] = {{\n")
for a in acc[:-1]:
    dataTowrite.append(f"\t{{{a[0]}, {a[1]}, {a[2]}}},\n")

dataTowrite.append(f"\t{{{acc[-1][0]}, {acc[-1][1]}, {acc[-1][2]}}}\n")
dataTowrite.append("};\n")

f.close()
f = open("IMUdata.c", 'a')
f.writelines(dataTowrite)
f.close()