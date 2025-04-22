import subprocess

#input_file = "HG003_BBBC20_HGSVC3.vcf.gz"
input_file = "HG005_BBBC5_HGSVC3.vcf.gz"
subprocess.run(["gunzip", input_file])

fi = open(input_file[:-3], "r")
fo = open("tmp.vcf", "w")

flag_skip = False
for line in fi:
    fields = line.split("\t")
    if fields[0] == "chr22":
        if fields[1] == "12977506":
            flag_skip = True
        elif fields[1] == "15154394":
            flag_skip = False
        if flag_skip:
            continue
        else:
            fo.write(line)
    else:
        fo.write(line)
fo.close()
fi.close()

subprocess.run(["bgzip", "tmp.vcf"])
subprocess.run(["mv", "tmp.vcf.gz", input_file])
subprocess.run(["tabix", "-p", "vcf", input_file])
