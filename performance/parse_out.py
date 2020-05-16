#BASE_FREQ = 1.6 * 10**9

# parse timing output
c = []
print("\ncycles:")
with open("cycles.txt", "r") as f:
    for line in f:
        if line.startswith("RDTSC"):
            cycles = line.split(" ")[2]
            #time = int(cycles) / BASE_FREQ
            print(cycles)
            c.append(int(cycles))
print()

# parse inctrumentation output
o = []
print("flops:")
with open("counts.txt", "r") as f:
    for line in f:
        if line.startswith("ADDS"):
            adds = line.split(" ")[0].split("=")[1]
            mults = line.split(" ")[1].split("=")[1]
            ops = int(adds) + int(mults)
            print(ops)
            o.append(ops)

print()

# calculate performance [flop/cycle

perf = []
print("performances [flops/cycle]")
for cnt, op in zip(c, o):
    perf = op / cnt
    print(perf)
            
