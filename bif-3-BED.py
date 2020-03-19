import sys

for line in sys.stdin:
    parts = line.split("\t")
    length = int(parts[2]) - int(parts[1])
    if length > 100000:
        print(line, end="")
