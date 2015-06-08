import sys

def calc_results(rsem_out, t_out):
    with open(rsem_out, 'r') as rsem_in:
        rsem_data = {}
        for (i, elem) in enumerate(rsem_in.readlines()[1:]):
            elem = elem.split("\t")
            rsem_data[elem[0]] = float(elem[-3])

    with open(t_out, 'r') as data_in:
        test_data = {}
        for (i, elem) in enumerate(data_in.readlines()):
            if elem[0] == ">":
                name = elem.split(" ")[0][1:]
            else:
                test_data[name] = float(elem.split("\t")[-1])

    with open("tpm-test.out", "w") as out:
        for k in rsem_data.keys():
            out.write(str(abs(rsem_data[k] - test_data[k])))
            out.write("\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("./test_tpm rsem-out your-out")
        sys.exit(1)

    calc_results(sys.argv[1], sys.argv[2])
