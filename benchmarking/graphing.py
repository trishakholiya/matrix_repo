import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("bench.csv", skiprows=4)

# Extract group name, library, and N
df["group"] = df["name"].apply(lambda x: x.split("_")[0])
df["library"] = df["name"].apply(lambda x: x.split("_")[1].split("/")[0])
df["N"] = df["name"].apply(lambda x: int(x.split("/")[-1]))

# Functions We Benchmarked
groups_to_plot = ["Addition", "Accessor", "Multiplication", "Transpose", "InitializationRandom", "InitializationZeros", "EigSym"]

for g in groups_to_plot:
    sub_group = df[df["group"] == g]

    mc = sub_group[sub_group["library"] == "MatrixClass"]
    arma = sub_group[sub_group["library"] == "Armadillo"]

    plt.figure(figsize=(7,5))

    plt.plot(mc["N"], mc["real_time"], marker="o", label="MatrixClass", color = 'red')
    plt.plot(arma["N"], arma["real_time"], marker="o", label="Armadillo", color = 'blue')

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("log Matrix size N")
    plt.ylabel("log Runtime (ns)")
    plt.title(f"{g} Benchmark")

    plt.legend()

    plt.savefig(f"{g}_benchmark.png")