import tadasets
import petls
folder = f"./torus_L"

def get_complex(replicate):
    replicate = 0
    # n=500 is a good number for the actual test, but not for debugging
    torus = tadasets.torus(n=500, c=3, a=1, noise=0.0, seed=replicate)
    complex = petls.Alpha(torus, max_dim=3)
    return complex

def get_filtrations(delta):
    start = 0.0
    end = 5.0
    filtrations = [start+i*delta for i in range(int((end-start)/delta) + 1)]
    return filtrations

def output_torus(replicate):

    delta = 0.1
    complex = get_complex(replicate)
    filtrations = get_filtrations(delta)
    dims = [0, 1, 2]

    for filtration in filtrations:
        for dim in dims:
            # C++ code will just read in the boundary matrix from a file            
            complex.store_L(dim, a=filtration, b=filtration+delta, prefix=f"{folder}/dim{dim}_a{filtration:.2f}_b{filtration+delta:.2f}_r{replicate}.mkt")

# for replicate in range(100):
#     output_torus(replicate=replicate)

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        replicate = int(sys.argv[1])
        output_torus(replicate=replicate)
    else:
        print("Please provide a replicate number as an argument.")