import gudhi
from petls import Complex

def test_simplex_tree():
    st = gudhi.SimplexTree()
    st.insert([0, 1])
    st.insert([0, 1, 2], filtration=4.0)
    pl = Complex(simplex_tree=st)
    print(pl.spectra())

if __name__=="__main__":
    test_simplex_tree()