# *Che*Sweet

[![Build Status](https://travis-ci.org/BIOS-IMASL/chesweet.svg?branch=master)](https://travis-ci.org/BIOS-IMASL/chesweet)

*Che*Sweet: chemical shifts for glycans  

The knowledge of the tridimensional structure of glycans is necessary to understand in detail, at atomic level, the molecular processes in which they are involved. Chemical Shifts (CS) are observables obtained from Nuclear Magnetic Resonance experiments that are highly sensitive probes to sense conformational changes. CS can be calculated accurately using quantum chemical tools, although these computations are very demanding for routine computations of more than a few conformations. For that reason we have developed *Che*Sweet, a Python module for accurate and fast computation of CS.

## Installing

### Requirements

*Che*Sweet requires:
- Python (>=3.3)
- Numpy (>= 1.8.2)
- Scipy (>= 0.13.3)

This packages will be automatically installed if you don't have them yet.

### Installation

#### Version 1

1. Download this repository and unzip it. 

2. Execute the following command in your terminal:
 
```
pip install chesweet-master/
```

#### Version 2

1. Clone this repository

2. Execute the following command in your terminal:

```
pip install chesweet/
```

## Testing

After installation is done you can run the automatic tests (optional), you will need to have `pytest` package installed, depending on the way that you downloaded the folder run:

```
pytest chesweet-master
```
**or**

```
pytest chesweet
```

This will check that all functions of *Che*Sweet work in your installation. You can also use the `nose` package instead of `pytest`.

## Using *Che*Sweet

### Import module
```python
import chesweet as chsw
```

### Load look-up tables (lut)
One of the main argument of `load` is `full`. Passing `False` is useful when working with the glycosidic torsional angles $\phi$, $\psi$ and $\omega$, `True` allows you to work with the mentioned torsional angles and also included the $\chi$ torsionals. In the image below you can see the mentioned torsional angles for $\alpha$-D-Glcp-(1-4)-$\alpha$-D-Glcp.  

<img src="img/maltose_all_torsional.png" width=500>

The names of the disaccharides used by *Che*Sweet follow the IUPAC syntax, but without the parentheses in the bond, e.g. for maltose the IUPAC name is a-D-Glcp-(1-4)-a-D-Glcp, *Che*Sweet uses a-D-Glcp-1-4-a-D-Glcp.

Examples:  
Load **all calculated dissacharides** in *Che*Sweet:  
```python
# Reduced version of the datasets (using only phi$, psi and omega)
disaccharides_red = chsw.CheSweet()

# Using full=True (using all torsionals)
disaccharides_full = chsw.CheSweet(full=True)
```

Load **one dissacharide**:  
```python
# Reduced version of the datasets
maltose_red = chsw.CheSweet(disaccharides = ['a-D-Glcp-1-4-a-D-Glcp'])

# Using full=True
maltose_full = chsw.CheSweet(disaccharides = ['a-D-Glcp-1-4-a-D-Glcp'], full=True)
```

Load **a number of dissacharides**:  
```python
disaccharides_list = ['a-D-Glcp-1-4-a-D-Glcp', 'a-D-Glcp-1-1-a-D-Glcp', 'a-D-Galp-1-3-b-D-Galp', 'b-D-Galp-1-6-b-D-Galp']

# Reduced version of the datasets
disaccharides_red = chsw.CheSweet(disaccharides = disaccharides_list)

# Using full=True
disaccharides_full = CheSweet(disaccharides = disaccharides_list, full=True)
```

Once you have loaded the lool-up tables it is possible to calculate the chemical shifts or the torsional that you want.  
In the next examples we show how to use these functions for maltose \[$\alpha$-D-Glcp-(1-4)-$\alpha$-D-Glcp]

### Calculate chemical shifts from torsional values (using `compute_cs()` function)

Using the reduced look-up table you have to provide the name of the disaccharide and the values of $\phi$ and $\psi$ torsionals:  

```python
chemC1_red, chemCx_red = maltose_red.compute_cs('a-D-Glcp-1-4-a-D-Glcp', 85.3, 76.8)
print(chemC1_red, chemCx_red)
102.852047 76.693767
```

When using the full look-up table you have to also provide the $\chi$ torsional angles:

```python
# Using full=True
chemC1_full, chemCx_full = maltose_full.compute_cs('a-D-Glcp-1-4-a-D-Glcp', 85.3, 76.8, 46.8, 78.6, 173.1)
print(chemC1_full, chemCx_full)
101.85525 72.19692
```

You can see in these examples that depending on the datasets used (reduced or full) the obtained result can have little differences (Garay et. al 2014).

### Obtain torsional list from CS values (using `compute_tors()` function)

This example is the inverse of the previous one. We are now passing the CS for the carbons involved in the glycosidic bond and we are getting the compatible torsional angles given a tolerance `eps`, by default `eps=0.5`. The first CS should be C1 and the second CS should be the second carbon in the glycosidic bond.

```python
# Reduced version of the datasets
tors_array_red = maltose_red.compute_tors('a-D-Glcp-1-4-a-D-Glcp', 109.52, 87.68)
print(tors_array_red)
[[  90.  120.]
 [ 110.  100.]
 [ 130.  130.]
 [ 140.  100.]]
```
In this case, we obtain information only about $\phi$ and $\psi$, first and second column respectively.
You can obtain more values for these torsionals by increasing the value of the parameter `eps` (by default 0.5):  

```python
# augmented eps
tors_array_red2 = maltose_red.compute_tors('a-D-Glcp-1-4-a-D-Glcp', 109.52, 87.68, eps=0.6)
print(tors_array_red2)
[[  90.  120.]
 [ 100.  110.]
 [ 110.  100.]
 [ 130.   90.]
 [ 130.  130.]
 [ 140.  100.]]
```

Alternatively, you can use the full datasets.  

```python
# Using full=True
tors_array_full = maltose_full.compute_tors('a-D-Glcp-1-4-a-D-Glcp', 109.52, 87.68)
print(tors_array_full)
[[  80.  120.  -60.  180.  -60.]
 [  90.  110.  180.  180.  -60.]
 [  90.  110.  -60.   60.   60.]
 [  90.  120.  -60.  -60.  -60.]
 [  90.  130.  180.   60.  -60.]
 [ 100.  100.  180.  180.  -60.]
 [ 100.  130.   60.  180.  180.]
 [ 110.   90.  180.  180.  -60.]
 [ 110.  100.  -60.  -60.  -60.]
 [ 110.  130.   60.  180.  180.]
 [ 110.  160.   60.   60.   60.]
 [ 110.  160.  180.   60.   60.]
 [ 120.  100.   60.   60.  -60.]
 [ 130.  100.   60.   60.  -60.]
 [ 130.  140.  180.   60.   60.]
 [ 140.  130.  -60.  -60.  -60.]]
```

You can see in these examples that depending on the datasets used (reduced or full) the obtained result can have little differences and be using the full version of the dataset you have additional information about the $\chi$ angles.

## References

Garay, P. G., Martin, O. A., Scheraga, H. A., & Vila, J. A. (2014). Factors affecting the computation of the $^{13}$C shielding in disaccharides. Journal of Computational Chemistry, 35(25), 1854–1864. https://doi.org/10.1002/jcc.23697

## Contributing

All contributions are welcome!  
You can contribute via GitHub, send your proposed changes with your Pull Requests and the errors/improvements/comments must be reported by Issues.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
