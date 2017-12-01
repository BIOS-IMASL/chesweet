---
title: 'CheSweet: An application to predict glycan’s chemicals shifts'
tags:
  - glycans
  - chemical shifts
  - python
  - carbohydrates
authors:
 - name: Pablo G Garay
   orcid: 0000-0001-6604-0329
   affiliation: 1
 - name: Jorge A Vila
   affiliation: 1
 - name: Osvaldo A Martin
   orcid: 0000-0001-7419-8978
   affiliation: 1
affiliations:
 - name: Instituto de Matemática San Luis - CONICET, Universidad Nacional de San Luis, San Luis, Argentina
   index: 1
date: 1 December 2017
bibliography: paper.bib
---

# Summary

Glycans are the most abundant and structurally diverse biomolecules in nature. The knowledge of the tridimensional structure of these molecules is necessary to understand in detail, at atomic level, the molecular processes in which they are involved. Chemical shifts (CS) are observables obtained from Nuclear Magnetic Resonance experiments and are highly sensitive probes to sense conformational changes (@Swalina2001, @Martin2012, @Garay2016). Here we present CheSweet a Python module to compute CS for glycans. The core of CheSweet is the fast calculation of CS, based on the pre-calculated values of CS at DFT level of theory, from the values of the torsional angles ($\phi$, $\psi$, $\omega$ and $\chi$, depending on its availability) of the glycosidic bond (@Garay2014). This calculation of CS is done through a linear interpolation (based on the existence of 4 points, two for $\phi$ and two for $\psi$), or by nearest neighbor interpolation if there is less information available. Is possible to have less than 4 values because the lookup table, from were CheSweet do the interpolation, have been constructed from low-energy conformations of disaccharides, thus not all $\phi$-$\psi$ pairs exist. If the number of dihedral angles is less than 4 (3, 2 or 1) is not possible to do a linear interpolation, in such a case, CheSweet returns the CS value of the nearest neighbor. If all the torsional values in the lookup table are at a distance larger than 10$^{\circ}$, from the provided torsionals, it is considered that there are no neighbors, then CheSweet returns \texttt{INF} (positive infinite). In a similar fashion CheSweet can be used to solve the inverse problem, i.e. to compute a set of torsional angles compatible with the provided CS values, from the carbons of the glycosidic bond. CheSweet has the potential to be used as part of more complex methods to predict, validate and refine glycan’s structures.


# References

