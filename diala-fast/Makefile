default: run fes.dat

run:
	python diala.py

fes.dat: HILLS
	plumed sum_hills --hills HILLS

fes.png: fes.dat
	gnuplot fes.gp


clean:
	rm bck.* HILLS colvar fes.dat fes.png plumed.log output.coor output.vel output.dcd *.xml
