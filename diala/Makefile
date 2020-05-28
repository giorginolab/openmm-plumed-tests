default: run fes.dat

run:
	acemd acemd_input

fes.dat:
	plumed sum_hills --hills HILLS

fes.png: fes.dat
	gnuplot fes.gp


clean:
	rm bck.* HILLS colvar fes.dat plumed.log output.coor output.vel output.dcd
