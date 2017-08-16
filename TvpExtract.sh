#!/bin/bash
if test -e log.lammps 
	then 
                echo "Start scraping for relevent data."

		grep "Temp" log.lammps | tr "T" "\n" | grep emp | replace emp Temp | replace "=" "" > tvp.csv
		echo "Finish scraping Temperature"

		grep Volume log.lammps | tr "=" " " >> tvp.csv
		echo "Finish scraping volume."

		grep "Press"  log.lammps | tr "P" "\n" | grep ress | replace ress Press | replace "=" "" >> tvp.csv
                echo "Finish scraping pressure." 
	else 
		echo "The file log.lammps is not here."
fi

